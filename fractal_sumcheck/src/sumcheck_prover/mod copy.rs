use std::{convert::TryInto, marker::PhantomData};

use crypto::ElementHasher;
use fractal_utils::polynomial_utils::*;
use fri::{DefaultProverChannel, FriOptions};
use math::{fft, FieldElement, StarkField};

use fractal_proofs::{OracleQueries, SumcheckProof};
#[cfg(test)]
mod tests;

pub struct RationalSumcheckProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
> {
    // this is p(x) as in the API for sumcheck in the paper
    summing_poly_numerator: Vec<E::BaseField>,
    // this is q(x)
    summing_poly_denominator: Vec<E::BaseField>,
    // this is \sigma as in the sumcheck i.e. the desired sum
    sigma: E,
    // For lincheck this domain is K
    summing_domain: Vec<E::BaseField>,
    summing_domain_twiddles: Vec<B>,
    // Eval domain is always L
    evaluation_domain: Vec<E::BaseField>,
    fri_options: FriOptions,
    pub channel: DefaultProverChannel<B, E, H>,
    _h: PhantomData<H>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    RationalSumcheckProver<B, E, H>
{
    pub fn new(
        summing_poly_numerator: Vec<B>,
        summing_poly_denominator: Vec<B>,
        sigma: E,
        summing_domain: Vec<B>,
        evaluation_domain: Vec<B>,
        fri_options: FriOptions,
        num_queries: usize,
    ) -> Self {
        let summing_domain_twiddles = fft::get_twiddles(summing_domain.len());
        let channel = DefaultProverChannel::new(evaluation_domain.len(), num_queries);
        RationalSumcheckProver {
            summing_poly_numerator,
            summing_poly_denominator,
            sigma,
            summing_domain,
            summing_domain_twiddles,
            evaluation_domain,
            fri_options,
            channel,
            _h: PhantomData,
        }
    }

    pub fn generate_proof(&mut self) -> SumcheckProof<B, E, H> {
        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        let mut summing_poly_numerator_evals = self.summing_poly_numerator.clone();
        let mut eval_domain_twiddles = fft::get_twiddles(self.summing_domain.len());

        println!("summing_poly_evals len = {:?}", summing_poly_numerator_evals.len());
        // let size_num_evals = summing_poly_numerator_evals.len().next_power_of_two();
        let size_num_evals = self.summing_domain.len();
        println!("Numerator evals = {}", size_num_evals);
        pad_with_zeroes(&mut summing_poly_numerator_evals, size_num_evals * 2);
        
        println!("Numerator evals = {}", summing_poly_numerator_evals.len());
        println!("Num twiddles = {}", eval_domain_twiddles.len());
        
        fft::evaluate_poly(
            &mut summing_poly_numerator_evals,
            &mut eval_domain_twiddles,
        );
        
        let mut summing_poly_denominator_evals = self.summing_poly_denominator.clone();
        // let size_denom_evals = summing_poly_denominator_evals.len().next_power_of_two();
        let size_denom_evals = self.evaluation_domain.len();
        println!("Denominator evals = {}", summing_poly_denominator_evals.len());
        pad_with_zeroes(&mut summing_poly_denominator_evals, size_denom_evals);
        println!("Denominator evals = {}", summing_poly_denominator_evals.len());
        fft::evaluate_poly(
            &mut summing_poly_denominator_evals,
            &mut eval_domain_twiddles,
        );
        println!("Denominator evals = {}", summing_poly_denominator_evals.len());

        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        // compute the polynomial e such that e = (Sigma(g, sigma) - summing_poly)/v_H over the summing domain H.
        let mut g_summing_domain_evals: Vec<E> = Vec::new();
        let mut e_summing_domain_evals: Vec<E> = Vec::new();
        let _sigma_inv = self.sigma.inv();
        for i in 0..summing_poly_numerator_evals.len() {
            let summing_poly_eval = B::div(
                summing_poly_numerator_evals[i],
                summing_poly_denominator_evals[i],
            );
            let g_val = self
                .compute_g_poly_on_val(E::from(self.evaluation_domain[i]), E::from(summing_poly_eval));
            g_summing_domain_evals.push(g_val);
            let e_val = self.compute_e_poly_on_val(
                E::from(self.evaluation_domain[i]),
                g_val,
                E::from(summing_poly_numerator_evals[i]),
                E::from(summing_poly_denominator_evals[i]),
            );
            e_summing_domain_evals.push(e_val);
        }
        let inv_twiddles_eval_domain: Vec<B> = fft::get_inv_twiddles(self.evaluation_domain.len());
        let mut g_poly = g_summing_domain_evals.clone();
        let mut e_poly = e_summing_domain_evals.clone();
        fft::interpolate_poly(&mut g_poly, &inv_twiddles_eval_domain);
        fft::interpolate_poly(&mut e_poly, &inv_twiddles_eval_domain);

        let twiddles_evaluation_domain: Vec<B> = fft::get_twiddles(self.evaluation_domain.len());
        let mut g_eval_domain_evals = g_poly.clone();
        let mut e_eval_domain_evals = e_poly.clone();
        fft::evaluate_poly(&mut g_eval_domain_evals, &twiddles_evaluation_domain);
        fft::evaluate_poly(&mut e_eval_domain_evals, &twiddles_evaluation_domain);
        // let mut channel = DefaultProverChannel::new(self.evaluation_domain.len(), self.num_queries);
        let query_positions = self.channel.draw_query_positions();
        let queried_positions = query_positions.clone();

        // Build proofs for the polynomial g
        let mut fri_prover =
            fri::FriProver::<B, E, DefaultProverChannel<B, E, H>, H>::new(self.fri_options.clone());
        fri_prover.build_layers(&mut self.channel, g_eval_domain_evals.clone());
        let fri_proof_g = fri_prover.build_proof(&query_positions);
        let g_queried_evaluations = query_positions
            .iter()
            .map(|&p| g_eval_domain_evals[p])
            .collect::<Vec<_>>();
        let g_commitments = self.channel.layer_commitments().to_vec();

        // reset to build proofs for the polynomial e
        fri_prover.reset();
        fri_prover.build_layers(&mut self.channel, e_eval_domain_evals.clone());
        let fri_proof_e = fri_prover.build_proof(&query_positions);
        let e_queried_evaluations = query_positions
            .iter()
            .map(|&p| e_summing_domain_evals[p])
            .collect::<Vec<_>>();
        let e_commitments = self.channel.layer_commitments().to_vec();

        SumcheckProof {
            options: self.fri_options.clone(),
            num_evaluations: self.evaluation_domain.len(),
            queried_positions,
            g_proof: fri_proof_g,
            g_queried: OracleQueries::new(g_queried_evaluations, vec![g_commitments]),
            g_max_degree: self.summing_poly_numerator.len() - self.summing_poly_denominator.len(),
            e_proof: fri_proof_e,
            e_queried: OracleQueries::new(e_queried_evaluations, vec![e_commitments]),
            e_max_degree: self.summing_poly_numerator.len()
                - self.summing_poly_denominator.len()
                - self.summing_domain.len()
                + 1,
        }
    }

    // SIGMA(g, sigma)(x) = f(x) = p(x)/q(x)
    // SIGMA(g, sigma) = x*g(x) + sigma*|summing_domain|^-1
    // g(x) = x^-1*(f(x) - sigma*|summing_domain|^-1)
    pub fn compute_g_poly_on_val(&self, x_val: E, f_x_val: E) -> E {
        let dividing_factor_for_sigma: u64 = self.summing_domain.len().try_into().unwrap();
        let subtracting_factor = self.sigma * E::from(dividing_factor_for_sigma).inv();
        let dividing_factor = x_val.inv();
        dividing_factor * (f_x_val - subtracting_factor)
    }

    pub fn compute_sigma_function_on_val(&self, x_val: E, g_val: E) -> E {
        let dividing_factor: u64 = self.summing_domain.len().try_into().unwrap();
        x_val * g_val + (self.sigma * E::from(dividing_factor).inv())
    }

    pub fn compute_e_poly_on_val(
        &self,
        x_val: E,
        g_val: E,
        summing_poly_numerator_val: E,
        summing_poly_denominator_val: E,
    ) -> E {
        let sigma_function = self.compute_sigma_function_on_val(x_val, g_val);
        let sigma_minus_f =
            sigma_function * summing_poly_denominator_val - summing_poly_numerator_val;
        let vanishing_on_x =
            vanishing_poly_for_mult_subgroup(x_val, self.summing_domain.len().try_into().unwrap());
        sigma_minus_f * vanishing_on_x.inv()
    }
}
