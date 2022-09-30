use std::{convert::TryInto, marker::PhantomData};

use fractal_utils::polynomial_utils::*;
use low_degree::low_degree_prover::LowDegreeProver;
use winter_crypto::ElementHasher;
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_math::{fft, FieldElement, StarkField};
use crate::log::debug;

use fractal_proofs::{OracleQueries, SumcheckProof, polynom};
#[cfg(test)]
mod tests;

pub struct RationalSumcheckProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
> {
    // this is p(x) as in the API for sumcheck in the paper
    numerator_coeffs: Vec<B>,
    // this is q(x)
    denominator_coeffs: Vec<B>,
    // this is \sigma as in the sumcheck i.e. the desired sum
    sigma: B,
    // For lincheck this domain is K
    summing_domain: Vec<E::BaseField>,
    eta: B,
    #[allow(dead_code)]
    summing_domain_twiddles: Vec<B>,
    // Eval domain is always L
    evaluation_domain: Vec<E::BaseField>,
    g_degree: usize,
    e_degree: usize,
    fri_options: FriOptions,
    pub channel: DefaultProverChannel<B, E, H>,
    _h: PhantomData<H>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    RationalSumcheckProver<B, E, H>
{
    pub fn new(
        numerator_coeffs: Vec<B>,
        denominator_coeffs: Vec<B>,
        sigma: B,
        summing_domain: Vec<B>,
        eta: B,
        evaluation_domain: Vec<B>,
        g_degree: usize,
        e_degree: usize,
        fri_options: FriOptions,
        num_queries: usize,
    ) -> Self {
        let summing_domain_twiddles = fft::get_twiddles(summing_domain.len());
        let channel = DefaultProverChannel::new(evaluation_domain.len(), num_queries);
        RationalSumcheckProver {
            numerator_coeffs,
            denominator_coeffs,
            sigma,
            summing_domain,
            eta,
            summing_domain_twiddles,
            evaluation_domain,
            g_degree,
            e_degree,
            fri_options,
            channel,
            _h: PhantomData,
        }
    }

    pub fn generate_proof(&mut self) -> SumcheckProof<B, E, H> {
        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        // compute the polynomial e such that e = (Sigma(g, sigma) - summing_poly)/v_H over the summing domain H.
        debug!("Starting a sumcheck proof");
        let _sigma_inv = self.sigma.inv();
        

        //might be faster to eval_many
        let f_hat_evals: Vec<B> = self.summing_domain.iter().map(|x| polynom::eval(&self.numerator_coeffs, *x) / polynom::eval(&self.denominator_coeffs, *x)).collect();

        let summing_domain_e: Vec<E> = self.summing_domain.iter().map(|f| E::from(*f) ).collect();
        let f_hat_coeffs = polynom::interpolate(&self.summing_domain, &f_hat_evals, true);
        let x_coeffs = vec![B::ZERO, B::ONE];
        let sub_factor = self.sigma / B::from(self.summing_domain.len() as u64);
        let f_hat_minus_sub_factor = polynom::sub(&f_hat_coeffs, &vec![sub_factor]);
        assert_eq!(f_hat_minus_sub_factor[0], B::ZERO);
        let g_hat_coeffs = polynom::div(&f_hat_minus_sub_factor, &x_coeffs);
        


        let eval_domain_e: Vec<E> = self.evaluation_domain.iter().map(|f| E::from(*f) ).collect();

        debug!("self.evaluation_domain.len(): {:?}", &self.evaluation_domain.len());
        


        ////g_hat test
        let dividing_factor_for_sigma: u64 = self.summing_domain.len().try_into().unwrap();
        let subtracting_factor = self.sigma * B::from(dividing_factor_for_sigma).inv();

        let g_eval_domain_evals = polynom::eval_many(&g_hat_coeffs, &self.evaluation_domain);

        let p_eval_domain_evals = polynom::eval_many(&self.numerator_coeffs, &self.evaluation_domain);
        let q_eval_domain_evals = polynom::eval_many(&self.denominator_coeffs, &self.evaluation_domain);

        let mut e_eval_domain_evals: Vec<B> = Vec::new();
        for i in 0..self.evaluation_domain.len() {
            let e_val = self.compute_e_poly_on_val(
                self.evaluation_domain[i],
                g_eval_domain_evals[i],
                p_eval_domain_evals[i],
                q_eval_domain_evals[i],
                self.eta,
            );
            e_eval_domain_evals.push(e_val);
        }

        let e_hat_coeffs = polynom::interpolate(&self.evaluation_domain, &e_eval_domain_evals, true);
        debug!("degree of e: {}", polynom::degree_of(&e_hat_coeffs));
        
        let query_positions = self.channel.draw_query_positions();
        let queried_positions = query_positions.clone();

        // Build proofs for the polynomial g
        let g_prover = LowDegreeProver::<B, E, H>::from_polynomial(&g_hat_coeffs, &self.evaluation_domain, self.g_degree, self.fri_options.clone());
        let g_proof = g_prover.generate_proof(&mut self.channel);

        // Build proofs for the polynomial e
        let e_prover = LowDegreeProver::<B, E, H>::from_polynomial(&e_hat_coeffs, &self.evaluation_domain, self.e_degree, self.fri_options.clone());
        let e_proof = e_prover.generate_proof(&mut self.channel);

        SumcheckProof {
            options: self.fri_options.clone(),
            num_evaluations: self.evaluation_domain.len(),
            queried_positions,
            g_proof: g_proof,
            g_max_degree: self.g_degree,
            e_proof: e_proof,
            e_max_degree: self.e_degree,
        }
    }

    // SIGMA(g, sigma)(x) = f(x) = p(x)/q(x)
    // SIGMA(g, sigma) = x*g(x) + sigma*|summing_domain|^-1
    // g(x) = x^-1*(f(x) - sigma*|summing_domain|^-1)
    pub fn compute_g_poly_on_val(&self, x_val: E, f_x_val: E) -> E {
        let dividing_factor_for_sigma: u64 = self.summing_domain.len().try_into().unwrap();
        let subtracting_factor = self.sigma * B::from(dividing_factor_for_sigma).inv();
        let dividing_factor = x_val.inv();
        dividing_factor * (f_x_val - E::from(subtracting_factor))
    }

    pub fn compute_sigma_function_on_val(&self, x_val: B, g_val: B) -> B {
        let dividing_factor: u64 = self.summing_domain.len().try_into().unwrap();
        x_val * g_val + (B::from(self.sigma) * B::from(dividing_factor).inv())
    }

    pub fn compute_e_poly_on_val(
        &self,
        x_val: B,
        g_val: B,
        summing_poly_numerator_val: B,
        summing_poly_denominator_val: B,
        eta: B,
    ) -> B {
        let sigma_function = self.compute_sigma_function_on_val(x_val, g_val);
        let sigma_minus_f =
            sigma_function * summing_poly_denominator_val - summing_poly_numerator_val;
        let vanishing_on_x = compute_vanishing_poly(x_val, eta, self.summing_domain.len());
            //vanishing_poly_for_mult_subgroup(x_val, self.summing_domain.len().try_into().unwrap());
        sigma_minus_f * vanishing_on_x.inv()
    }
}
