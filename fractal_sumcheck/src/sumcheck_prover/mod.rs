use std::{convert::TryInto, marker::PhantomData};

use fractal_utils::polynomial_utils::*;
use winter_crypto::ElementHasher;
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_math::{fft, FieldElement, StarkField};

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
        // let mut summing_poly_numerator_evals = self.summing_poly_numerator.clone();
        // let mut eval_domain_twiddles = fft::get_twiddles(self.summing_domain.len());

        // println!("summing_poly_evals len = {:?}", summing_poly_numerator_evals.len());
        // // let size_num_evals = summing_poly_numerator_evals.len().next_power_of_two();
        // let size_num_evals = self.summing_domain.len();
        // println!("Numerator evals = {}", size_num_evals);
        // pad_with_zeroes(&mut summing_poly_numerator_evals, size_num_evals * 2);
        
        // println!("Numerator evals = {}", summing_poly_numerator_evals.len());
        // println!("Num twiddles = {}", eval_domain_twiddles.len());
        
        // fft::evaluate_poly(
        //     &mut summing_poly_numerator_evals,
        //     &mut eval_domain_twiddles,
        // );
        
        // let mut summing_poly_denominator_evals = self.summing_poly_denominator.clone();
        // // let size_denom_evals = summing_poly_denominator_evals.len().next_power_of_two();
        // let size_denom_evals = self.evaluation_domain.len();
        // println!("Denominator evals = {}", summing_poly_denominator_evals.len());
        // pad_with_zeroes(&mut summing_poly_denominator_evals, size_denom_evals);
        // println!("Denominator evals = {}", summing_poly_denominator_evals.len());
        // fft::evaluate_poly(
        //     &mut summing_poly_denominator_evals,
        //     &mut eval_domain_twiddles,
        // );
        // println!("Denominator evals = {}", summing_poly_denominator_evals.len());

        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        // compute the polynomial e such that e = (Sigma(g, sigma) - summing_poly)/v_H over the summing domain H.
        println!("Starting a sumcheck proof");
        let mut g_eval_domain_evals: Vec<E> = Vec::new();
        let mut e_eval_domain_evals: Vec<E> = Vec::new();
        let mut f_hat_evals: Vec<E> = Vec::new();
        let _sigma_inv = self.sigma.inv();
        /*for i in 0..self.summing_poly_numerator_evals.len() {
            let summing_poly_eval = B::div(
                self.summing_poly_numerator_evals[i],
                self.summing_poly_denominator_evals[i],
            );
            f_hat_evals.push(E::from(summing_poly_eval));
            /*let g_val = self
                .compute_g_poly_on_val(E::from(self.evaluation_domain[i]), E::from(summing_poly_eval));
            g_eval_domain_evals.push(g_val);
            let e_val = self.compute_e_poly_on_val(
                E::from(self.evaluation_domain[i]),
                g_val,
                E::from(self.summing_poly_numerator_evals[i]),
                E::from(self.summing_poly_denominator_evals[i]),
            );
            e_eval_domain_evals.push(e_val);*/
        }*/

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
        //let g_coeffs = polynom::interpolate(&eval_domain_e, &g_eval_domain_evals, true);
        println!("self.evaluation_domain.len(): {:?}", &self.evaluation_domain.len());
        //println!("degree of g_coeffs {}", polynom::degree_of(&g_coeffs));
        //let summing_poly_coeffs = polynom::interpolate(&eval_domain_e, &summing_poly_evals, true);
        //println!("degree of summing_poly_coeffs {}", polynom::degree_of(&summing_poly_coeffs));
        //let g_eval_domain_evals2: Vec<E> = polynom::eval_many(g_coeffs.clone().as_slice(), eval_domain_e.clone().as_slice());// Vec::new();
        //println!("old evals: {:?}, new evals: {:?}", &g_eval_domain_evals, &g_eval_domain_evals2);


        ////g_hat test
        let dividing_factor_for_sigma: u64 = self.summing_domain.len().try_into().unwrap();
        let subtracting_factor = self.sigma * B::from(dividing_factor_for_sigma).inv();
        println!("sigma: {}", &self.sigma);
        println!("subtracting factor: {}", &subtracting_factor);
        //println!("f_hat(x): {:?}", &summing_poly_coeffs);
        /// 


        //let e_coeffs = polynom::interpolate(&eval_domain_e, &e_eval_domain_evals, true);
        //println!("degree of e_coeffs {}", polynom::degree_of(&e_coeffs));
        //println!("e_eval_domain_evals {:?}", e_eval_domain_evals); //all 0's

        //let g_comp_coeffs = get_complementary_poly::<E>(polynom::degree_of(&g_coeffs), 64);//self.max_degree - 1);
        //let new_g = polynom::mul(&g_coeffs, &g_comp_coeffs);
        //let g_evals = polynom::eval_many(&new_g, &eval_domain_e);
        //g_eval_domain_evals = g_evals;

        let g_eval_domain_evals = polynom::eval_many(&g_hat_coeffs, &eval_domain_e);

        let p_eval_domain_evals = polynom::eval_many(&self.numerator_coeffs, &eval_domain_e);
        let q_eval_domain_evals = polynom::eval_many(&self.denominator_coeffs, &eval_domain_e);

        let mut e_eval_domain_evals: Vec<E> = Vec::new();
        for i in 0..self.evaluation_domain.len() {
            let e_val = self.compute_e_poly_on_val(
                E::from(self.evaluation_domain[i]),
                g_eval_domain_evals[i],
                p_eval_domain_evals[i],
                q_eval_domain_evals[i],
                self.eta,
            );
            e_eval_domain_evals.push(e_val);
        }

        println!("degree of e: {}", polynom::degree_of(&polynom::interpolate(&eval_domain_e, &e_eval_domain_evals, true)));
        
        //let inv_twiddles_eval_domain: Vec<B> = fft::get_inv_twiddles(self.evaluation_domain.len());
        //let mut g_poly = g_eval_domain_evals.clone(); //g_summing_domain_evals.clone();
        //let mut e_poly = g_eval_domain_evals.clone(); //e_summing_domain_evals.clone();
        //fft::interpolate_poly(&mut g_poly, &inv_twiddles_eval_domain);
        //fft::interpolate_poly(&mut e_poly, &inv_twiddles_eval_domain);
        //println!("g_len = {}", g_poly.len());
        //println!("e_len = {}", e_poly.len());
        //print!("eval_domain_len = {}", self.evaluation_domain.len());

        // let twiddles_evaluation_domain: Vec<B> = fft::get_twiddles(self.evaluation_domain.len());
        // let mut g_eval_domain_evals = g_poly.clone();
        // let mut e_eval_domain_evals = e_poly.clone();
        // fft::evaluate_poly(&mut g_eval_domain_evals, &twiddles_evaluation_domain);
        // fft::evaluate_poly(&mut e_eval_domain_evals, &twiddles_evaluation_domain);
        // let mut channel = DefaultProverChannel::new(self.evaluation_domain.len(), self.num_queries);
        let query_positions = self.channel.draw_query_positions();
        let queried_positions = query_positions.clone();

        // Build proofs for the polynomial g
        let mut fri_prover =
            winter_fri::FriProver::<B, E, DefaultProverChannel<B, E, H>, H>::new(self.fri_options.clone());
        fri_prover.build_layers(&mut self.channel, g_eval_domain_evals.clone());
        let fri_proof_g = fri_prover.build_proof(&query_positions);
        let g_queried_evaluations = query_positions.clone()
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
            .map(|&p| e_eval_domain_evals[p])
            .collect::<Vec<_>>();
        //todo: consider being less hacky
        let e_commitments = self.channel.layer_commitments()[self.channel.layer_commitments().len()/2..].to_vec();
        println!("@@@@@@@@@@@@Prover's queried positions {:?} ", &queried_positions);

        SumcheckProof {
            options: self.fri_options.clone(),
            num_evaluations: self.evaluation_domain.len(),
            queried_positions,
            g_proof: fri_proof_g,
            g_queried: OracleQueries::new(g_queried_evaluations, vec![g_commitments]),
            g_max_degree: self.g_degree,
            e_proof: fri_proof_e,
            e_queried: OracleQueries::new(e_queried_evaluations, vec![e_commitments]),
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

    pub fn compute_sigma_function_on_val(&self, x_val: E, g_val: E) -> E {
        let dividing_factor: u64 = self.summing_domain.len().try_into().unwrap();
        x_val * g_val + (E::from(self.sigma) * E::from(dividing_factor).inv())
    }

    pub fn compute_e_poly_on_val(
        &self,
        x_val: E,
        g_val: E,
        summing_poly_numerator_val: E,
        summing_poly_denominator_val: E,
        eta: B,
    ) -> E {
        let sigma_function = self.compute_sigma_function_on_val(x_val, g_val);
        let sigma_minus_f =
            sigma_function * summing_poly_denominator_val - summing_poly_numerator_val;
        let vanishing_on_x = compute_vanishing_poly(x_val, E::from(eta), self.summing_domain.len());
            //vanishing_poly_for_mult_subgroup(x_val, self.summing_domain.len().try_into().unwrap());
        sigma_minus_f * vanishing_on_x.inv()
    }
}
