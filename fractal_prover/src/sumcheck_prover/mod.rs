use std::{convert::TryInto, marker::PhantomData};

use crate::accumulator::Accumulator;
use crate::errors::ProverError;
use crate::low_degree_batch_prover::LowDegreeBatchProver;
use crate::low_degree_prover::LowDegreeProver;
use crate::{channel::DefaultFractalProverChannel, log::debug};
use crate::{FractalOptions, LayeredProver};
use fractal_utils::polynomial_utils::*;
use winter_crypto::ElementHasher;
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_math::{fft, FieldElement, StarkField};

use fractal_proofs::{polynom, OracleQueries, SumcheckProof};
#[cfg(test)]
mod tests;

pub struct RationalSumcheckProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
> {
    // this is p(x) as in the API for sumcheck in the paper
    numerator_coeffs: Vec<E>,
    // this is q(x)
    denominator_coeffs: Vec<E>,
    // this is \sigma as in the sumcheck i.e. the desired sum
    sigma: E,
    // For lincheck this domain is K
    summing_domain: Vec<E::BaseField>,
    eta: B,
    #[allow(dead_code)]
    summing_domain_twiddles: Vec<B>,
    // Eval domain is always L
    g_degree: usize,
    e_degree: usize,
    fractal_options: FractalOptions<B>,
    _h: PhantomData<H>,
    current_layer: usize,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    RationalSumcheckProver<B, E, H>
{
    pub fn new(
        numerator_coeffs: Vec<E>,
        denominator_coeffs: Vec<E>,
        sigma: E,
        summing_domain: Vec<B>,
        eta: B,
        g_degree: usize,
        e_degree: usize,
        fractal_options: FractalOptions<B>,
    ) -> Self {
        let summing_domain_twiddles = fft::get_twiddles(summing_domain.len());

        RationalSumcheckProver {
            numerator_coeffs,
            denominator_coeffs,
            sigma,
            summing_domain,
            eta,
            summing_domain_twiddles,
            g_degree,
            e_degree,
            fractal_options,
            _h: PhantomData,
            current_layer: 0,
        }
    }

    pub fn sumcheck_layer_one(
        &self,
        //query: E,
        accumulator: &mut Accumulator<B, E, H>,
        //channel: &mut DefaultFractalProverChannel<B, E, H>,
        //initial_queries: Vec<usize>,
    ) {
        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        // compute the polynomial e such that e = (Sigma(g, sigma) - summing_poly)/v_H over the summing domain H.
        debug!("Starting a sumcheck proof");
        let _sigma_inv = self.sigma.inv();

        //might be faster to eval_many
        let f_hat_evals: Vec<E> = self
            .summing_domain
            .iter()
            .map(|x| {
                polynom::eval(&self.numerator_coeffs, E::from(*x))
                    / polynom::eval(&self.denominator_coeffs, E::from(*x))
            })
            .collect();

        let summing_domain_e: Vec<E> = self.summing_domain.iter().map(|f| E::from(*f)).collect();
        let f_hat_coeffs = polynom::interpolate(
            &self
                .summing_domain
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
            &f_hat_evals,
            true,
        );
        let x_coeffs = vec![E::ZERO, E::ONE];
        let sub_factor = self.sigma / E::from(self.summing_domain.len() as u64);
        let f_hat_minus_sub_factor = polynom::sub(&f_hat_coeffs, &vec![E::from(sub_factor)]);
        assert_eq!(f_hat_minus_sub_factor[0], E::ZERO);
        let g_hat_coeffs = polynom::div(&f_hat_minus_sub_factor, &x_coeffs);

        let eval_domain_e: Vec<E> = self
            .fractal_options
            .evaluation_domain
            .iter()
            .map(|f| E::from(*f))
            .collect();

        debug!(
            "self.evaluation_domain.len(): {:?}",
            &self.fractal_options.evaluation_domain.len()
        );

        ////g_hat test
        let dividing_factor_for_sigma: u64 = self.summing_domain.len().try_into().unwrap();
        let subtracting_factor = self.sigma * E::from(dividing_factor_for_sigma).inv();

        let g_eval_domain_evals = polynom::eval_many(
            &g_hat_coeffs,
            &self
                .fractal_options
                .evaluation_domain
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
        );

        let p_eval_domain_evals = polynom::eval_many(
            &self.numerator_coeffs,
            &self
                .fractal_options
                .evaluation_domain
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
        );
        let q_eval_domain_evals = polynom::eval_many(
            &self.denominator_coeffs,
            &self
                .fractal_options
                .evaluation_domain
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
        );

        let mut g_summing_domain_evals: Vec<E> = Vec::new();
        for i in 0..self.fractal_options.summing_domain.len() {
            let g_val = self.compute_g_poly_on_val(
                E::from(self.fractal_options.summing_domain[i]),
                E::from(f_hat_evals[i]),
            );
            g_summing_domain_evals.push(g_val);
        }
        
        let mut e_eval_domain_evals: Vec<E> = Vec::new();
        for i in 0..self.fractal_options.evaluation_domain.len() {
            let e_val = self.compute_e_poly_on_val(
                E::from(self.fractal_options.evaluation_domain[i]),
                E::from(g_eval_domain_evals[i]),
                p_eval_domain_evals[i],
                q_eval_domain_evals[i],
                self.eta,
            );
            e_eval_domain_evals.push(e_val);
        }
        // let mut e_eval_domain_evals: Vec<E> = Vec::new();
        // for i in 0..self.fractal_options.evaluation_domain.len() {
        //     let e_val = self.compute_e_poly_on_val(
        //         E::from(self.fractal_options.evaluation_domain[i]),
        //         E::from(g_eval_domain_evals[i]),
        //         p_eval_domain_evals[i],
        //         q_eval_domain_evals[i],
        //         self.eta,
        //     );
        //     e_eval_domain_evals.push(e_val);
        // }

        let e_hat_coeffs = polynom::interpolate(
            &self
                .fractal_options
                .evaluation_domain
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
            &e_eval_domain_evals,
            true,
        );
        debug!("degree of e: {}", polynom::degree_of(&e_hat_coeffs));

        accumulator.add_polynomial_e(g_hat_coeffs, self.g_degree);
        accumulator.add_polynomial_e(e_hat_coeffs, self.e_degree);
    }

    // SIGMA(g, sigma)(x) = f(x) = p(x)/q(x)
    // SIGMA(g, sigma) = x*g(x) + sigma*|summing_domain|^-1
    // g(x) = x^-1*(f(x) - sigma*|summing_domain|^-1)
    pub fn compute_g_poly_on_val(&self, x_val: E, f_x_val: E) -> E {
        let dividing_factor_for_sigma: u64 = self.summing_domain.len().try_into().unwrap();
        let subtracting_factor = self.sigma * E::from(dividing_factor_for_sigma).inv();
        let dividing_factor = x_val.inv();
        dividing_factor * (f_x_val - E::from(subtracting_factor))
    }

    pub fn compute_sigma_function_on_val(&self, x_val: E, g_val: E) -> E {
        let dividing_factor: u64 = self.summing_domain.len().try_into().unwrap();
        (x_val * g_val) + (E::from(self.sigma) * E::from(dividing_factor).inv())
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
        sigma_minus_f * vanishing_on_x.inv()
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredProver<B, E, H> for RationalSumcheckProver<B, E, H>
{
    fn run_next_layer(
        &mut self,
        _query: E,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError> {
        if self.get_current_layer() == 0 {
            self.sumcheck_layer_one(accumulator);
            self.current_layer += 1;
        }
        Ok(())
    }
    fn get_current_layer(&self) -> usize {
        self.current_layer
    }
    fn get_num_layers(&self) -> usize {
        1
    }
    fn get_fractal_options(&self) -> crate::FractalOptions<B> {
        self.fractal_options.clone()
    }
}
