use std::{convert::TryInto, marker::PhantomData};

use crate::errors::ProverError;
use crate::LayeredSubProver;
use fractal_accumulator::accumulator::Accumulator;
use fractal_proofs::batch_inversion;
use fractal_utils::channel::DefaultFractalProverChannel;
use fractal_utils::polynomial_utils::*;
use fractal_utils::FractalProverOptions;
use log::debug;
use low_degree_prover::low_degree_batch_prover::LowDegreeBatchProver;
use low_degree_prover::low_degree_prover::LowDegreeProver;
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
    // summing_domain: Vec<E::BaseField>,
    eta: B,
    #[allow(dead_code)]
    // summing_domain_twiddles: Vec<B>,
    // summing_domain_inv_twiddles: Vec<B>,
    // Eval domain is always L
    g_degree: usize,
    e_degree: usize,
    // fractal_options: FractalProverOptions<B>,
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
        // summing_domain: Vec<B>,
        eta: B,
        g_degree: usize,
        e_degree: usize,
        // fractal_options: FractalProverOptions<B>,
    ) -> Self {
        // let summing_domain_twiddles = fft::get_twiddles(summing_domain.len());
        // let summing_domain_inv_twiddles = fft::get_inv_twiddles(summing_domain.len());
        RationalSumcheckProver {
            numerator_coeffs,
            denominator_coeffs,
            sigma,
            // summing_domain,
            eta,
            // summing_domain_twiddles,
            // summing_domain_inv_twiddles,
            g_degree,
            e_degree,
            // fractal_options,
            _h: PhantomData,
            current_layer: 0,
        }
    }

    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    pub fn sumcheck_layer_one(
        &self,
        //query: E,
        accumulator: &mut Accumulator<B, E, H>,
        summing_domain: &Vec<B>,
        summing_domain_inv_twiddles: &Vec<B>,
        //channel: &mut DefaultFractalProverChannel<B, E, H>,
        //initial_queries: Vec<usize>,
    ) {
        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        // compute the polynomial e such that e = (Sigma(g, sigma) - summing_poly)/v_H over the summing domain H.
        debug!("Starting a sumcheck proof");

        let _sigma_inv = self.sigma.inv();
        let summing_domain_len = summing_domain.len();
        let summing_domain_e: Vec<E> = summing_domain.iter().map(|x| E::from(*x)).collect();
        //let mut numerator_vals = self.numerator_coeffs.clone();
        //let mut denominator_vals = self.denominator_coeffs.clone();

        // fft::evaluate_poly_with_offset(&mut numerator_vals, &summing_domain_twiddles, self.eta, 1);
        // fft::evaluate_poly_with_offset(&mut denominator_vals, &summing_domain_twiddles, self.eta, 1);
        let numerator_vals = polynom::eval_many(&self.numerator_coeffs, &summing_domain_e);
        let mut denominator_vals = polynom::eval_many(&self.denominator_coeffs, &summing_domain_e);
        denominator_vals = batch_inversion(&denominator_vals);
        let f_hat_evals: Vec<E> = (0..summing_domain_len)
            .into_iter()
            .map(|i| numerator_vals[i] * denominator_vals[i])
            .collect();

        let mut f_hat_coeffs = f_hat_evals;
        pad_with_zeroes(&mut f_hat_coeffs, summing_domain.len());
        fft::interpolate_poly_with_offset(
            &mut f_hat_coeffs,
            &summing_domain_inv_twiddles,
            self.eta,
        );

        let x_coeffs = vec![E::ZERO, E::ONE];
        let sub_factor = self.sigma / E::from(summing_domain.len() as u64);
        let f_hat_minus_sub_factor = polynom::sub(&f_hat_coeffs, &vec![E::from(sub_factor)]);
        assert_eq!(f_hat_minus_sub_factor[0], E::ZERO);
        let g_hat_coeffs = polynom::div(&f_hat_minus_sub_factor, &x_coeffs);

        let e_hat_coeffs = self.compute_e_poly(
            &g_hat_coeffs,
            &self.numerator_coeffs,
            &self.denominator_coeffs,
            self.eta,
            summing_domain_len,
        );
        accumulator.add_polynomial_e(g_hat_coeffs, self.g_degree);
        accumulator.add_polynomial_e(e_hat_coeffs, self.e_degree);
    }

    // SIGMA(g, sigma)(x) = f(x) = p(x)/q(x)
    // SIGMA(g, sigma) = x*g(x) + sigma*|summing_domain|^-1
    // g(x) = x^-1*(f(x) - sigma*|summing_domain|^-1)
    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    pub fn compute_g_poly_on_val(&self, x_val: E, f_x_val: E, summing_domain_len: usize) -> E {
        let dividing_factor_for_sigma: u64 = summing_domain_len.try_into().unwrap();
        let subtracting_factor = self.sigma * E::from(dividing_factor_for_sigma).inv();
        let dividing_factor = x_val.inv();
        dividing_factor * (f_x_val - E::from(subtracting_factor))
    }

    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    pub fn compute_sigma_function_on_val(
        &self,
        x_val: E,
        g_val: E,
        summing_domain_len: usize,
    ) -> E {
        let dividing_factor: u64 = summing_domain_len.try_into().unwrap();
        (x_val * g_val) + (E::from(self.sigma) * E::from(dividing_factor).inv())
    }

    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    pub fn compute_e_poly_on_val(
        &self,
        x_val: E,
        g_val: E,
        summing_poly_numerator_val: E,
        summing_poly_denominator_val: E,
        eta: B,
        summing_domain_len: usize,
    ) -> E {
        let sigma_function = self.compute_sigma_function_on_val(x_val, g_val, summing_domain_len);
        let sigma_minus_f =
            sigma_function * summing_poly_denominator_val - summing_poly_numerator_val;
        let vanishing_on_x = compute_vanishing_poly(x_val, E::from(eta), summing_domain_len);
        sigma_minus_f * vanishing_on_x.inv()
    }

    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    pub fn compute_e_poly(
        &self,
        g_hat_coeffs: &Vec<E>,
        summing_poly_numerator: &Vec<E>,
        summing_poly_denominator: &Vec<E>,
        eta: B,
        summing_domain_len: usize,
    ) -> Vec<E> {
        let dividing_factor: u64 = summing_domain_len.try_into().unwrap();
        let x_func = [E::ZERO, E::ONE];
        let mut sigma_function = polynom::mul(&x_func, &g_hat_coeffs);
        sigma_function[0] =
            sigma_function[0] + (E::from(self.sigma) * E::from(dividing_factor).inv());
        let sigma_minus_f = polynom::sub(
            &polynom::mul(&sigma_function, &summing_poly_denominator),
            &summing_poly_numerator,
        );
        let vanishing_on_x = get_vanishing_poly(E::from(eta), summing_domain_len);
        polynom::div(&sigma_minus_f, &vanishing_on_x)
    }

    fn get_current_layer(&self) -> usize {
        self.current_layer
    }
    fn get_num_layers(&self) -> usize {
        1
    }

    pub fn run_next_layer(
        &mut self,
        _query: E,
        accumulator: &mut Accumulator<B, E, H>,
        summing_domain: &Vec<B>,
        summing_domain_inv_twiddles: &Vec<B>,
    ) -> Result<(), ProverError> {
        if self.get_current_layer() == 0 {
            self.sumcheck_layer_one(accumulator, summing_domain, summing_domain_inv_twiddles);
            self.current_layer += 1;
        }
        Ok(())
    }
}

// impl<
//         B: StarkField,
//         E: FieldElement<BaseField = B>,
//         H: ElementHasher + ElementHasher<BaseField = B>,
//     > LayeredSubProver<B, E, H> for RationalSumcheckProver<B, E, H>
// {
//     fn run_next_layer(
//         &mut self,
//         _query: E,
//         accumulator: &mut Accumulator<B, E, H>,
//         options: &FractalProverOptions<B>,
//     ) -> Result<(), ProverError> {
//         if self.get_current_layer() == 0 {
//             self.sumcheck_layer_one(accumulator);
//             self.current_layer += 1;
//         }
//         Ok(())
//     }
//     fn get_current_layer(&self) -> usize {
//         self.current_layer
//     }
//     fn get_num_layers(&self) -> usize {
//         1
//     }
//     // fn get_fractal_options(&self) -> FractalProverOptions<B> {
//     //     self.fractal_options.clone()
//     // }
// }
