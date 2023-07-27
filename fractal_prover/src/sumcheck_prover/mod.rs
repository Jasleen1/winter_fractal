use std::cmp::max;
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
use winter_math::{fft, log2, FieldElement, StarkField};

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
    eta: B,
    #[allow(dead_code)]
    g_degree: usize,
    e_degree: usize,
    _h: PhantomData<H>,
    current_layer: usize,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    RationalSumcheckProver<B, E, H>
{
    /// Generate a new rational sumcheck prover for fractal
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

    /// This function computes the first layer of the fractal sumcheck
    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    pub fn sumcheck_layer_one(
        &mut self,
        accumulator: &mut Accumulator<B, E, H>,
        domain: &Vec<B>,
        options: &FractalProverOptions<B>,
    ) {
        // compute the polynomial g such that Sigma(g, sigma) = summing_poly
        // compute the polynomial e such that e = (Sigma(g, sigma) - summing_poly)/v_H over the summing domain H.
        debug!("Starting a sumcheck proof");

        let _sigma_inv = self.sigma.inv();

        //todo: don't need to recompute these here. You could try searching options for something the right size?
        let inv_twiddles = fft::get_inv_twiddles(domain.len());

        // the following fft code could be replaced with:
        // let domain_e: Vec<E> = domain.iter().map(|x| E::from(*x)).collect();
        // numerator_vals = polynom::eval_many(&self.numerator_coeffs, &domain_e);
        // denominator_vals = polynom::eval_many(&self.denominator_coeffs, &domain_e);
        // ffts are used for efficiency, even though more evaluations are calculated than necessary sometimes.
        let numerator_vals: Vec<E>;
        let mut denominator_vals: Vec<E>;

        let num_factor = max(
            1,
            self.numerator_coeffs.len().next_power_of_two() / domain.len(),
        );
        // println!("Num factor = {:?}", num_factor);
        // println!("Original = {:?}", self.numerator_coeffs.len());
        pad_with_zeroes(&mut self.numerator_coeffs, num_factor * domain.len());
        let num_twiddles = fft::get_twiddles(num_factor * domain.len());
        let numerator_more_vals =
            fft::evaluate_poly_with_offset(&self.numerator_coeffs, &num_twiddles, self.eta, 1);
        numerator_vals = (0..domain.len())
            .into_iter()
            .map(|i| numerator_more_vals[num_factor * i])
            .collect();

        // save some work if the denominator happens to be a constant
        if self.denominator_coeffs.len() == 1 {
            denominator_vals = vec![self.denominator_coeffs[0]; domain.len()];
        } else {
            let denom_factor = max(
                1,
                self.denominator_coeffs.len().next_power_of_two() / domain.len(),
            );
            pad_with_zeroes(&mut self.denominator_coeffs, denom_factor * domain.len());
            let denom_twiddles = fft::get_twiddles(denom_factor * domain.len());
            let denominator_more_vals = fft::evaluate_poly_with_offset(
                &self.denominator_coeffs,
                &denom_twiddles,
                self.eta,
                1,
            );
            denominator_vals = (0..domain.len())
                .into_iter()
                .map(|i| denominator_more_vals[denom_factor * i])
                .collect();
        }

        // invert all denominator values at once for much cheaper
        let inv_denominator_vals = batch_inversion(&denominator_vals);
        let f_hat_evals: Vec<E> = (0..domain.len())
            .into_iter()
            .map(|i| numerator_vals[i] * inv_denominator_vals[i])
            .collect();

        let mut sum_val = E::ZERO;
        for term in f_hat_evals.clone() {
            sum_val = sum_val + term;
        }

        // println!("sum_val = {:?}", sum_val);
        // println!("sigma = {:?}", self.sigma);

        let mut f_hat_coeffs = f_hat_evals;
        pad_with_zeroes(&mut f_hat_coeffs, domain.len());
        fft::interpolate_poly_with_offset(&mut f_hat_coeffs, &inv_twiddles, self.eta);
        // println!("f_hat degree = {:?}", polynom::degree_of(&f_hat_coeffs));

        let x_coeffs = vec![E::ZERO, E::ONE];
        let sub_factor = self.sigma / E::from(domain.len() as u64);
        let f_hat_minus_sub_factor = polynom::sub(&f_hat_coeffs, &vec![E::from(sub_factor)]);
        assert_eq!(f_hat_minus_sub_factor[0], E::ZERO);
        let g_hat_coeffs = polynom::div(&f_hat_minus_sub_factor, &x_coeffs);

        // let e_hat_coeffs = self.compute_e_poly(
        //     &g_hat_coeffs,
        //     &self.numerator_coeffs,
        //     &self.denominator_coeffs,
        //     self.eta,
        //     domain.len(),
        // );

        // println!("Domain size = {}", domain.len());
        let mut numerator = numerator_vals.clone();
        fft::interpolate_poly_with_offset(&mut numerator, &inv_twiddles, self.eta);

        let mut denominator = denominator_vals.clone();
        fft::interpolate_poly_with_offset(&mut denominator, &inv_twiddles, self.eta);

        let e_hat_coeffs = self.compute_e_poly(
            &g_hat_coeffs,
            &self.numerator_coeffs,
            &self.denominator_coeffs,
            self.eta,
            domain.len(),
        );

        // println!("e actual degree = {:?}", polynom::degree_of(&e_hat_coeffs));
        // println!("e expected degree = {:?}", self.e_degree);
        // println!("g actual degree = {:?}", polynom::degree_of(&g_hat_coeffs));
        // println!("g expected degree = {:?}", self.g_degree);

        accumulator.add_polynomial_e(g_hat_coeffs, self.g_degree);
        accumulator.add_polynomial_e(e_hat_coeffs, self.e_degree);
    }

    // SIGMA(g, sigma)(x) = f(x) = p(x)/q(x)
    // SIGMA(g, sigma) = x*g(x) + sigma*|summing_domain|^-1
    // g(x) = x^-1*(f(x) - sigma*|summing_domain|^-1)
    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    fn compute_g_poly_on_val(&self, x_val: E, f_x_val: E, summing_domain_len: usize) -> E {
        let dividing_factor_for_sigma: u64 = summing_domain_len.try_into().unwrap();
        let subtracting_factor = self.sigma * E::from(dividing_factor_for_sigma).inv();
        let dividing_factor = x_val.inv();
        dividing_factor * (f_x_val - E::from(subtracting_factor))
    }

    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    fn compute_sigma_function_on_val(&self, x_val: E, g_val: E, summing_domain_len: usize) -> E {
        let dividing_factor: u64 = summing_domain_len.try_into().unwrap();
        (x_val * g_val) + (E::from(self.sigma) * E::from(dividing_factor).inv())
    }

    #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    fn compute_e_poly_on_val(
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
    fn compute_e_poly(
        &self,
        g_hat_coeffs: &Vec<E>,
        summing_poly_numerator: &Vec<E>,
        summing_poly_denominator: &Vec<E>,
        eta: B,
        summing_domain_len: usize,
    ) -> Vec<E> {
        let dividing_factor: u64 = summing_domain_len.try_into().unwrap();
        let x_func = [E::ZERO, E::ONE];
        flame::start("mul");
        let mut sigma_function = polynom::mul(&x_func, &g_hat_coeffs);
        flame::end("mul");
        flame::start("inv");
        sigma_function[0] += E::from(self.sigma) * E::from(dividing_factor).inv();
        flame::end("inv");
        flame::start("submul");
        let mut sigma_minus_f = polynom::sub(
            &fft_mul(&sigma_function, &summing_poly_denominator),
            &summing_poly_numerator,
        );
        flame::end("submul");
        divide_by_vanishing_in_place(&mut sigma_minus_f, E::from(eta), summing_domain_len);
        // Now we want to lower the degree to sigma_minus_f to where we want it
        // let summing_domain_base = B::get_root_of_unity(log2(summing_domain_len));
        // let summing_domain = winter_math::get_power_series_with_offset(summing_domain_base, eta, summing_domain_len);
        // let summing_domain_e = summing_domain.iter().map(|x| E::from(*x)).collect::<Vec<E>>();

        // // let mut sigma_minus_f_evals = polynom::eval_many(&sigma_minus_f, &summing_domain_e);
        // // let inv_twiddles = fft::get_inv_twiddles(summing_domain_len);
        // // fft::interpolate_poly_with_offset(&mut sigma_minus_f_evals, &inv_twiddles, eta);

        sigma_minus_f
        // let vanishing_on_x = get_vanishing_poly(E::from(eta), summing_domain_len);
        // polynom::div(&sigma_minus_f, &vanishing_on_x)
    }

    // #[cfg_attr(feature = "flame_it", flame("sumcheck_prover"))]
    // fn compute_e_poly_with_evals(
    //     &self,
    //     g_hat_coeffs: &Vec<E>,
    //     summing_poly_numerator: &Vec<E>,
    //     summing_poly_denominator: &Vec<E>,
    //     eta: B,
    //     summing_domain_len: usize,
    // ) -> Vec<E> {
    //     let dividing_factor: u64 = summing_domain_len.try_into().unwrap();
    //     let x_func = [E::ZERO, E::ONE];
    //     flame::start("mul");
    //     let mut sigma_function = polynom::mul(&x_func, &g_hat_coeffs);
    //     flame::end("mul");
    //     flame::start("inv");
    //     sigma_function[0] += E::from(self.sigma) * E::from(dividing_factor).inv();
    //     flame::end("inv");
    //     // flame::start("submul");
    //     let mut sigma_minus_f = polynom::sub(
    //         &fft_mul(&sigma_function, &summing_poly_denominator),
    //         &summing_poly_numerator,
    //     );
    //     // println!("denom times summing = {:?}", fft_mul(&sigma_function, &summing_poly_denominator));
    //     // println!("neg num =  {:?}", polynom::mul_by_scalar(&summing_poly_numerator, E::ONE.neg()));
    //     println!("Sigma minus f = {:?}", polynom::add(&fft_mul(&sigma_function, &summing_poly_denominator), &polynom::mul_by_scalar(&summing_poly_numerator, E::ONE.neg())));
    //     // let twiddles = fft::get_twiddles(summing_domain_len);
    //     let inv_twiddles = fft::get_inv_twiddles(summing_domain_len);
    //     let summing_domain_base = B::get_root_of_unity(log2(summing_domain_len));
    //     let summing_domain = winter_math::get_power_series_with_offset(summing_domain_base, eta, summing_domain_len);

    //     // fractal_utils::polynomial_utils::pad_with_zeroes(&mut sigma_minus_f, summing_domain_len);
    //     // fft::evaluate_poly(&mut sigma_minus_f, &twiddles);
    //     let summing_domain_e = summing_domain.iter().map(|x| E::from(*x)).collect::<Vec<E>>();

    //     // let sigma_minus_f_evals = polynom::eval_many(&sigma_minus_f, &summing_domain_e);

    //     divide_by_vanishing_in_place(&mut sigma_minus_f, E::from(eta), summing_domain_len);

    //     let mut sigma_minus_f_evals = polynom::eval_many(&sigma_minus_f, &summing_domain_e);

    //     // let div_evals = (0..summing_domain_len)
    //     // .into_iter()
    //     // .map(|i| sigma_minus_f_evals[i] / E::from(vanishing_poly_evals[i]))
    //     // .collect::<Vec::<E>>();

    //     println!("summing domain len = {:?}", summing_domain_e.len());
    //     // println!("div_evals = {:?}", div_evals);
    //     fft::interpolate_poly_with_offset(&mut sigma_minus_f_evals, &inv_twiddles, eta);
    //     // let div_coeffs = polynom::interpolate(&summing_domain_e, &sigma_minus_f_evals, true);

    //     // div_coeffs.to_vec()
    //     sigma_minus_f_evals

    //     // flame::end("submul");
    //     // divide_by_vanishing_in_place(&mut sigma_minus_f, E::from(eta), summing_domain_len);
    //     // sigma_minus_f
    //     //let vanishing_on_x = get_vanishing_poly(E::from(eta), summing_domain_len);
    //     //polynom::div(&sigma_minus_f, &vanishing_on_x)
    // }

    fn get_current_layer(&self) -> usize {
        self.current_layer
    }
    fn get_num_layers(&self) -> usize {
        1
    }

    /// Run the sumcheck next layer
    pub fn run_next_layer(
        &mut self,
        _query: E,
        accumulator: &mut Accumulator<B, E, H>,
        domain: &Vec<B>,
        options: &FractalProverOptions<B>,
    ) -> Result<(), ProverError> {
        if self.get_current_layer() == 0 {
            self.sumcheck_layer_one(accumulator, domain, options);
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
