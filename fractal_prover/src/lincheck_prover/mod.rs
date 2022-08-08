use std::{marker::PhantomData, usize};

use fractal_indexer::{hash_values, snark_keys::*};
use fractal_utils::polynomial_utils::*;

use fractal_sumcheck::sumcheck_prover::*;

use fractal_proofs::{fft, polynom, LincheckProof, OracleQueries, TryInto};

use winter_crypto::{ElementHasher, MerkleTree};
use winter_fri::ProverChannel;
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use crate::{errors::LincheckError, FractalOptions};

const n: usize = 1;
// TODO: Will need to ask Irakliy whether a channel should be passed in here
pub struct LincheckProver<
    'a,
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    alpha: B,
    prover_matrix_index: &'a ProverMatrixIndex<H, B>,
    f_1_poly_coeffs: Vec<B>,
    f_2_poly_coeffs: Vec<B>,
    options: &'a FractalOptions<B>,
    _h: PhantomData<H>,
    _e: PhantomData<E>,
}

impl<
        'a,
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LincheckProver<'a, B, E, H>
{
    pub fn new(
        alpha: B,
        prover_matrix_index: &'a ProverMatrixIndex<H, B>,
        f_1_poly_coeffs: Vec<B>,
        f_2_poly_coeffs: Vec<B>,
        options: &'a FractalOptions<B>,
    ) -> Self {
        LincheckProver {
            alpha,
            prover_matrix_index,
            f_1_poly_coeffs,
            f_2_poly_coeffs,
            options,
            _h: PhantomData,
            _e: PhantomData,
        }
    }

    /// The polynomial t_alpha(X) = u_M(X, alpha). 
    /// We also know that u_M(X, alpha) = M_star(X, alpha).
    /// Further, M_star(X, Y) = 
    /// sum_{k in summing domain} (v_H(X)/ (X - row(k))) * (v_H(Y)/ (Y - col(k))) * val(k).
    /// Fixing Y = alpha, this gives us t_alpha(X) = sum_k (v_H(X)/ (X - row(k))) * (v_H(alpha)/ (alpha - col(k))) * val(k).
    /// = v_H(alpha) * sum_k (v_H(X)/ (X - row(k))) * (val(k)/ (alpha - col(k)))
    pub fn generate_t_alpha_evals(&self) -> Vec<B> {
        // Lets get the coefficients (val(k)/ (alpha - col(k)) 
        // for all values of k, since these don't change with X.
        let mut coefficient_values = Vec::new();
        for id in 0..self.options.summing_domain.len() {
            let summing_elt = self.options.summing_domain[id];
            let denom_term = self.alpha - self.prover_matrix_index.get_col_eval(summing_elt);
            let inv_denom_term = denom_term.inv();
            // This computes the term val(k) / (alpha - col(k))
            let k_term = self.prover_matrix_index.get_val_eval(summing_elt) * inv_denom_term;
            coefficient_values.push(k_term)
        }
        // This is the v_h(alpha) term, which only needs to be computed once.
        let v_h_alpha = vanishing_poly_for_mult_subgroup(self.alpha, self.options.size_subgroup_h);
        // Now we compute the terms sum_k (v_H(X)/ (X - row(k))) * (val(k)/ (alpha - col(k)))
        // over the eval domain.
        let mut t_evals = Vec::new();
        for x_val_id in 0..self.options.evaluation_domain.len() {
            let x_val = self.options.evaluation_domain[x_val_id];
            
            // Getting sum_k (1/ (X - row(k))) * (val(k)/ (alpha - col(k)))
            let mut sum_without_vs = B::ZERO;
            for id in 0..self.options.summing_domain.len() {
                let summing_elt = self.options.summing_domain[id];
                let denom_term = x_val - self.prover_matrix_index.get_row_eval(summing_elt);
                let prod_term = coefficient_values[id] * denom_term.inv();
                sum_without_vs = sum_without_vs + prod_term;
            }
            // This is v_H(X).
            let v_h_x = vanishing_poly_for_mult_subgroup(x_val, self.options.size_subgroup_h);
            // This is finally v_H(X) * v_H(alpha) * sum_K (1/ (X - row(k))) * (val(k)/ (alpha - col(k)))
            let sum_with_vs = (sum_without_vs * v_h_x) * v_h_alpha;
            t_evals.push(sum_with_vs);
        }
        // println!("t_alpha init evals {:?}", t_evals);
        t_evals
    }

    pub fn generate_t_alpha_on_h(&self, t_evals: Vec<B>) -> Vec<B> {
        let mut t_alpha_h_domain_poly: Vec<B> = t_evals.clone();
        // println!("t_alpha on h is {:?}", t_alpha_h_domain_poly);
        // let twiddles_evaluation_domain: Vec<B> =
        //     fft::get_twiddles(self.options.evaluation_domain.len());
        let twiddles_h_domain: Vec<B> =
            fft::get_inv_twiddles(self.options.h_domain.len());
        println!("h twidz len = {:?}", twiddles_h_domain.len());
        println!("t_alpha len = {:?}", t_alpha_h_domain_poly.len());
        fft::interpolate_poly(&mut t_alpha_h_domain_poly, &twiddles_h_domain);
        // println!("t_alpha on h is {:?}", t_alpha_h_domain_poly);
        t_alpha_h_domain_poly
    }

    pub fn generate_t_alpha(&self, t_evals: Vec<B>) -> Vec<B> {
        let mut t_alpha_eval_domain_poly: Vec<B> = t_evals.clone()[0..self.options.h_domain.len()].to_vec();
        let twiddles_evaluation_domain: Vec<B> =
            fft::get_inv_twiddles(self.options.h_domain.len());
        fft::interpolate_poly(&mut t_alpha_eval_domain_poly, &twiddles_evaluation_domain);
        println!("t_alpha = {:?}", t_alpha_eval_domain_poly);
        // fractal_utils::polynomial_utils::get_to_degree_size(&mut t_alpha_eval_domain_poly);
        t_alpha_eval_domain_poly
    }

    pub fn generate_poly_prod(&self, t_alpha_coeffs: &Vec<B>) -> Vec<B> {
        // This function needs to compute the polynomial
        // u_H(X, alpha)*f_1 - t_alpha*f_2
        // here are the steps to this:
        // 1. find out how polynomials are represented and get u_H(X, alpha) = (X^|H| - alpha)/(X - alpha)
        // 2. Polynom includes a mul and a sub function, use these to do the respective ops
        let mut u_numerator = vec![B::ZERO; (self.options.size_subgroup_h).try_into().unwrap()];
        u_numerator[0] = self.alpha.neg();
        u_numerator.push(B::ONE);
        let u_denominator = vec![self.alpha.neg(), B::ONE];
        let mut u_alpha_coeffs = polynom::div(&u_numerator, &u_denominator);
        fractal_utils::polynomial_utils::get_to_degree_size(&mut u_alpha_coeffs);
        println!("u_alpha_len = {}", u_alpha_coeffs.len());
        println!("f_1_len = {}", self.f_1_poly_coeffs.len());
        println!("f_2_len = {}", self.f_2_poly_coeffs.len());
        let mut poly = polynom::sub(
            &polynom::mul(&u_alpha_coeffs, &self.f_1_poly_coeffs),
            &polynom::mul(t_alpha_coeffs, &self.f_2_poly_coeffs),
        );
        fractal_utils::polynomial_utils::get_to_degree_size(&mut poly);
        println!("Poly size = {:?}", poly.len());
        println!("Poly = {:?}", poly);
        poly
    }

    pub fn generate_poly_prod_evals(&self, t_alpha: &Vec<B>) -> Vec<B> {
        // This function needs to compute the polynomial
        // u_H(X, alpha)*f_1 - t_alpha*f_2
        // here are the steps to this:
        // 1. find out how polynomials are represented and get u_H(X, alpha) = (X^|H| - alpha)/(X - alpha)
        // 2. Polynom includes a mul and a sub function, use these to do the respective ops
        let mut u_numerator = vec![B::ZERO; (self.options.size_subgroup_h).try_into().unwrap()];
        u_numerator[0] = self.alpha.neg();
        u_numerator.push(B::ONE);
        let u_denominator = vec![self.alpha.neg(), B::ONE];
        let mut u_alpha = polynom::div(&u_numerator, &u_denominator);
        // fractal_utils::polynomial_utils::get_to_degree_size(&mut u_alpha_coeffs);
        println!("u_alpha_len = {}", u_alpha.len());
        println!("f_1_len = {}", self.f_1_poly_coeffs.len());
        println!("f_2_len = {}", self.f_2_poly_coeffs.len());
        let mut prod = Vec::<B>::new();
        let eval_twiddles = fft::get_twiddles(self.options.evaluation_domain.len());
        let mut f_1_eval = self.f_1_poly_coeffs.clone();
        fractal_utils::polynomial_utils::pad_with_zeroes(&mut f_1_eval, self.options.evaluation_domain.len());
        fft::evaluate_poly(&mut f_1_eval, &mut eval_twiddles.clone());
        let mut f_2_eval = self.f_2_poly_coeffs.clone();
        fractal_utils::polynomial_utils::pad_with_zeroes(&mut f_2_eval, self.options.evaluation_domain.len());
        fft::evaluate_poly(&mut f_2_eval, &mut eval_twiddles.clone());
        fractal_utils::polynomial_utils::pad_with_zeroes(&mut u_alpha, self.options.evaluation_domain.len());
        fft::evaluate_poly(&mut u_alpha, &mut eval_twiddles.clone());
        for pos in 0..self.options.evaluation_domain.len() {
            let next = (u_alpha[pos] * f_1_eval[pos]) - (t_alpha[pos] * f_2_eval[pos]);
            prod.push(next);
        }
        prod
    }

    pub fn generate_lincheck_proof(&self) -> Result<LincheckProof<B, E, H>, LincheckError> {
        let t_alpha_evals = self.generate_t_alpha_evals();
        let mut t_alpha = self.generate_t_alpha(t_alpha_evals.clone());
        // let eval_twiddles = fft::get_twiddles::<B>(self.options.evaluation_domain.len());
        // let mut f_1_eval = self.f_1_poly_coeffs.clone();
        // fractal_utils::polynomial_utils::get_to_degree_size(f_1_eval, self.options.evaluation_domain.len());
        // fft::evaluate_poly(&mut f_1_eval, eval_twiddles);
        // // println!("t_alpha = {:?}", t_alpha);
        // get_to_degree_size(&mut t_alpha);
        // println!("t_alpha_size = {}", t_alpha.len());
        let poly_prod = self.generate_poly_prod_evals(&t_alpha_evals);
        // Next use poly_beta in a sumcheck proof but
        // the sumcheck domain is H, which isn't included here
        // Use that to produce the sumcheck proof.
        println!("Poly prod len = {}", poly_prod.len());
        // println!("Poly prod= {:?}", poly_prod);

        let denom_eval = vec![B::ONE; self.options.evaluation_domain.len()];
        let g_degree = self.options.h_domain.len() - 2;
        let e_degree = self.options.h_domain.len() - 1;
        let g_max_degree = g_degree.next_power_of_two();
        let e_max_degree = e_degree.next_power_of_two();
        let mut product_sumcheck_prover = RationalSumcheckProver::<B, E, H>::new(
            poly_prod,
            denom_eval,
            E::ZERO,
            self.options.h_domain.clone(),
            self.options.evaluation_domain.clone(),
            g_max_degree, 
            e_max_degree,
            self.options.fri_options.clone(),
            self.options.num_queries,
        );
        let products_sumcheck_proof = product_sumcheck_prover.generate_proof();
        let beta =
            FieldElement::as_base_elements(&[product_sumcheck_prover.channel.draw_fri_alpha()])[0];
        let gamma = polynom::eval(&t_alpha, beta);
        let matrix_proof_numerator = polynom::mul_by_scalar(
            &self.prover_matrix_index.val_poly.polynomial,
            compute_vanishing_poly(self.alpha, B::ONE, self.options.size_subgroup_h)
                * compute_vanishing_poly(beta, B::ONE, self.options.size_subgroup_h),
        );
        let mut alpha_minus_row =
            polynom::mul_by_scalar(&self.prover_matrix_index.row_poly.polynomial, -B::ONE);
        alpha_minus_row[0] = alpha_minus_row[0] + self.alpha;
        let mut beta_minus_col =
            polynom::mul_by_scalar(&self.prover_matrix_index.col_poly.polynomial, -B::ONE);
        beta_minus_col[0] = beta_minus_col[0] + beta;
        let matrix_proof_denominator = polynom::mul(&alpha_minus_row, &beta_minus_col);
        // NEXT TODO: Make sure to write the evals for the two polynomials for the next sumcheck
        let matrix_proof_numerator_evals = polynom::eval_many(&matrix_proof_numerator, &self.options.evaluation_domain);
        let matrix_proof_denominator_evals = polynom::eval_many(&matrix_proof_denominator, &self.options.evaluation_domain);
        // println!("Denom deg = {}, num degree = {}, eval_len = {}", matrix_proof_denominator.len(), matrix_proof_numerator.len());
        // println!("Divided poly = {:?}", polynom::div(&matrix_proof_numerator, &matrix_proof_denominator));
        let mut matrix_sumcheck_prover = RationalSumcheckProver::<B, E, H>::new(
            matrix_proof_numerator_evals,
            matrix_proof_denominator_evals,
            E::from(gamma),
            self.options.summing_domain.clone(),
            self.options.evaluation_domain.clone(),
            self.options.summing_domain.len() - 2,
            2 * self.options.summing_domain.len() - 3,
            self.options.fri_options.clone(),
            self.options.num_queries,
        );
        let matrix_sumcheck_proof = matrix_sumcheck_prover.generate_proof();

        let queried_positions = matrix_sumcheck_proof.queried_positions.clone();

        let row_queried_evaluations = queried_positions
            .iter()
            .map(|&p| E::from(self.prover_matrix_index.row_poly.evaluations[p]))
            .collect::<Vec<_>>();
        let row_proofs_results = queried_positions
            .iter()
            .map(|&p| self.prover_matrix_index.row_poly.tree.prove(p))
            .collect::<Vec<_>>();
        let mut row_proofs = Vec::new();
        for row_proof in row_proofs_results {
            if !row_proof.is_ok() {
                println!("row problem: {:?}", row_proof);
            }
            row_proofs.push(row_proof?);
        }
        let row_queried = OracleQueries::<B, E, H>::new(row_queried_evaluations, row_proofs);

        let col_queried_evaluations = queried_positions
            .iter()
            .map(|&p| E::from(self.prover_matrix_index.col_poly.evaluations[p]))
            .collect::<Vec<_>>();
        let col_proofs_results = queried_positions
            .iter()
            .map(|&p| self.prover_matrix_index.col_poly.tree.prove(p))
            .collect::<Vec<_>>();
        let mut col_proofs = Vec::new();
        for col_proof in col_proofs_results {
            if !col_proof.is_ok() {
                println!("col problem: {:?}", col_proof);
            }
            col_proofs.push(col_proof?);
        }
        let col_queried = OracleQueries::<B, E, H>::new(col_queried_evaluations, col_proofs);

        let val_queried_evaluations = queried_positions
            .iter()
            .map(|&p| E::from(self.prover_matrix_index.val_poly.evaluations[p]))
            .collect::<Vec<_>>();
        let val_proofs_results = queried_positions
            .iter()
            .map(|&p| self.prover_matrix_index.val_poly.tree.prove(p))
            .collect::<Vec<_>>();
        let mut val_proofs = Vec::new();
        for val_proof in val_proofs_results {
            if !val_proof.is_ok() {
                println!("val problem: {:?}", val_proof);
            }
            val_proofs.push(val_proof?);
        }
        let val_queried = OracleQueries::<B, E, H>::new(val_queried_evaluations, val_proofs);

        let t_alpha_transposed_evaluations = transpose_slice::<_, { n }>(&t_alpha_evals.clone());
        let hashed_evaluations = hash_values::<H, B, { n }>(&t_alpha_transposed_evaluations);
        let t_alpha_tree = MerkleTree::<H>::new(hashed_evaluations)?;
        let t_alpha_commitment = *t_alpha_tree.root();
        let t_alpha_queried_evaluations = queried_positions
            .iter()
            .map(|&p| E::from(t_alpha_evals[p]))
            .collect::<Vec<_>>();
        let t_alpha_proofs_results = queried_positions
            .iter()
            .map(|&p| t_alpha_tree.prove(p))
            .collect::<Vec<_>>();
        let mut t_alpha_proofs = Vec::new();
        for t_alpha_proof in t_alpha_proofs_results {
            if !t_alpha_proof.is_ok() {
                println!("T alpha problem: {:?}", t_alpha_proof);
            }
            t_alpha_proofs.push(t_alpha_proof?);
        }
        let t_alpha_queried =
            OracleQueries::<B, E, H>::new(t_alpha_queried_evaluations, t_alpha_proofs);
        Ok(LincheckProof::<B, E, H> {
            options: self.options.fri_options.clone(),
            num_evaluations: self.options.evaluation_domain.len(),
            alpha: self.alpha,
            beta,
            t_alpha_commitment,
            t_alpha_queried,
            products_sumcheck_proof,
            gamma,
            row_queried,
            col_queried,
            val_queried,
            matrix_sumcheck_proof,
            _e: PhantomData,
        })
        // unimplemented!()
    }
}
