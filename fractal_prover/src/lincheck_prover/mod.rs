use std::{marker::PhantomData, usize, sync::Arc, collections::HashMap, hash::BuildHasherDefault};

use fractal_indexer::{hash_values, index::IndexParams, snark_keys::*};
use fractal_utils::polynomial_utils::*;
use models::r1cs::Matrix;
use nohash_hasher::NoHashHasher;
use rustc_hash::FxHashMap;

use crate::{errors::ProverError, sumcheck_prover::*, LayeredProver, LayeredSubProver};
use fractal_accumulator::accumulator::Accumulator;
use fractal_utils::channel::DefaultFractalProverChannel;

use fractal_proofs::{
    batch_inversion, fft, polynom, LayeredLincheckProof, LincheckProof, OracleQueries,
    TopLevelProof, TryInto,
};

use fractal_utils::FractalProverOptions;
use winter_crypto::{
    BatchMerkleProof, ElementHasher, Hasher, MerkleTree, MerkleTreeError, RandomCoin,
};
use winter_fri::ProverChannel;
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use crate::{errors::LincheckError, log::debug};

const n: usize = 1;
/// This is the modular prover for Fractal's Lincheck.
pub struct LincheckProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    prover_matrix_index: Arc<ProverMatrixIndex<B, E>>,
    f_1_poly_coeffs: Vec<B>,
    f_2_poly_coeffs: Vec<B>,
    _h: PhantomData<H>,
    _e: PhantomData<E>,
    current_layer: usize,
    t_alpha: Option<Vec<E>>,
    alpha: Option<E>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LincheckProver<B, E, H>
{
    /// Create a new fractal lincheck prover
    pub fn new(
        prover_matrix_index: Arc<ProverMatrixIndex<B, E>>,
        f_1_poly_coeffs: Vec<B>,
        f_2_poly_coeffs: Vec<B>,
    ) -> Self {
        LincheckProver {
            prover_matrix_index: prover_matrix_index,
            f_1_poly_coeffs,
            f_2_poly_coeffs,
            _h: PhantomData,
            _e: PhantomData,
            current_layer: 0,
            t_alpha: None,
            alpha: None,
        }
    }

    #[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    fn lincheck_layer_one(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
        options: &FractalProverOptions<B>,
    ) -> Result<(), ProverError> {
        self.alpha = Some(query);
        let t_alpha = self.generate_t_alpha(query, options);
        debug!("t_alpha degree: {}", &t_alpha.len() - 1);
        accumulator.add_polynomial_e(t_alpha.clone(), options.size_subgroup_h - 1);
        self.t_alpha = Some(t_alpha.clone());

        let poly_prod_coeffs = self.generate_poly_prod(query, &t_alpha, options);
        debug!(
            "poly_prod_coeffs degree {}",
            polynom::degree_of(&poly_prod_coeffs)
        );

        let g_degree = options.h_domain.len() - 2;
        let e_degree = options.h_domain.len() - 1;

        let mut product_sumcheck_prover = RationalSumcheckProver::<B, E, H>::new(
            poly_prod_coeffs.clone(),
            vec![E::ONE],
            E::ZERO,
            options.eta,
            g_degree,
            e_degree,
        );

        product_sumcheck_prover.run_next_layer(query, accumulator, &options.h_domain, options)?;
        Ok(())
    }

    #[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    fn lincheck_layer_two(
        &self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
        options: &FractalProverOptions<B>,
    ) {
        let beta = query;
        let alpha = self.alpha.unwrap();
        debug!(
            "Alpha and beta in lincheck prover l2 = {:?}, {:?}",
            alpha, beta
        );
        // t_alpha is the only state we need to retain from layer 1
        // if we wanted to be really fancy, we could extract this from the accumulator...
        let gamma = polynom::eval(&self.t_alpha.as_ref().unwrap(), beta);
        let matrix_proof_numerator = polynom::mul_by_scalar(
            &self
                .prover_matrix_index
                .val_poly
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
            compute_vanishing_poly(alpha, E::from(options.eta), options.size_subgroup_h)
                * compute_vanishing_poly(beta, E::from(options.eta), options.size_subgroup_h),
        );
        let mut alpha_minus_row =
            polynom::mul_by_scalar(&self.prover_matrix_index.row_poly, -B::ONE)
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>();
        alpha_minus_row[0] += alpha;
        let mut beta_minus_col =
            polynom::mul_by_scalar(&self.prover_matrix_index.col_poly, -B::ONE)
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>();
        beta_minus_col[0] += beta;

        let mut alpha_minus_col =
            polynom::mul_by_scalar(&self.prover_matrix_index.col_poly, -B::ONE)
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>();
        alpha_minus_col[0] += alpha;
        let mut beta_minus_row =
            polynom::mul_by_scalar(&self.prover_matrix_index.row_poly, -B::ONE)
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>();
        beta_minus_row[0] += beta;

        //let matrix_proof_denominator = polynom::mul(&alpha_minus_row, &beta_minus_col);
        let matrix_proof_denominator = fft_mul(&alpha_minus_col, &beta_minus_row);

        //matrix_proof_numerator/matrix_proof_denominator should evaluate to gamma when summed over K. Let's double check this
        // let mut mat_sum = E::ZERO;
        // for k in self.options.summing_domain.iter() {
        //     let temp = polynom::eval(&matrix_proof_numerator, E::from(*k))
        //         / polynom::eval(&matrix_proof_denominator, E::from(*k));
        //     mat_sum += temp;
        // }

        let mut matrix_sumcheck_prover = RationalSumcheckProver::<B, E, H>::new(
            matrix_proof_numerator,
            matrix_proof_denominator,
            gamma,
            options.eta_k,
            options.summing_domain.len() - 2,
            2 * options.summing_domain.len() - 3,
        );

        matrix_sumcheck_prover
            .run_next_layer(query, accumulator, &options.summing_domain, options)
            .unwrap();
    }

    pub(crate) fn retrieve_gamma(&self, beta: E) -> Result<E, LincheckError> {
        let t_alpha = self
            .t_alpha
            .clone()
            .ok_or(LincheckError::GammaCompErr("t_alpha not set".to_string()))?;
        Ok(polynom::eval(&t_alpha, beta))
    }

    /*#[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    pub fn generate_t_alpha(&self, t_evals: Vec<E>, options: &FractalProverOptions<B>) -> Vec<E> {
        let mut t_alpha_eval_domain_poly: Vec<E> = t_evals.clone().to_vec();
        // let twiddles_evaluation_domain: Vec<B> =
        //     fft::get_inv_twiddles(options.evaluation_domain.len());
        fft::interpolate_poly(
            &mut t_alpha_eval_domain_poly,
            &options.l_domain_inv_twiddles,
        );
        t_alpha_eval_domain_poly
    }*/

    /*#[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    pub fn generate_t_alpha_evals(&self, alpha: E, options: &FractalProverOptions<B>) -> Vec<E> {
        let poly = self.generate_t_alpha2(alpha, options);
        let evaluation_domain_e: Vec<E> = options
             .evaluation_domain
             .iter()
             .map(|i| E::from(*i))
             .collect();
        polynom::eval_many(&poly, &evaluation_domain_e)
    }*/

    /// The polynomial t_alpha(X) = u_M(X, alpha).
    /// We also know that u_M(X, alpha) = M_star(X, alpha).
    /// Further, M_star(X, Y) =
    /// sum_{k in summing domain} (v_H(X)/ (X - row(k))) * (v_H(Y)/ (Y - col(k))) * val(k).
    /// Fixing Y = alpha, this gives us t_alpha(X) = sum_k (v_H(X)/ (X - row(k))) * (v_H(alpha)/ (alpha - col(k))) * val(k).
    /// = v_H(alpha) * sum_k (v_H(X)/ (X - row(k))) * (val(k)/ (alpha - col(k)))
    #[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    fn generate_t_alpha(&self, alpha: E, options: &FractalProverOptions<B>) -> Vec<E> {
        let v_h_alpha =
            compute_vanishing_poly(alpha.clone(), E::from(options.eta), options.size_subgroup_h);
        let v_h_x = get_vanishing_poly(options.eta, options.size_subgroup_h);

        let summing_twiddles = fft::get_twiddles(options.summing_domain.len());

        let col_evals = fft::evaluate_poly_with_offset(
            &self.prover_matrix_index.col_poly,
            &summing_twiddles,
            options.eta_k,
            1,
        );
        let val_evals = fft::evaluate_poly_with_offset(
            &self.prover_matrix_index.val_poly,
            &summing_twiddles,
            options.eta_k,
            1,
        );
        let row_evals = fft::evaluate_poly_with_offset(
            &self.prover_matrix_index.row_poly,
            &summing_twiddles,
            options.eta_k,
            1,
        );

        let mut denom_terms: Vec<E> = col_evals
            .iter()
            .map(|col_eval| alpha - E::from(*col_eval))
            .collect();
        denom_terms = batch_inversion(&denom_terms);
        // This computes the term val(k) / (alpha - col(k))
        let coefficient_values: Vec<E> = (0..options.summing_domain.len())
            .into_iter()
            .map(|id| E::from(val_evals[id]) * denom_terms[id])
            .collect();


        // For efficiency, we compute t_alpha as a evaluations over the H domain, as this allows us to skip most of the computation
        // Below is a commented, slow version of what we're trying to do (but which is maybe easier to understand)
        /*let mut evals_h = vec![];
        for h_elt in options.h_domain.iter(){
            let mut val = E::ZERO;
            for k_idx in 0..options.summing_domain.len(){
                if h_elt != &row_evals[k_idx]{
                    continue;
                }
                val += E::from(compute_derivative_on_single_val(h_elt.clone(), options.h_domain.len() as u128)) * coefficient_values[k_idx];
            }
            evals_h.push(val);
        }*/

        // Instead of a double loop, use a hashmap to be able to look up which h_domain element a given row_poly evaluation is equal to
        // As E doesn't implement Hash, we need to hash its bytes representation instead
        let mut locations = FxHashMap::<&[u8], usize>::default();
        let _: Vec<_> = options.h_domain.iter().enumerate().map(|(i, h)| locations.insert(h.as_bytes(), i)).collect();

        let mut evals_h = vec![E::ZERO; options.h_domain.len()];
        for k_idx in 0..options.summing_domain.len(){
            let h_idx = *locations.get(row_evals[k_idx].as_bytes()).unwrap();
            evals_h[h_idx] += E::from(compute_derivative_on_single_val(row_evals[k_idx], options.h_domain.len() as u128)) * coefficient_values[k_idx];
        }

        fft::interpolate_poly_with_offset(&mut evals_h, &options.h_domain_inv_twiddles, options.eta);
        polynom::mul_by_scalar(&evals_h, v_h_alpha)
    }

    /*#[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    pub fn generate_t_alpha_evals(&self, alpha: E, options: &FractalProverOptions<B>) -> Vec<E> {
        // Lets get the coefficients (val(k)/ (alpha - col(k))
        // for all values of k, since these don't change with X.
        let col_evals =
            polynom::eval_many(&self.prover_matrix_index.col_poly, &options.summing_domain);
        let val_evals =
            polynom::eval_many(&self.prover_matrix_index.val_poly, &options.summing_domain);
        let mut denom_terms: Vec<E> = col_evals
            .iter()
            .map(|col_eval| alpha - E::from(*col_eval))
            .collect();
        denom_terms = batch_inversion(&denom_terms);
        // This computes the term val(k) / (alpha - col(k))
        let coefficient_values: Vec<E> = (0..options.summing_domain.len())
            .into_iter()
            .map(|id| E::from(val_evals[id]) * denom_terms[id])
            .collect();

        // This is the v_h(alpha) term, which only needs to be computed once.
        let v_h_alpha =
            compute_vanishing_poly(alpha.clone(), E::from(options.eta), options.size_subgroup_h);
        //let v_h_alpha = vanishing_poly_for_mult_subgroup(self.alpha, self.options.size_subgroup_h);
        // Now we compute the terms sum_k (v_H(X)/ (X - row(k))) * (val(k)/ (alpha - col(k)))
        // over the eval domain.
        let mut t_evals = Vec::new();
        let row_evals =
            polynom::eval_many(&self.prover_matrix_index.row_poly, &options.summing_domain);

        flame::start("double loop");
        for x_val_id in 0..options.evaluation_domain.len() {
            let x_val = options.evaluation_domain[x_val_id];

            // Getting sum_k (1/ (X - row(k))) * (val(k)/ (alpha - col(k)))
            let mut sum_without_vs = E::ZERO;
            let mut denom_terms: Vec<B> =
                row_evals.iter().map(|row_eval| x_val - *row_eval).collect();
            denom_terms = batch_inversion(&denom_terms);
            for id in 0..options.summing_domain.len() {
                let prod_term = coefficient_values[id] * E::from(denom_terms[id]);
                sum_without_vs += prod_term;
            }
            // This is v_H(X).
            let v_h_x = compute_vanishing_poly(x_val, options.eta, options.size_subgroup_h);
            //let v_h_x = vanishing_poly_for_mult_subgroup(x_val, self.options.size_subgroup_h);
            // This is finally v_H(X) * v_H(alpha) * sum_K (1/ (X - row(k))) * (val(k)/ (alpha - col(k)))
            let sum_with_vs = (sum_without_vs * E::from(v_h_x)) * v_h_alpha;
            t_evals.push(sum_with_vs);
        }
        flame::end("double loop");
        t_evals
    }*/

    /// This function needs to compute the polynomial
    /// u_H(X, alpha)*f_1 - t_alpha*f_2
    #[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    fn generate_poly_prod(
        &self,
        alpha: E,
        t_alpha_coeffs: &Vec<E>,
        options: &FractalProverOptions<B>,
    ) -> Vec<E> {
        // here are the steps to this:
        // 1. find out how polynomials are represented and get u_H(X, alpha) = (X^|H| - alpha)/(X - alpha)
        // 2. Polynom includes a mul and a sub function, use these to do the respective ops
        //eta_to_h_size = eta.exp(B::PositiveInteger::from(self.options.size_subgroup_h));
        let alpha_to_h_size = alpha.exp(E::PositiveInteger::from(options.size_subgroup_h as u64));
        debug!("alpha_to_h_size: {}", &alpha_to_h_size);
        let mut u_numerator = vec![E::ZERO; (options.size_subgroup_h).try_into().unwrap()];
        u_numerator[0] = alpha_to_h_size.neg();
        u_numerator.push(E::ONE);
        /*let u_alpha_evals: Vec<E> = self
        .options
        .evaluation_domain
        .iter()
        .map(|e| polynom::eval(&u_numerator, E::from(*e)) / polynom::eval(&u_denominator, E::from(*e)))
        .collect();*/
        //let u_alpha_coeffs2 =
        //polynom::interpolate(&self.options.evaluation_domain, &u_alpha_evals, true);
        //let u_alpha_coeffs = polynom::div(&u_numerator, &u_denominator);

        // equivalent to polynom::div(&u_numerator, &vec![alpha.neg(), E::ONE]), but faster
        let u_alpha_coeffs = polynom::syn_div(&u_numerator, 1, alpha);
        //let reconstituted = polynom::mul(&u_alpha_coeffs, &u_denominator);

        flame::start("submul");
        let mut poly = polynom::sub(
            &fft_mul(
                &u_alpha_coeffs,
                &self
                    .f_1_poly_coeffs
                    .iter()
                    .map(|i| E::from(*i))
                    .collect::<Vec<E>>(),
            ),
            &fft_mul(
                t_alpha_coeffs,
                &self
                    .f_2_poly_coeffs
                    .iter()
                    .map(|i| E::from(*i))
                    .collect::<Vec<E>>(),
            ),
        );
        flame::end("submul");

        fractal_utils::polynomial_utils::get_to_degree_size(&mut poly);

        poly
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredSubProver<B, E, H> for LincheckProver<B, E, H>
{
    fn run_next_layer(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
        options: &FractalProverOptions<B>,
    ) -> Result<(), ProverError> {
        match self.get_current_layer() {
            0 => {
                self.lincheck_layer_one(query, accumulator, options)?;
            }
            1 => {
                self.lincheck_layer_two(query, accumulator, options);
            }
            _ => (),
        };
        self.current_layer += 1;
        Ok(())
    }
    fn get_current_layer(&self) -> usize {
        self.current_layer.clone()
    }

    fn get_num_layers(&self) -> usize {
        2
    }

    fn get_max_degree_constraint(
        num_input_variables: usize,
        num_non_zero: usize,
        num_constraints: usize,
    ) -> usize {
        let summing_domain_len = num_non_zero;
        let h_domain_len = std::cmp::max(num_input_variables, num_constraints);
        let v = vec![
            h_domain_len - 2,           //product sumcheck g_degree
            summing_domain_len - 2,     //matrix sumcheck g_degree
            2 * summing_domain_len - 3, //matrix sumcheck e_degree
        ];
        v.iter().max().unwrap().next_power_of_two()
    }

    // fn get_fractal_options(&self) -> FractalProverOptions<B> {
    //     self.options.clone()
    // }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredProver<B, E, H, LayeredLincheckProof<B, E>> for LincheckProver<B, E, H>
{
    #[cfg_attr(feature = "flame_it", flame("lincheck_prover"))]
    fn generate_proof(
        &mut self,
        prover_key: &Option<ProverKey<B, E, H>>,
        public_inputs_bytes: Vec<u8>,
        options: &FractalProverOptions<B>,
    ) -> Result<TopLevelProof<B, E, H>, ProverError> {
        // let options = self.get_fractal_options();
        let mut coin = RandomCoin::<B, H>::new(&public_inputs_bytes);

        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            options.evaluation_domain.len(),
            options.num_queries,
            public_inputs_bytes.clone(),
        );
        let mut acc = Accumulator::<B, E, H>::new(
            options.evaluation_domain.len(),
            B::ONE,
            options.evaluation_domain.clone(),
            options.num_queries,
            options.fri_options.clone(),
            public_inputs_bytes,
            prover_key.as_ref().unwrap().params.max_degree,
        );

        acc.add_unchecked_polynomial(self.f_2_poly_coeffs.clone());
        acc.add_unchecked_polynomial(self.f_1_poly_coeffs.clone());
        let initial_commitment = acc.commit_layer()?;

        let mut layer_commitments = vec![];
        let mut local_queries = Vec::<E>::new();

        for i in 0..self.get_num_layers() {
            let previous_commit = acc.get_layer_commitment(i + 1)?;
            channel.commit_fractal_iop_layer(previous_commit);
            coin.reseed(previous_commit);

            let query = coin.draw().expect("failed to draw FRI alpha"); //channel.draw_fri_alpha();
            local_queries.push(query);
            self.run_next_layer(query, &mut acc, options)?;
            layer_commitments.push(acc.commit_layer()?); //todo: do something with this
        }

        let queries = acc.draw_query_positions()?;

        let initial_decommitment = acc.decommit_layer_with_queries(1, &queries)?;
        let layer_decommits = vec![
            acc.decommit_layer_with_queries(2, &queries)?,
            acc.decommit_layer_with_queries(3, &queries)?,
        ];
        let preprocessing_decommitment = prover_key
            .as_ref()
            .unwrap()
            .accumulator
            .decommit_layer_with_queries(1, &queries)?;

        let beta = local_queries[1];

        println!("Prover alpha?, beta: {}, {}", &local_queries[0], &beta);
        let gammas = vec![self.retrieve_gamma(beta)?];

        let low_degree_proof = acc.create_fri_proof()?;

        let proof = TopLevelProof {
            preprocessing_decommitment,
            layer_commitments: layer_commitments.to_vec(),
            layer_decommitments: layer_decommits,
            initial_commitment,
            initial_decommitment,
            unverified_misc: gammas,
            low_degree_proof,
        };
        Ok(proof)
    }
}
