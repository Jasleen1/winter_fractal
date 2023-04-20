use std::{marker::PhantomData, usize};

use fractal_indexer::{hash_values, snark_keys::*};
use fractal_utils::polynomial_utils::*;
use models::r1cs::Matrix;

use crate::{
    errors::ProverError,
    sumcheck_prover::*, LayeredSubProver, LayeredProver,
};
use fractal_accumulator::accumulator::Accumulator;
use fractal_utils::channel::DefaultFractalProverChannel;

use fractal_proofs::{fft, polynom, LincheckProof, OracleQueries, TryInto, LayeredLincheckProof, TopLevelProof};

use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree, MerkleTreeError, RandomCoin, Hasher};
use winter_fri::ProverChannel;
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;
use fractal_utils::FractalOptions;

use crate::{errors::LincheckError, log::debug};

const n: usize = 1;
// TODO: Will need to ask Irakliy whether a channel should be passed in here
pub struct LincheckProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    prover_matrix_index: ProverMatrixIndex<B,E>,
    f_1_poly_coeffs: Vec<B>,
    f_2_poly_coeffs: Vec<B>,
    options: FractalOptions<B>,
    evaluation_domain_e: Vec<E>,
    _h: PhantomData<H>,
    _e: PhantomData<E>,
    current_layer: usize,
    // layered prover state
    t_alpha: Option<Vec<E>>,
    alpha: Option<E>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LincheckProver<B, E, H>
{
    pub fn new(
        prover_matrix_index: ProverMatrixIndex<B,E>,
        f_1_poly_coeffs: Vec<B>,
        f_2_poly_coeffs: Vec<B>,
        options: &FractalOptions<B>,
    ) -> Self {
        let evaluation_domain_e = options
            .evaluation_domain
            .iter()
            .map(|i| E::from(*i))
            .collect();
        LincheckProver {
            prover_matrix_index: prover_matrix_index,
            f_1_poly_coeffs,
            f_2_poly_coeffs,
            options: options.clone(),
            evaluation_domain_e,
            _h: PhantomData,
            _e: PhantomData,
            current_layer: 0,
            t_alpha: None,
            alpha: None,
        }
    }

    fn lincheck_layer_one(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError> {
        self.alpha = Some(query);
        let t_alpha_evals = self.generate_t_alpha_evals(query);
        let t_alpha = self.generate_t_alpha(t_alpha_evals.clone());
        debug!("t_alpha degree: {}", &t_alpha.len() - 1);
        accumulator.add_polynomial_e(t_alpha.clone(), self.options.size_subgroup_h - 1);
        self.t_alpha = Some(t_alpha.clone());

        let poly_prod = self.generate_poly_prod_evals(query, &t_alpha_evals);
        let poly_prod_coeffs = self.generate_poly_prod(query, &t_alpha);
        debug!(
            "poly_prod_coeffs degree {}",
            polynom::degree_of(&poly_prod_coeffs)
        );

        //poly_prod_coeffs should evaluate to 0 when summed over H. Let's double check this
        let mut pp_sum = E::ZERO;
        for h in self.options.h_domain.iter() {
            let temp = polynom::eval(&poly_prod_coeffs, E::from(*h));
            pp_sum += temp;
        }
        debug_assert!(
            pp_sum == E::ZERO,
            "Sum of product polynomials over h domain is not 0"
        );

        // Next use poly_beta in a sumcheck proof but
        // the sumcheck domain is H, which isn't included here
        // Use that to produce the sumcheck proof.
        debug!("Poly prod len = {}", poly_prod.len());

        //let denom_eval = vec![B::ONE; self.options.evaluation_domain.len()];
        let denom_eval = vec![B::ONE; self.options.h_domain.len()];

        // use h_domain rather than eval_domain
        // let poly_prod = polynom::eval_many(
        //     &poly_prod_coeffs,
        //     &self
        //         .options
        //         .h_domain
        //         .iter()
        //         .map(|i| E::from(*i))
        //         .collect::<Vec<E>>(),
        // );

        let g_degree = self.options.h_domain.len() - 2;
        let e_degree = self.options.h_domain.len() - 1;

        let mut product_sumcheck_prover = RationalSumcheckProver::<B, E, H>::new(
            poly_prod_coeffs.clone(),
            vec![E::ONE],
            E::ZERO,
            self.options.h_domain.clone(),
            self.options.eta,
            g_degree,
            e_degree,
            self.options.clone(),
        );
        //if this needs a channel... problem
        product_sumcheck_prover.run_next_layer(query, accumulator)?;
        Ok(())
    }

    fn lincheck_layer_two(&self, query: E, accumulator: &mut Accumulator<B, E, H>) {
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
            compute_vanishing_poly(
                alpha,
                E::from(self.options.eta),
                self.options.size_subgroup_h,
            ) * compute_vanishing_poly(
                beta,
                E::from(self.options.eta),
                self.options.size_subgroup_h,
            ),
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
        let matrix_proof_denominator = polynom::mul(&alpha_minus_col, &beta_minus_row);

        //matrix_proof_numerator/matrix_proof_denominator should evaluate to gamma when summed over K. Let's double check this
        let mut mat_sum = E::ZERO;
        for k in self.options.summing_domain.iter() {
            let temp = polynom::eval(&matrix_proof_numerator, E::from(*k))
                / polynom::eval(&matrix_proof_denominator, E::from(*k));
            mat_sum += temp;
        }

        let mut matrix_sumcheck_prover = RationalSumcheckProver::<B, E, H>::new(
            matrix_proof_numerator,
            matrix_proof_denominator,
            gamma,
            self.options.summing_domain.clone(),
            self.options.eta_k,
            self.options.summing_domain.len() - 2,
            2 * self.options.summing_domain.len() - 3,
            self.options.clone(),
        );

        matrix_sumcheck_prover
            .run_next_layer(query, accumulator)
            .unwrap();
    }

    pub fn retrieve_gamma(&self, beta: E) -> Result<E, LincheckError> {
        let t_alpha = self
            .t_alpha
            .clone()
            .ok_or(LincheckError::GammaCompErr("t_alpha not set".to_string()))?;
        Ok(polynom::eval(&t_alpha, beta))
    }

    /// The polynomial t_alpha(X) = u_M(X, alpha).
    /// We also know that u_M(X, alpha) = M_star(X, alpha).
    /// Further, M_star(X, Y) =
    /// sum_{k in summing domain} (v_H(X)/ (X - row(k))) * (v_H(Y)/ (Y - col(k))) * val(k).
    /// Fixing Y = alpha, this gives us t_alpha(X) = sum_k (v_H(X)/ (X - row(k))) * (v_H(alpha)/ (alpha - col(k))) * val(k).
    /// = v_H(alpha) * sum_k (v_H(X)/ (X - row(k))) * (val(k)/ (alpha - col(k)))
    pub fn generate_t_alpha_evals(&self, alpha: E) -> Vec<E> {
        // Lets get the coefficients (val(k)/ (alpha - col(k))
        // for all values of k, since these don't change with X.
        let mut coefficient_values = Vec::new();
        for id in 0..self.options.summing_domain.len() {
            let summing_elt = E::from(self.options.summing_domain[id]);
            let denom_term = alpha - self.prover_matrix_index.get_col_eval(summing_elt);
            let inv_denom_term = denom_term.inv();
            // This computes the term val(k) / (alpha - col(k))
            // Why does this type as B instead of E?
            let k_term =
                E::from(self.prover_matrix_index.get_val_eval(summing_elt)) * inv_denom_term;
            coefficient_values.push(k_term)
        }
        // This is the v_h(alpha) term, which only needs to be computed once.
        let v_h_alpha = compute_vanishing_poly(
            alpha.clone(),
            E::from(self.options.eta),
            self.options.size_subgroup_h,
        );
        //let v_h_alpha = vanishing_poly_for_mult_subgroup(self.alpha, self.options.size_subgroup_h);
        // Now we compute the terms sum_k (v_H(X)/ (X - row(k))) * (val(k)/ (alpha - col(k)))
        // over the eval domain.
        let mut t_evals = Vec::new();
        for x_val_id in 0..self.options.evaluation_domain.len() {
            let x_val = self.options.evaluation_domain[x_val_id];

            // Getting sum_k (1/ (X - row(k))) * (val(k)/ (alpha - col(k)))
            let mut sum_without_vs = E::ZERO;
            for id in 0..self.options.summing_domain.len() {
                //summing \n summing
                let summing_elt = E::from(self.options.summing_domain[id]);
                let denom_term =
                    E::from(x_val) - self.prover_matrix_index.get_row_eval(summing_elt);
                let prod_term = coefficient_values[id] * denom_term.inv();
                sum_without_vs = sum_without_vs + prod_term;
            }
            // This is v_H(X).
            let v_h_x =
                compute_vanishing_poly(x_val, self.options.eta, self.options.size_subgroup_h);
            //let v_h_x = vanishing_poly_for_mult_subgroup(x_val, self.options.size_subgroup_h);
            // This is finally v_H(X) * v_H(alpha) * sum_K (1/ (X - row(k))) * (val(k)/ (alpha - col(k)))
            let sum_with_vs = (sum_without_vs * E::from(v_h_x)) * v_h_alpha;
            t_evals.push(sum_with_vs);
        }
        t_evals
    }

    pub fn generate_t_alpha_on_h(&self, t_evals: Vec<B>) -> Vec<B> {
        let mut t_alpha_h_domain_poly: Vec<B> = t_evals.clone();

        let twiddles_h_domain: Vec<B> = fft::get_inv_twiddles(self.options.h_domain.len());

        fft::interpolate_poly(&mut t_alpha_h_domain_poly, &twiddles_h_domain);

        t_alpha_h_domain_poly
    }

    pub fn generate_t_alpha(&self, t_evals: Vec<E>) -> Vec<E> {
        let mut t_alpha_eval_domain_poly: Vec<E> =
            t_evals.clone()[0..self.options.h_domain.len()].to_vec();
        let twiddles_evaluation_domain: Vec<B> = fft::get_inv_twiddles(self.options.h_domain.len());
        polynom::interpolate(
            &self
                .options
                .evaluation_domain
                .iter()
                .map(|i| E::from(*i))
                .collect::<Vec<E>>(),
            &t_evals.to_vec(),
            true,
        )
    }

    pub fn generate_poly_prod(&self, alpha: E, t_alpha_coeffs: &Vec<E>) -> Vec<E> {
        // This function needs to compute the polynomial
        // u_H(X, alpha)*f_1 - t_alpha*f_2
        // here are the steps to this:
        // 1. find out how polynomials are represented and get u_H(X, alpha) = (X^|H| - alpha)/(X - alpha)
        // 2. Polynom includes a mul and a sub function, use these to do the respective ops
        //eta_to_h_size = eta.exp(B::PositiveInteger::from(self.options.size_subgroup_h));
        let alpha_to_h_size = alpha.exp(E::PositiveInteger::from(
            self.options.size_subgroup_h as u64,
        ));
        debug!("alpha_to_h_size: {}", &alpha_to_h_size);
        let mut u_numerator = vec![E::ZERO; (self.options.size_subgroup_h).try_into().unwrap()];
        u_numerator[0] = alpha_to_h_size.neg();
        u_numerator.push(E::ONE);
        let u_denominator = vec![alpha.neg(), E::ONE];
        /*let u_alpha_evals: Vec<E> = self
        .options
        .evaluation_domain
        .iter()
        .map(|e| polynom::eval(&u_numerator, E::from(*e)) / polynom::eval(&u_denominator, E::from(*e)))
        .collect();*/
        //let u_alpha_coeffs2 =
        //polynom::interpolate(&self.options.evaluation_domain, &u_alpha_evals, true);
        let u_alpha_coeffs = polynom::div(&u_numerator, &u_denominator);
        //let reconstituted = polynom::mul(&u_alpha_coeffs, &u_denominator);

        let mut poly = polynom::sub(
            &polynom::mul(
                &u_alpha_coeffs,
                &self
                    .f_1_poly_coeffs
                    .iter()
                    .map(|i| E::from(*i))
                    .collect::<Vec<E>>(),
            ),
            &polynom::mul(
                t_alpha_coeffs,
                &self
                    .f_2_poly_coeffs
                    .iter()
                    .map(|i| E::from(*i))
                    .collect::<Vec<E>>(),
            ),
        );

        fractal_utils::polynomial_utils::get_to_degree_size(&mut poly);

        poly
    }

    pub fn generate_poly_prod_evals(&self, alpha: E, t_alpha: &Vec<E>) -> Vec<E> {
        // This function needs to compute the polynomial
        // u_H(X, alpha)*f_1 - t_alpha*f_2
        // here are the steps to this:
        // 1. find out how polynomials are represented and get u_H(X, alpha) = (X^|H| - alpha)/(X - alpha)
        // 2. Polynom includes a mul and a sub function, use these to do the respective ops
        // botttom of page 29
        let alpha_to_h_size = alpha.exp(E::PositiveInteger::from(
            self.options.size_subgroup_h as u64,
        ));
        let mut u_numerator = vec![E::ZERO; (self.options.size_subgroup_h).try_into().unwrap()];
        u_numerator[0] = alpha_to_h_size.neg();
        u_numerator.push(E::ONE);
        let u_denominator = vec![alpha.neg(), E::ONE];
        let mut u_alpha = polynom::div(&u_numerator, &u_denominator);

        let mut prod = Vec::<E>::new();
        let eval_twiddles = fft::get_twiddles(self.options.evaluation_domain.len());
        let mut f_1_eval = self
            .f_1_poly_coeffs
            .iter()
            .map(|i| E::from(*i))
            .collect::<Vec<E>>();
        fractal_utils::polynomial_utils::pad_with_zeroes(
            &mut f_1_eval,
            self.options.evaluation_domain.len(),
        );

        fft::evaluate_poly(&mut f_1_eval, &mut eval_twiddles.clone());
        let mut f_2_eval = self
            .f_2_poly_coeffs
            .iter()
            .map(|i| E::from(*i))
            .collect::<Vec<E>>();
        fractal_utils::polynomial_utils::pad_with_zeroes(
            &mut f_2_eval,
            self.options.evaluation_domain.len(),
        );
        fft::evaluate_poly(&mut f_2_eval, &mut eval_twiddles.clone());
        fractal_utils::polynomial_utils::pad_with_zeroes(
            &mut u_alpha,
            self.options.evaluation_domain.len(),
        );
        fft::evaluate_poly(&mut u_alpha, &mut eval_twiddles.clone());
        //none of that fft nonsense, let's do this the lagrange way
        //f_1_eval = polynom::eval_many(&self.f_1_poly_coeffs, &self.options.evaluation_domain).iter().map(|i| E::from(*i)).collect::<Vec<E>>();
        //f_2_eval = polynom::eval_many(&self.f_2_poly_coeffs, &self.options.evaluation_domain).iter().map(|i| E::from(*i)).collect::<Vec<E>>();
        // u_alpha = polynom::eval_many(&u_alpha, &self.evaluation_domain_e);
        for pos in 0..self.options.evaluation_domain.len() {
            let next = (u_alpha[pos] * f_1_eval[pos]) - (t_alpha[pos] * f_2_eval[pos]);
            prod.push(next);
        }
        prod
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
    ) -> Result<(), ProverError> {
        match self.get_current_layer() {
            0 => {
                self.lincheck_layer_one(query, accumulator);
            }
            1 => {
                self.lincheck_layer_two(query, accumulator);
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

    fn get_fractal_options(&self) -> FractalOptions<B> {
        self.options.clone()
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredProver<B, E, H, LayeredLincheckProof<B, E>> for LincheckProver<B, E, H>
{
    fn generate_proof(
        &mut self,
        prover_key: &Option<ProverKey<B,E,H>>,
        public_inputs_bytes: Vec<u8>,
    ) -> Result<TopLevelProof<B, E, H>, ProverError> {
        let options = self.get_fractal_options();
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
            public_inputs_bytes
        );
        //let mut layer_commitments = [<H as Hasher>::hash(&[0u8]); 3];
        let mut layer_commitments = vec![];
        let mut local_queries = Vec::<E>::new();

        // BEGIN JANK
        /*let inv_twiddles_h = fft::get_inv_twiddles(self.variable_assignment.len());
        let mut z_coeffs = &mut self.variable_assignment.clone(); // evals
        fft::interpolate_poly_with_offset(
            &mut z_coeffs,
            &inv_twiddles_h,
            prover_key.as_ref().unwrap().params.eta,
        ); // coeffs

        let f_az_coeffs = compute_matrix_mul_poly_coeffs(
            &prover_key.as_ref().unwrap().matrix_a_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            prover_key.as_ref().unwrap().params.eta,
        )?;

        let f_bz_coeffs = compute_matrix_mul_poly_coeffs(
            &prover_key.as_ref().unwrap().matrix_b_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            prover_key.as_ref().as_ref().unwrap().params.eta,
        )?;

        let f_cz_coeffs = compute_matrix_mul_poly_coeffs(
            &prover_key.as_ref().unwrap().matrix_c_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            prover_key.as_ref().unwrap().params.eta,
        )?;

        //TODO: Put in correct degree constraints
        acc.add_unchecked_polynomial(z_coeffs.to_vec());
        acc.add_unchecked_polynomial(f_az_coeffs.to_vec());
        acc.add_unchecked_polynomial(f_bz_coeffs.to_vec());
        acc.add_unchecked_polynomial(f_cz_coeffs.to_vec());*/

        // END JANK

        for i in 0..self.get_num_layers() {
            // local_queries.push(query);
            // Doing this rn to make sure prover and verifier sample identically
            if i > 0 {
                let previous_commit = acc.get_layer_commitment(i)?;
                channel.commit_fractal_iop_layer(previous_commit);
                coin.reseed(previous_commit);
            }
            let query = coin.draw().expect("failed to draw FRI alpha"); //channel.draw_fri_alpha();
            local_queries.push(query);
            self.run_next_layer(query, &mut acc)?;
            layer_commitments.push(acc.commit_layer()?); //todo: do something with this
        }

        let queries = acc.draw_query_positions()?;

        let beta = local_queries[1];

        let layer_decommits = vec![
            acc.decommit_layer_with_queries(1, &queries)?,
            acc.decommit_layer_with_queries(2, &queries)?,
        ];
        let gammas = vec![
            self.retrieve_gamma(beta)?,
        ];

        let preprocessing_decommitment = prover_key.as_ref().unwrap().accumulator.decommit_layer_with_queries(1, &queries)?;
        
        let low_degree_proof = acc.create_fri_proof()?;

        let proof = TopLevelProof {
            preprocessing_decommitment,
            layer_commitments: layer_commitments.to_vec(),
            layer_decommitments: layer_decommits,
            unverified_misc: gammas,
            low_degree_proof
        };
        Ok(proof)
    }
}

fn compute_matrix_mul_poly_coeffs<B: StarkField>(
    matrix: &Matrix<B>,
    vec: &Vec<B>,
    inv_twiddles: &[B],
    eta: B,
) -> Result<Vec<B>, ProverError> {
    let mut product = matrix.dot(vec); // as evals
    fft::interpolate_poly_with_offset(&mut product, inv_twiddles, eta); // as coeffs
    Ok(product) // as coeffs
}