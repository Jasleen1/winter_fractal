use std::marker::PhantomData;

use fractal_indexer::snark_keys::*;
use fractal_proofs::{
    fft, polynom, DefaultProverChannel, FractalProof, LincheckProof, MultiEval, TryInto,
};
use models::r1cs::Matrix;

use winter_crypto::{ElementHasher, MerkleTree, RandomCoin};
use winter_fri::ProverChannel;
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use crate::{
    channel::DefaultFractalProverChannel, errors::ProverError, lincheck_prover::LincheckProver,
    rowcheck_prover::RowcheckProver, FractalOptions,
};

pub struct FractalProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    prover_key: ProverKey<H, B>,
    options: FractalOptions<B>,
    witness: Vec<B>,
    variable_assignment: Vec<B>,
    pub_input_bytes: Vec<u8>,
    _e: PhantomData<E>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > FractalProver<B, E, H>
{
    pub fn new(
        prover_key: ProverKey<H, B>,
        options: FractalOptions<B>,
        witness: Vec<B>,
        variable_assignment: Vec<B>,
        pub_input_bytes: Vec<u8>,
    ) -> Self {
        FractalProver {
            prover_key,
            options,
            witness,
            variable_assignment,
            pub_input_bytes,
            _e: PhantomData,
        }
    }

    pub fn generate_proof(&mut self) -> Result<FractalProof<B, E, H>, ProverError> {
        let mut channel = &mut DefaultFractalProverChannel::<B, E, H>::new(
            self.options.evaluation_domain.len(),
            self.options.num_queries,
            self.pub_input_bytes.clone(),
        );

        // This is the less efficient version and assumes only dealing with the var assignment,
        // not z = (x, w)
        // TODO the function used here needs to be modified.
        let alpha = channel.draw_random_b_pt();
        let inv_twiddles_h = fft::get_inv_twiddles(self.variable_assignment.len());

        // 1. Generate lincheck proofs for the A,B,C matrices.
        let mut z_coeffs = &mut self.variable_assignment.clone(); // evals
        fft::interpolate_poly_with_offset(
            &mut z_coeffs,
            &inv_twiddles_h,
            self.prover_key.params.eta,
        ); // coeffs
        let mut f_az_coeffs = &mut self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.matrix_a_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            self.prover_key.params.eta,
        )?;

        let mut f_bz_coeffs = &mut self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.matrix_b_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            self.prover_key.params.eta,
        )?;

        let mut f_cz_coeffs = &mut self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.matrix_c_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            self.prover_key.params.eta,
        )?;

        let coefficients = vec![
            z_coeffs.clone(),
            f_az_coeffs.clone(),
            f_bz_coeffs.clone(),
            f_cz_coeffs.clone(),
        ];
        let mut initial_vector_polys =
            MultiEval::<B, E, H>::new(coefficients, self.options.evaluation_domain.len(), B::ONE);

        // let eval_twiddles = fft::get_twiddles(self.options.evaluation_domain.len());
        // let mut f_z_eval = z_coeffs.clone();
        // fractal_utils::polynomial_utils::pad_with_zeroes(
        //     &mut f_z_eval,
        //     self.options.evaluation_domain.len(),
        // );

        // fft::evaluate_poly(&mut f_z_eval, &mut eval_twiddles.clone());

        // let mut f_az_eval = f_az_coeffs.clone();
        // fractal_utils::polynomial_utils::pad_with_zeroes(
        //     &mut f_az_eval,
        //     self.options.evaluation_domain.len(),
        // );

        // fft::evaluate_poly(&mut f_az_eval, &mut eval_twiddles.clone());

        // let mut f_bz_eval = f_bz_coeffs.clone();
        // fractal_utils::polynomial_utils::pad_with_zeroes(
        //     &mut f_bz_eval,
        //     self.options.evaluation_domain.len(),
        // );

        // fft::evaluate_poly(&mut f_bz_eval, &mut eval_twiddles.clone());

        // let mut f_cz_eval = f_cz_coeffs.clone();
        // fractal_utils::polynomial_utils::pad_with_zeroes(
        //     &mut f_cz_eval,
        //     self.options.evaluation_domain.len(),
        // );

        // fft::evaluate_poly(&mut f_cz_eval, &mut eval_twiddles.clone());

        // let transposed_evaluations = transpose_slice(evaluations);
        // let hashed_evaluations = hash_values::<H, E, N>(&transposed_evaluations);
        // let evaluation_tree =
        //     MerkleTree::<H>::new(hashed_evaluations).expect("failed to construct FRI layer tree");

        let lincheck_a = self.create_lincheck_proof(
            alpha,
            &self.prover_key.matrix_a_index,
            &z_coeffs.clone(),
            &f_az_coeffs,
            channel,
        )?;

        let lincheck_b = self.create_lincheck_proof(
            alpha,
            &self.prover_key.matrix_b_index,
            &z_coeffs.clone(),
            &f_bz_coeffs,
            channel,
        )?;

        let lincheck_c = self.create_lincheck_proof(
            alpha,
            &self.prover_key.matrix_c_index,
            &z_coeffs.clone(),
            &f_cz_coeffs,
            channel,
        )?;

        println!("Done with linchecks");

        // 2. Generate the rowcheck proof.

        // Evaluate the Az, Bz, Cz polynomials.
        // let eval_twiddles = fft::get_twiddles(self.options.evaluation_domain.len());

        // let mut f_az_evals = f_az_coeffs.clone();
        // let f_az_evals = polynom::eval_many(&f_az_coeffs.clone(), &self.options.evaluation_domain);
        // fft::evaluate_poly(&mut f_az_evals, &eval_twiddles);

        // let f_bz_evals = polynom::eval_many(&f_bz_coeffs.clone(), &self.options.evaluation_domain);
        // fft::evaluate_poly(&mut f_bz_evals, &eval_twiddles);

        // let f_cz_evals = polynom::eval_many(&f_cz_coeffs.clone(), &self.options.evaluation_domain);
        // fft::evaluate_poly(&mut f_cz_evals, &eval_twiddles);

        // Issue a rowcheck proof.
        let rowcheck_prover = RowcheckProver::<B, E, H>::new(
            f_az_coeffs.clone(),
            f_bz_coeffs.clone(),
            f_cz_coeffs.clone(),
            self.options.degree_fs,
            self.options.size_subgroup_h.try_into().unwrap(),
            self.options.evaluation_domain.clone(),
            self.options.fri_options.clone(),
            self.options.num_queries,
            self.prover_key.params.max_degree,
            self.prover_key.params.eta,
        );
        let rowcheck_proof = rowcheck_prover.generate_proof(channel)?;
        println!("Done with rowcheck");
        // 3. Build and return an overall fractal proof.
        Ok(FractalProof {
            rowcheck_proof,
            lincheck_a,
            lincheck_b,
            lincheck_c,
        })
    }

    // Multiply a matrix times a vector of evaluations, then interpolate a poly and return its coeffs.
    fn compute_matrix_mul_poly_coeffs(
        &self,
        matrix: &Matrix<B>,
        vec: &Vec<B>,
        inv_twiddles: &[B],
        eta: B,
    ) -> Result<Vec<B>, ProverError> {
        let mut product = matrix.dot(vec); // as evals
        fft::interpolate_poly_with_offset(&mut product, inv_twiddles, eta); // as coeffs
        Ok(product) // as coeffs
    }

    // Indexed matrix; variable assignments as polynomial evaluation points.
    fn create_lincheck_proof(
        &self,
        alpha: B,
        matrix_index: &ProverMatrixIndex<H, B>,
        z_coeffs: &Vec<B>,
        prod_m_z_coeffs: &Vec<B>,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
    ) -> Result<LincheckProof<B, E, H>, ProverError> {
        let lincheck_prover = LincheckProver::<B, E, H>::new(
            alpha,
            &matrix_index,
            prod_m_z_coeffs.to_vec(),
            z_coeffs.to_vec(),
            &self.options,
        );
        let lincheck_proof = lincheck_prover.generate_lincheck_proof(channel)?;
        Ok(lincheck_proof)
    }
}
