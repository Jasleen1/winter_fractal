use std::marker::PhantomData;

use fractal_indexer::snark_keys::*;
use fractal_proofs::{fft, polynom, FractalProof, LincheckProof, TryInto};
use models::r1cs::Matrix;

use winter_crypto::{ElementHasher, RandomCoin};
use winter_math::{FieldElement, StarkField};

use crate::{
    errors::ProverError,
    lincheck_prover::LincheckProver,
    rowcheck_prover::RowcheckProver,
    FractalOptions,
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
    public_coin: RandomCoin<B, H>,
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
        pub_inputs_bytes: Vec<u8>,
    ) -> Self {
        let coin_seed = pub_inputs_bytes;
        FractalProver {
            prover_key,
            options,
            witness,
            variable_assignment,
            public_coin: RandomCoin::new(&coin_seed),
            _e: PhantomData,
        }
    }

    pub fn generate_proof(&mut self) -> Result<FractalProof<B, E, H>, ProverError> {
        // This is the less efficient version and assumes only dealing with the var assignment,
        // not z = (x, w)
        let alpha = self.public_coin.draw().expect("failed to draw OOD point");
        let inv_twiddles_h = fft::get_inv_twiddles(self.variable_assignment.len());

        // 1. Generate lincheck proofs for the A,B,C matrices.
        let mut z_coeffs = &mut self.variable_assignment.clone();  // evals
        fft::interpolate_poly_with_offset(&mut z_coeffs, &inv_twiddles_h, self.prover_key.params.eta);  // coeffs
        let f_az_coeffs = self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.matrix_a_index.matrix, 
            &self.variable_assignment.clone(), 
            &inv_twiddles_h,
            self.prover_key.params.eta)?;
        let lincheck_a = self.create_lincheck_proof(
            alpha,
            &self.prover_key.matrix_a_index,
            &z_coeffs.clone(),
            &f_az_coeffs)?;

        let f_bz_coeffs = self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.matrix_b_index.matrix, 
            &self.variable_assignment.clone(), 
            &inv_twiddles_h,
            self.prover_key.params.eta)?;
        let lincheck_b = self.create_lincheck_proof(
            alpha,
            &self.prover_key.matrix_b_index,
            &z_coeffs.clone(),
            &f_bz_coeffs)?;

        let f_cz_coeffs = self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.matrix_c_index.matrix, 
            &self.variable_assignment.clone(), 
            &inv_twiddles_h,
            self.prover_key.params.eta)?;
        let lincheck_c = self.create_lincheck_proof(
            alpha,
            &self.prover_key.matrix_c_index,
            &z_coeffs.clone(),
            &f_cz_coeffs)?;
        
        println!("Done with linchecks");
        
        // 2. Generate the rowcheck proof.

        // Evaluate the Az, Bz, Cz polynomials.
        // let eval_twiddles = fft::get_twiddles(self.options.evaluation_domain.len());

        // let mut f_az_evals = f_az_coeffs.clone();
        let f_az_evals = polynom::eval_many(&f_az_coeffs.clone(), &self.options.evaluation_domain);
        // fft::evaluate_poly(&mut f_az_evals, &eval_twiddles);

        let f_bz_evals = polynom::eval_many(&f_bz_coeffs.clone(), &self.options.evaluation_domain);
        // fft::evaluate_poly(&mut f_bz_evals, &eval_twiddles);

        let f_cz_evals = polynom::eval_many(&f_cz_coeffs.clone(), &self.options.evaluation_domain);
        // fft::evaluate_poly(&mut f_cz_evals, &eval_twiddles);
        
        // Issue a rowcheck proof.
        let rowcheck_prover = RowcheckProver::<B, E, H>::new(
            f_az_coeffs,
            f_bz_coeffs,
            f_cz_coeffs,
            self.options.degree_fs,
            self.options.size_subgroup_h.try_into().unwrap(),
            self.options.evaluation_domain.clone(),
            self.options.fri_options.clone(),
            self.options.num_queries,
            self.prover_key.params.max_degree,
            self.prover_key.params.eta,
        );
        let rowcheck_proof = rowcheck_prover.generate_proof()?;
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
        let mut product = matrix.dot(vec);  // as evals
        fft::interpolate_poly_with_offset(&mut product, inv_twiddles, eta);  // as coeffs
        Ok(product)  // as coeffs
    }

    // Indexed matrix; variable assignments as polynomial evaluation points.
    fn create_lincheck_proof(
        &self,
        alpha: B,
        matrix_index: &ProverMatrixIndex<H, B>,
        z_coeffs: &Vec<B>,
        prod_m_z_coeffs: &Vec<B>) -> Result<LincheckProof<B, E, H>, ProverError> {

        let lincheck_prover = LincheckProver::<B, E, H>::new(
            alpha,
            &matrix_index,
            prod_m_z_coeffs.to_vec(),
            z_coeffs.to_vec(),
            &self.options,
        );
        let lincheck_proof = lincheck_prover.generate_lincheck_proof()?;
        Ok(lincheck_proof)
        // let mut matrix = Matrix::new(matrix_label, Vec::<Vec<B>>::new())?;
        // match matrix_label {
        //     "a" => {
        //         matrix = self.prover_key.matrix_a_index.matrix.clone();
        //     }
        //     "b" => {
        //         matrix = self.prover_key.matrix_b_index.matrix.clone();
        //     }
        //     "c" => {
        //         matrix = self.prover_key.matrix_c_index.matrix.clone();
        //     }
        //     _ => {}
        // }
        // if matrix.mat.len() == 0 {
        //     return Err(ProverError::InvalidMatrixName(matrix_label.to_string()));
        // }
        // let mut f_2_vals = matrix.dot(self.variable_assignment.clone());
        // println!(
        //     "Matrix rows = {}, matrix cols = {}",
        //     matrix.dims.0, matrix.dims.1
        // );
        // println!("Len f_2 = {}", f_2_vals.len());
        // println!("Len twiddles = {}", inv_twiddles.len());
        // fft::interpolate_poly(&mut f_2_vals, inv_twiddles);
        // Ok(f_2_vals)
    }
}
