use std::marker::PhantomData;

use fractal_indexer::{index::IndexParams, snark_keys::*};
use fractal_proofs::{
    fft, polynom, DefaultProverChannel, FractalProof, FriOptions, InitialPolyProof,
    LayeredFractalProof, LincheckProof, LowDegreeBatchProof, MultiEval, MultiPoly, TryInto,
};
use models::r1cs::Matrix;

use winter_crypto::{ElementHasher, Hasher, MerkleTree, RandomCoin};
use winter_fri::ProverChannel;
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use crate::{
    accumulator::Accumulator, channel::DefaultFractalProverChannel, errors::ProverError,
    lincheck_prover::LincheckProver, rowcheck_prover::RowcheckProver, FractalOptions,
    LayeredProver, LayeredSubProver, FRACTAL_LAYERS,
};

pub struct FractalProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    prover_key: Option<ProverKey<B, E, H>>,
    options: FractalOptions<B>,
    witness: Vec<B>,
    variable_assignment: Vec<B>,
    pub_input_bytes: Vec<u8>,
    _e: PhantomData<E>,
    current_layer: usize,
    // state variables
    f_az_coeffs: Vec<B>,
    f_bz_coeffs: Vec<B>,
    f_cz_coeffs: Vec<B>,
    z_coeffs: Vec<B>,
    lincheck_provers: Vec<LincheckProver<B, E, H>>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > FractalProver<B, E, H>
{
    pub fn new(
        prover_key: ProverKey<B, E, H>,
        options: FractalOptions<B>,
        witness: Vec<B>,
        variable_assignment: Vec<B>,
        pub_input_bytes: Vec<u8>,
    ) -> Self {
        FractalProver {
            prover_key: Some(prover_key),
            options,
            witness,
            variable_assignment,
            pub_input_bytes,
            _e: PhantomData,
            current_layer: 0,
            f_az_coeffs: Vec::new(),
            f_bz_coeffs: Vec::new(),
            f_cz_coeffs: Vec::new(),
            z_coeffs: Vec::new(),
            lincheck_provers: Vec::new(),
        }
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

    fn fractal_layer_one(
        &mut self,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError> {
        let inv_twiddles_h = fft::get_inv_twiddles(self.variable_assignment.len());
        // 1. Generate lincheck proofs for the A,B,C matrices.
        let mut z_coeffs = &mut self.variable_assignment.clone(); // evals
        fft::interpolate_poly_with_offset(
            &mut z_coeffs,
            &inv_twiddles_h,
            self.prover_key.as_ref().unwrap().params.eta,
        ); // coeffs
        let f_az_coeffs = &mut self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.as_ref().unwrap().matrix_a_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            self.prover_key.as_ref().unwrap().params.eta,
        )?;

        let f_bz_coeffs = &mut self.compute_matrix_mul_poly_coeffs(
            &self
                .prover_key
                .as_ref()
                .as_ref()
                .unwrap()
                .matrix_b_index
                .matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            self.prover_key.as_ref().as_ref().unwrap().params.eta,
        )?;

        let f_cz_coeffs = &mut self.compute_matrix_mul_poly_coeffs(
            &self.prover_key.as_ref().unwrap().matrix_c_index.matrix,
            &self.variable_assignment.clone(),
            &inv_twiddles_h,
            self.prover_key.as_ref().unwrap().params.eta,
        )?;

        self.f_az_coeffs = f_az_coeffs.to_vec();
        self.f_bz_coeffs = f_bz_coeffs.to_vec();
        self.f_cz_coeffs = f_cz_coeffs.to_vec();
        self.z_coeffs = z_coeffs.to_vec();

        //TODO: Put in correct degree constraints
        accumulator.add_unchecked_polynomial(z_coeffs.to_vec());
        accumulator.add_unchecked_polynomial(f_az_coeffs.to_vec());
        accumulator.add_unchecked_polynomial(f_bz_coeffs.to_vec());
        accumulator.add_unchecked_polynomial(f_cz_coeffs.to_vec());
        Ok(())
    }

    fn fractal_layer_two(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError> {
        // 1. Generate the rowcheck proof.
        // Evaluate the Az, Bz, Cz polynomials.
        let mut rowcheck_prover = RowcheckProver::<B, E, H>::new(
            self.f_az_coeffs.clone(),
            self.f_bz_coeffs.clone(),
            self.f_cz_coeffs.clone(),
            self.options.clone(),
        );

        //hacky way to avoid lifetimes: move prover_key contents to LincheckProvers in this step
        let prover_key = std::mem::replace(&mut self.prover_key, None).unwrap();
        self.prover_key = None;
        let a_index = prover_key.matrix_a_index;
        let b_index = prover_key.matrix_b_index;
        let c_index = prover_key.matrix_c_index;

        let mut lincheck_prover_a = LincheckProver::<B, E, H>::new(
            a_index,
            self.f_az_coeffs.to_vec(),
            self.z_coeffs.to_vec(),
            &self.options,
        );
        let mut lincheck_prover_b = LincheckProver::<B, E, H>::new(
            b_index,
            self.f_bz_coeffs.to_vec(),
            self.z_coeffs.to_vec(),
            &self.options,
        );
        let mut lincheck_prover_c = LincheckProver::<B, E, H>::new(
            c_index,
            self.f_cz_coeffs.to_vec(),
            self.z_coeffs.to_vec(),
            &self.options,
        );

        rowcheck_prover.run_next_layer(query, accumulator)?;
        lincheck_prover_a.run_next_layer(query, accumulator)?;
        lincheck_prover_b.run_next_layer(query, accumulator)?;
        lincheck_prover_c.run_next_layer(query, accumulator)?;
        self.lincheck_provers = vec![lincheck_prover_a, lincheck_prover_b, lincheck_prover_c];

        Ok(())
    }

    fn fractal_layer_three(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError> {
        for lincheck_prover in self.lincheck_provers.iter_mut() {
            lincheck_prover.run_next_layer(query, accumulator)?;
        }
        Ok(())
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredSubProver<B, E, H> for FractalProver<B, E, H>
{
    fn run_next_layer(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError> {
        match self.current_layer {
            0 => {
                self.fractal_layer_one(accumulator)?;
                self.current_layer += 1;
            }
            1 => {
                self.fractal_layer_two(query, accumulator)?;
                self.current_layer += 1;
            }
            2 => {
                self.fractal_layer_three(query, accumulator)?;
                self.current_layer += 1;
            }
            _ => (),
        };
        Ok(())
    }
    fn get_current_layer(&self) -> usize {
        self.current_layer
    }
    fn get_num_layers(&self) -> usize {
        FRACTAL_LAYERS
    }
    fn get_fractal_options(&self) -> FractalOptions<B> {
        self.options.clone()
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredProver<B, E, H, LayeredFractalProof<B, E, H>> for FractalProver<B, E, H>
{
    fn generate_proof(
        &mut self,
        public_input_bytes: Vec<u8>,
    ) -> Result<LayeredFractalProof<B, E, H>, ProverError> {
        let options = self.get_fractal_options();
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            options.evaluation_domain.len(),
            options.num_queries,
            public_input_bytes,
        );
        let mut acc = Accumulator::<B, E, H>::new(
            options.evaluation_domain.len(),
            B::ONE,
            options.evaluation_domain.clone(),
            options.num_queries,
            options.fri_options.clone(),
        );
        let mut layer_commitments = [<H as Hasher>::hash(&[0u8]); 3];
        let mut local_queries = Vec::<E>::new();

        let mut coin = RandomCoin::<B, H>::new(&vec![]);

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
            layer_commitments[i] = acc.commit_layer()?; //todo: do something with this
        }

        let queries = acc.draw_query_positions()?;

        let beta = local_queries[2];

        let preprocessing_decommits_a =
            self.lincheck_provers[0].decommit_proprocessing(&queries)?;
        let preprocessing_decommits_b =
            self.lincheck_provers[1].decommit_proprocessing(&queries)?;
        let preprocessing_decommits_c =
            self.lincheck_provers[2].decommit_proprocessing(&queries)?;
        let layer_decommits = [
            acc.decommit_layer_with_queries(1, &queries)?,
            acc.decommit_layer_with_queries(2, &queries)?,
            acc.decommit_layer_with_queries(3, &queries)?,
        ];
        let gammas = [
            self.lincheck_provers[0].retrieve_gamma(beta)?,
            self.lincheck_provers[1].retrieve_gamma(beta)?,
            self.lincheck_provers[2].retrieve_gamma(beta)?,
        ];

        let low_degree_proof = acc.create_fri_proof()?;
        let proof = LayeredFractalProof {
            preprocessing_decommits_a,
            preprocessing_decommits_b,
            preprocessing_decommits_c,
            layer_commitments,
            gammas,
            layer_decommits,
            low_degree_proof,
        };
        Ok(proof)
    }
}
