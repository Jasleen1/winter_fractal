use std::marker::PhantomData;

use fractal_indexer::{index::IndexParams, snark_keys::*};
use fractal_proofs::{
    fft, polynom, DefaultProverChannel, FractalProof, FriOptions, InitialPolyProof,
    LayeredFractalProof, LincheckProof, LowDegreeBatchProof, MultiEval, MultiPoly, TryInto, LayeredLincheckProof, LayeredRowcheckProof, TopLevelProof, IopData,
};
use models::r1cs::Matrix;

use winter_crypto::{ElementHasher, Hasher, MerkleTree, RandomCoin, BatchMerkleProof};
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
            &self.prover_key.as_ref().unwrap().matrix_b_index.matrix,
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
    > LayeredProver<B, E, H, LayeredFractalProof<B, E>> for FractalProver<B, E, H>
{
    fn generate_proof(
        &mut self,
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
        let mut layer_commitments = [<H as Hasher>::hash(&[0u8]); 3];
        let mut local_queries = Vec::<E>::new();

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
        let layer_decommits = vec![
            acc.decommit_layer_with_queries(1, &queries)?,
            acc.decommit_layer_with_queries(2, &queries)?,
            acc.decommit_layer_with_queries(3, &queries)?,
        ];
        let gammas = vec![
            self.lincheck_provers[0].retrieve_gamma(beta)?,
            self.lincheck_provers[1].retrieve_gamma(beta)?,
            self.lincheck_provers[2].retrieve_gamma(beta)?,
        ];

        let preprocessing_decommitments = [preprocessing_decommits_a, preprocessing_decommits_b, preprocessing_decommits_c];
        let low_degree_proof = acc.create_fri_proof()?;
        /*let subcomponents = parse_into_subcomponents(
            preprocessing_decommits_a,
            preprocessing_decommits_b,
            preprocessing_decommits_c,
            layer_commitments,
            gammas,
            layer_decommits,
        );
        let [lincheck_a, lincheck_b, lincheck_c]  = subcomponents.0;
        let iop_data = LayeredFractalProof{
            rowcheck: subcomponents.1,
            lincheck_a,
            lincheck_b,
            lincheck_c
        };*/
        let proof = TopLevelProof {
            preprocessing_decommitments,
            layer_commitments: layer_commitments.to_vec(),
            layer_decommitments: layer_decommits,
            unverified_misc: gammas,
            low_degree_proof
        };
        Ok(proof)
    }
}

/*pub struct LayeredFractalProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub preprocessing_decommits_a: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub preprocessing_decommits_b: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub preprocessing_decommits_c: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub layer_commitments: [H::Digest; 3],
    pub gammas: [E; 3],
    pub layer_decommits: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub low_degree_proof: LowDegreeBatchProof<B, E, H>,
}*/
/* 
fn parse_into_subcomponents<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    preprocessing_decommits_a: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    preprocessing_decommits_b: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    preprocessing_decommits_c: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    layer_commitments: [H::Digest; 3],
    gammas: [E; 3],
    layer_decommits: Vec<(Vec<Vec<E>>, BatchMerkleProof<H>)>,
) -> ([LayeredLincheckProof<B, E>; 3], LayeredRowcheckProof<B, E>) {
    // We also need a new error type for this function that can be passed through to the VerifierError type

    // Matrix A preprocessing
    let row_a = extract_vec_e(&preprocessing_decommits_a[0].0, 0);
    let col_a = extract_vec_e(&preprocessing_decommits_a[1].0, 0);
    let val_a = extract_vec_e(&preprocessing_decommits_a[2].0, 0);

    // Matrix B preprocessing
    let row_b = extract_vec_e(&preprocessing_decommits_b[0].0, 0);
    let col_b = extract_vec_e(&preprocessing_decommits_b[1].0, 0);
    let val_b = extract_vec_e(&preprocessing_decommits_b[2].0, 0);

    // Matrix C preprocessing
    let row_c = extract_vec_e(&preprocessing_decommits_c[0].0, 0);
    let col_c = extract_vec_e(&preprocessing_decommits_c[1].0, 0);
    let val_c = extract_vec_e(&preprocessing_decommits_c[2].0, 0);

    // get values from the first layer
    let f_z_vals = extract_vec_e(&layer_decommits[0].0, 0);
    let f_az_vals = extract_vec_e(&layer_decommits[0].0, 1);
    let f_bz_vals = extract_vec_e(&layer_decommits[0].0, 2);
    let f_cz_vals = extract_vec_e(&layer_decommits[0].0, 3);

    // get values from the second layer
    let s_vals = extract_vec_e(&layer_decommits[1].0, 0);
    let t_alpha_a_vals = extract_vec_e(&layer_decommits[1].0, 1);
    let product_sumcheck_a_vals = extract_sumcheck_vec_e(&layer_decommits[1].0, 2, 3);
    let t_alpha_b_vals = extract_vec_e(&layer_decommits[1].0, 4);
    let product_sumcheck_b_vals = extract_sumcheck_vec_e(&layer_decommits[1].0, 5, 6);
    let t_alpha_c_vals = extract_vec_e(&layer_decommits[1].0, 7);
    let product_sumcheck_c_vals = extract_sumcheck_vec_e(&layer_decommits[1].0, 8, 9);

    // get values from the third layer
    let matrix_sumcheck_a_vals = extract_sumcheck_vec_e(&layer_decommits[2].0, 0, 1);
    let matrix_sumcheck_b_vals = extract_sumcheck_vec_e(&layer_decommits[2].0, 2, 3);
    let matrix_sumcheck_c_vals = extract_sumcheck_vec_e(&layer_decommits[2].0, 4, 5);

    // Get the lincheck query values appropriately
    let mut coin = RandomCoin::<B, H>::new(&vec![]);
    coin.reseed(layer_commitments[0]);
    let alpha: E = coin.draw().expect("failed to draw FRI alpha");

    coin.reseed(layer_commitments[1]);
    let beta: E = coin.draw().expect("failed to draw FRI alpha");

    let lincheck_a_proof = LayeredLincheckProof {
        row_vals: row_a,
        col_vals: col_a,
        val_vals: val_a,
        f_z_vals: f_z_vals.clone(),
        f_mz_vals: f_az_vals.clone(),
        t_alpha_vals: t_alpha_a_vals,
        product_sumcheck_vals: product_sumcheck_a_vals,
        matrix_sumcheck_vals: matrix_sumcheck_a_vals,
        alpha,
        beta,
        gamma: gammas[0],
    };

    let lincheck_b_proof = LayeredLincheckProof {
        row_vals: row_b,
        col_vals: col_b,
        val_vals: val_b,
        f_z_vals: f_z_vals.clone(),
        f_mz_vals: f_bz_vals.clone(),
        t_alpha_vals: t_alpha_b_vals,
        product_sumcheck_vals: product_sumcheck_b_vals,
        matrix_sumcheck_vals: matrix_sumcheck_b_vals,
        alpha,
        beta,
        gamma: gammas[1],
    };

    let lincheck_c_proof = LayeredLincheckProof {
        row_vals: row_c,
        col_vals: col_c,
        val_vals: val_c,
        f_z_vals: f_z_vals.clone(),
        f_mz_vals: f_cz_vals.clone(),
        t_alpha_vals: t_alpha_c_vals,
        product_sumcheck_vals: product_sumcheck_c_vals,
        matrix_sumcheck_vals: matrix_sumcheck_c_vals,
        alpha,
        beta,
        gamma: gammas[2],
    };

    let rowcheck_proof = LayeredRowcheckProof {
        f_z_vals,
        f_az_vals,
        f_bz_vals,
        f_cz_vals,
        s_vals,
    };

    (
        [lincheck_a_proof, lincheck_b_proof, lincheck_c_proof],
        rowcheck_proof,
    )
}

fn extract_vec_e<B: StarkField, E: FieldElement<BaseField = B>>(
    vec_of_decommits: &Vec<Vec<E>>,
    position: usize,
) -> Vec<E> {
    vec_of_decommits
        .iter()
        .map(|x| x[position])
        .collect::<Vec<E>>()
}

fn extract_sumcheck_vec_e<B: StarkField, E: FieldElement<BaseField = B>>(
    vec_of_decommits: &Vec<Vec<E>>,
    position_g: usize,
    position_e: usize,
) -> Vec<(E, E)> {
    vec_of_decommits
        .iter()
        .map(|x| (x[position_g], x[position_e]))
        .collect::<Vec<(E, E)>>()
}
*/