//use core::num::dec2flt::parse;

use crate::{
    batched_lincheck_verifier::{
        self, parse_proofs_for_subroutines as other_parse_proofs, verify_layered_lincheck_proof,
    },
    errors::FractalVerifierError, //lincheck_verifier::verify_layered_lincheck_proof,
    rowcheck_verifier::verify_layered_rowcheck_proof,
    verifier,
};

use fractal_accumulator_verifier::accumulator_verifier::AccumulatorVerifier;
use fractal_indexer::snark_keys::*;
use fractal_proofs::{
    BatchedLayeredLincheckProof, FieldElement, FractalProof, IopData, LayeredFractalProof,
    LayeredLincheckProof, LayeredRowcheckProof, MultiEval, MultiPoly, StarkField, TopLevelProof,
};

use fractal_prover::batched_lincheck_full_prover;
use fractal_utils::channel::DefaultFractalProverChannel;
use fractal_utils::FractalOptions;
use log::debug;
use winter_crypto::{ElementHasher, RandomCoin};

use crate::{lincheck_verifier::verify_lincheck_proof, rowcheck_verifier::verify_rowcheck_proof};

/*pub fn verify_fractal_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: VerifierKey<B, E, H>,
    proof: FractalProof<B, E, H>,
    pub_inputs_bytes: Vec<u8>,
    options: FractalOptions<B>,
) -> Result<(), FractalVerifierError> {
    let mut public_coin = RandomCoin::<_, H>::new(&pub_inputs_bytes);
    let expected_alpha: B = public_coin.draw().expect("failed to draw OOD point");
    // let mut channel = DefaultFractalProverChannel::new();

    debug!(
        "Lincheck a indexes: {:?}",
        &proof.lincheck_a.products_sumcheck_proof.queried_positions
    );

    let initial_evals = proof.initial_poly_proof.evals.clone();
    public_coin.reseed(proof.initial_poly_proof.commitment);
    //let indices = public_coin.draw_integers(options.num_queries, options.evaluation_domain.len()).unwrap();
    let indices = proof.rowcheck_proof.s_proof.queried_positions.clone();

    //todo: add this back
    /*MultiEval::<B, E, H>::batch_verify_values_and_proofs_at(proof.initial_poly_proof.evals,
    &proof.initial_poly_proof.commitment,
    &proof.initial_poly_proof.proof, indices.clone())?;*/

    verify_lincheck_proof(
        &verifier_key,
        proof.lincheck_a,
        expected_alpha,
        &mut public_coin,
        options.num_queries,
    )?;
    println!("Lincheck a verified");
    verify_lincheck_proof(
        &verifier_key,
        proof.lincheck_b,
        expected_alpha,
        &mut public_coin,
        options.num_queries,
    )?;
    println!("Lincheck b verified");
    verify_lincheck_proof(
        &verifier_key,
        proof.lincheck_c,
        expected_alpha,
        &mut public_coin,
        options.num_queries,
    )?;
    println!("Lincheck c verified");
    verify_rowcheck_proof(
        &verifier_key,
        proof.rowcheck_proof,
        &mut public_coin,
        initial_evals,
        options.num_queries,
    )?;
    println!("Rowcheck verified");
    Ok(())
}*/

#[cfg_attr(feature = "flame_it", flame("fractal_verifier"))]
pub fn verify_layered_fractal_proof_from_top<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
    proof: &TopLevelProof<B, E, H>,
    pub_inputs_bytes: &Vec<u8>,
    options: &FractalOptions<B>,
) -> Result<(), FractalVerifierError> {
    let mut accumulator_verifier: AccumulatorVerifier<B, E, H> = AccumulatorVerifier::new(
        options.evaluation_domain.len(),
        options.eta,
        options.evaluation_domain.clone(),
        options.num_queries,
        options.fri_options.clone(),
        pub_inputs_bytes.clone(),
    );

    // draw queries using only the last iop layer commit and the public input.
    // this helps keep the rngs in sync, but proper chaining of layers needs to be checked elsewhere!
    let query_seed = proof.layer_commitments[1];
    let mut coin = RandomCoin::<B, H>::new(&pub_inputs_bytes);
    coin.reseed(query_seed);

    let query_indices = coin
        .draw_integers(options.num_queries, options.evaluation_domain.len())
        .expect("failed to draw query position");

    verify_decommitments(
        &verifier_key,
        &proof,
        &query_indices,
        &mut accumulator_verifier,
    )?;

    let lincheck_proof = other_parse_proofs(&verifier_key, &proof, &pub_inputs_bytes);
    let fractal_proof = parse_proofs_for_subroutines(&verifier_key, &proof, &pub_inputs_bytes);

    verify_layered_fractal_proof(
        &verifier_key,
        fractal_proof,
        lincheck_proof,
        query_indices,
        1,
        &mut accumulator_verifier,
    )?;

    accumulator_verifier.verify_fri_proof(
        proof.layer_commitments[1],
        &proof.low_degree_proof,
        &pub_inputs_bytes,
    )?;

    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("fractal_verifier"))]
pub fn verify_decommitments<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
    proof: &TopLevelProof<B, E, H>,
    query_indices: &Vec<usize>,
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
) -> Result<(), FractalVerifierError> {
    // Verify that the committed preprocessing was queried correctly
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.commitment,
        query_indices,
        &proof.preprocessing_decommitment.0,
        &proof.preprocessing_decommitment.1,
    )?;

    // Verifier that the initial layer was queried correctly
    accumulator_verifier.verify_layer_with_queries(
        proof.initial_commitment,
        query_indices,
        &proof.initial_decommitment.0,
        &proof.initial_decommitment.1,
    )?;

    // Verify that the committed layers were queried correctly
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[0],
        query_indices,
        &proof.layer_decommitments[0].0,
        &proof.layer_decommitments[0].1,
    )?;
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[1],
        query_indices,
        &proof.layer_decommitments[1].0,
        &proof.layer_decommitments[1].1,
    )?;
    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("fractal_verifier"))]
pub fn verify_layered_fractal_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
    proof: LayeredRowcheckProof<B, E>,
    batched_lincheck_proof: BatchedLayeredLincheckProof<B, E>,
    query_indices: Vec<usize>,
    starting_layer: usize,
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
) -> Result<(), FractalVerifierError> {
    verify_layered_rowcheck_proof(
        accumulator_verifier,
        verifier_key,
        &query_indices,
        &proof,
        starting_layer,
    )?;

    verify_layered_lincheck_proof(
        accumulator_verifier,
        verifier_key,
        &query_indices,
        &batched_lincheck_proof,
        starting_layer,
    )?;
    // verify_layered_lincheck_proof(
    //     accumulator_verifier,
    //     verifier_key,
    //     &query_indices,
    //     &proof.lincheck_b,
    //     starting_layer,
    // )?;
    // verify_layered_lincheck_proof(
    //     accumulator_verifier,
    //     verifier_key,
    //     &query_indices,
    //     &proof.lincheck_c,
    //     starting_layer,
    // )?;

    Ok(())
}

/// This function should take as input the full layered fractal proof and return proofs to be passed into the subroutines.
/// Correctness of decommitments should be checked elsewhere.
#[cfg_attr(feature = "flame_it", flame("fractal_verifier"))]
fn parse_proofs_for_subroutines<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
    proof: &TopLevelProof<B, E, H>,
    public_inputs_bytes: &Vec<u8>,
) -> LayeredRowcheckProof<B, E> {
    // Matrix A preprocessing
    let col_a = extract_vec_e(&proof.preprocessing_decommitment.0, 0);
    let row_a = extract_vec_e(&proof.preprocessing_decommitment.0, 1);
    let val_a = extract_vec_e(&proof.preprocessing_decommitment.0, 2);

    // Matrix B preprocessing
    let col_b = extract_vec_e(&proof.preprocessing_decommitment.0, 3);
    let row_b = extract_vec_e(&proof.preprocessing_decommitment.0, 4);
    let val_b = extract_vec_e(&proof.preprocessing_decommitment.0, 5);

    // Matrix C preprocessing
    let col_c = extract_vec_e(&proof.preprocessing_decommitment.0, 6);
    let row_c = extract_vec_e(&proof.preprocessing_decommitment.0, 7);
    let val_c = extract_vec_e(&proof.preprocessing_decommitment.0, 8);

    // get values from the initial layer
    let f_z_vals = extract_vec_e(&proof.initial_decommitment.0, 0);
    let f_az_vals = extract_vec_e(&proof.initial_decommitment.0, 1);
    let f_bz_vals = extract_vec_e(&proof.initial_decommitment.0, 2);
    let f_cz_vals = extract_vec_e(&proof.initial_decommitment.0, 3);

    // get values from the first layer
    let s_vals = extract_vec_e(&proof.layer_decommitments[0].0, 0);
    // let t_alpha_a_vals = extract_vec_e(&proof.layer_decommitments[0].0, 1);
    // let product_sumcheck_a_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[0].0, 2, 3);
    // let t_alpha_b_vals = extract_vec_e(&proof.layer_decommitments[0].0, 4);
    // let product_sumcheck_b_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[0].0, 5, 6);
    // let t_alpha_c_vals = extract_vec_e(&proof.layer_decommitments[0].0, 7);
    // let product_sumcheck_c_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[0].0, 8, 9);

    // // get values from the second layer
    // let matrix_sumcheck_a_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[1].0, 0, 1);
    // let matrix_sumcheck_b_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[1].0, 2, 3);
    // let matrix_sumcheck_c_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[1].0, 4, 5);

    // Sample our own alpha and beta to check the prover
    let mut coin = RandomCoin::<B, H>::new(&public_inputs_bytes);
    coin.reseed(verifier_key.commitment);
    let _: E = coin.draw().expect("failed to draw FRI alpha");
    coin.reseed(proof.initial_commitment);
    let alpha: E = coin.draw().expect("failed to draw FRI alpha");
    coin.reseed(proof.layer_commitments[0]);
    let beta: E = coin.draw().expect("failed to draw FRI alpha");

    //    coin.reseed(proof.layer_commitments[1]);

    // let gammas = &proof.unverified_misc;

    // let lincheck_a_proof = LayeredLincheckProof {
    //     row_vals: row_a,
    //     col_vals: col_a,
    //     val_vals: val_a,
    //     f_z_vals: f_z_vals.clone(),
    //     f_mz_vals: f_az_vals.clone(),
    //     t_alpha_vals: t_alpha_a_vals,
    //     product_sumcheck_vals: product_sumcheck_a_vals,
    //     matrix_sumcheck_vals: matrix_sumcheck_a_vals,
    //     alpha,
    //     beta,
    //     gamma: gammas[0],
    // };

    // let lincheck_b_proof = LayeredLincheckProof {
    //     row_vals: row_b,
    //     col_vals: col_b,
    //     val_vals: val_b,
    //     f_z_vals: f_z_vals.clone(),
    //     f_mz_vals: f_bz_vals.clone(),
    //     t_alpha_vals: t_alpha_b_vals,
    //     product_sumcheck_vals: product_sumcheck_b_vals,
    //     matrix_sumcheck_vals: matrix_sumcheck_b_vals,
    //     alpha,
    //     beta,
    //     gamma: gammas[1],
    // };

    // let lincheck_c_proof = LayeredLincheckProof {
    //     row_vals: row_c,
    //     col_vals: col_c,
    //     val_vals: val_c,
    //     f_z_vals: f_z_vals.clone(),
    //     f_mz_vals: f_cz_vals.clone(),
    //     t_alpha_vals: t_alpha_c_vals,
    //     product_sumcheck_vals: product_sumcheck_c_vals,
    //     matrix_sumcheck_vals: matrix_sumcheck_c_vals,
    //     alpha,
    //     beta,
    //     gamma: gammas[2],
    // };

    LayeredRowcheckProof {
        f_z_vals,
        f_az_vals,
        f_bz_vals,
        f_cz_vals,
        s_vals,
    }
    // LayeredRowcheckProof {

    // }
    // LayeredFractalProof {
    //     rowcheck: rowcheck_proof,
    //     lincheck_a: lincheck_a_proof,
    //     lincheck_b: lincheck_b_proof,
    //     lincheck_c: lincheck_c_proof,
    // }
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
