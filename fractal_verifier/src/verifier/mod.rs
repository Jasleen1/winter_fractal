use crate::{
    accumulator_verifier::{self, AccumulatorVerifier},
    errors::FractalVerifierError,
};

use fractal_indexer::snark_keys::*;
use fractal_proofs::{
    FieldElement, FractalProof, LayeredFractalProof, MultiEval, MultiPoly, StarkField,
};

use fractal_prover::{channel::DefaultFractalProverChannel, FractalOptions};
use log::debug;
use winter_crypto::{ElementHasher, RandomCoin};

use crate::{lincheck_verifier::verify_lincheck_proof, rowcheck_verifier::verify_rowcheck_proof};

pub fn verify_fractal_proof<
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
}

fn verify_layered_fractal_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: VerifierKey<B, E, H>,
    proof: LayeredFractalProof<B, E, H>,
    pub_inputs_bytes: H::Digest,
    options: FractalOptions<B>,
) -> Result<(), FractalVerifierError> {
    let mut accumulator_verifier: AccumulatorVerifier<B, E, H> = AccumulatorVerifier::new(
        options.evaluation_domain.len(),
        options.eta,
        options.evaluation_domain.clone(),
        options.num_queries,
        options.fri_options.clone(),
    );

    // Here, we check that the sent over proofs verify with respect to the sent commitments.
    // Step A: draw queries
    let query_seed = proof.layer_commitments[2];
    let mut coin = RandomCoin::<B, H>::new(&vec![]);
    coin.reseed(query_seed);
    let query_indices = coin
        .draw_integers(options.num_queries, options.evaluation_domain.len())
        .expect("failed to draw query position");
    // Step B: Verify that the preprocessing was queried correctly
    // Do everything for matrix A preprocessing
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_a_commitments.row_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_a[0].0,
        &proof.preprocessing_decommits_a[0].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_a_commitments.col_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_a[1].0,
        &proof.preprocessing_decommits_a[1].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_a_commitments.val_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_a[2].0,
        &proof.preprocessing_decommits_a[2].1,
    );
    // Do everything for matrix B preprocessing
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_b_commitments.row_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_b[0].0,
        &proof.preprocessing_decommits_b[0].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_b_commitments.col_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_b[1].0,
        &proof.preprocessing_decommits_b[1].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_b_commitments.val_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_b[2].0,
        &proof.preprocessing_decommits_b[2].1,
    );
    // Do everything for matrix C preprocessing
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_c_commitments.row_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_c[0].0,
        &proof.preprocessing_decommits_c[0].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_c_commitments.col_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_c[1].0,
        &proof.preprocessing_decommits_c[1].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.matrix_c_commitments.val_poly_commitment,
        &query_indices,
        &proof.preprocessing_decommits_c[2].0,
        &proof.preprocessing_decommits_c[2].1,
    );

    // Step C: Verify that the committed layers were queried correctly
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[0],
        &query_indices,
        &proof.layer_decommits[0].0,
        &proof.layer_decommits[0].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[1],
        &query_indices,
        &proof.layer_decommits[1].0,
        &proof.layer_decommits[1].1,
    );
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[2],
        &query_indices,
        &proof.layer_decommits[2].0,
        &proof.layer_decommits[2].1,
    );

    // Parse all the various proof fields into inputs for the 3 layered_lincheck_verifier instances and for the rowcheck instance.

    // Verify FRI proof.

    Ok(())
}
