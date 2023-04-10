use crate::errors::FractalVerifierError;

use fractal_indexer::snark_keys::*;
use fractal_proofs::{FieldElement, FractalProof, MultiEval, MultiPoly, StarkField};

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
