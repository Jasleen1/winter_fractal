use crate::errors::FractalVerifierError;

use fractal_indexer::snark_keys::*;
use fractal_proofs::{FieldElement, FractalProof, StarkField};

use fractal_prover::channel::DefaultFractalProverChannel;
use log::debug;
use winter_crypto::{ElementHasher, RandomCoin};

use crate::{lincheck_verifier::verify_lincheck_proof, rowcheck_verifier::verify_rowcheck_proof};

pub fn verify_fractal_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: VerifierKey<H, B>,
    proof: FractalProof<B, E, H>,
    pub_inputs_bytes: Vec<u8>,
) -> Result<(), FractalVerifierError> {
    let mut public_coin = RandomCoin::<_, H>::new(&pub_inputs_bytes);
    let expected_alpha: B = public_coin.draw().expect("failed to draw OOD point");
    // let mut channel = DefaultFractalProverChannel::new();

    debug!(
        "Lincheck a indexes: {:?}",
        &proof.lincheck_a.products_sumcheck_proof.queried_positions
    );
    verify_lincheck_proof(
        &verifier_key,
        proof.lincheck_a,
        expected_alpha,
        &mut public_coin,
    )?;
    println!("Lincheck a verified");
    verify_lincheck_proof(
        &verifier_key,
        proof.lincheck_b,
        expected_alpha,
        &mut public_coin,
    )?;
    println!("Lincheck b verified");
    verify_lincheck_proof(
        &verifier_key,
        proof.lincheck_c,
        expected_alpha,
        &mut public_coin,
    )?;
    println!("Lincheck c verified");
    verify_rowcheck_proof(&verifier_key, proof.rowcheck_proof, &mut public_coin)?;
    println!("Rowcheck verified");
    Ok(())
}
