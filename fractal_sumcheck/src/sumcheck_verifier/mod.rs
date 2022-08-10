use crate::errors::SumcheckVerifierError;

use fractal_proofs::{FieldElement, SumcheckProof};

use winter_crypto::{ElementHasher, RandomCoin};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

// pub struct SumcheckVerifier<E>  {
//     context: VerifierContext,
//     proof: SumcheckProof,
// }

pub fn verify_sumcheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    proof: SumcheckProof<B, E, H>,
) -> Result<(), SumcheckVerifierError> {
    // let mut public_coin_seed = Vec::new();
    // proof.write_into(&mut public_coin_seed);
    // let mut public_coin = RandomCoin::new(&public_coin_seed);
    let mut public_coin = RandomCoin::new(&[]);
    println!("proof.g_max_degree = {}", proof.g_max_degree);

    let mut g_channel = DefaultVerifierChannel::<E, H>::new(
        proof.g_proof,
        proof.g_queried.queried_proofs[0].clone(),
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;

    let g_verifier = FriVerifier::<B, E, DefaultVerifierChannel<E, H>, H>::new(
        &mut g_channel,
        &mut public_coin,
        proof.options.clone(),
        proof.g_max_degree - 1,
    )?;
    let g_queried_evals = proof.g_queried.queried_evals;
    g_verifier.verify(&mut g_channel, &g_queried_evals, &proof.queried_positions)?;

    let mut e_channel = DefaultVerifierChannel::<E, H>::new(
        proof.e_proof,
        proof.e_queried.queried_proofs[0].clone(),
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;

    let e_verifier = FriVerifier::<B, E, DefaultVerifierChannel<E, H>, H>::new(
        &mut e_channel,
        &mut public_coin,
        proof.options.clone(),
        proof.e_max_degree - 1,
    )?;

    let e_queried_evals = proof.e_queried.queried_evals;
    Ok(e_verifier.verify(&mut e_channel, &e_queried_evals, &proof.queried_positions)?)
}
