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
    println!("sumcheck verifier: proof.num_evaluations:{} ", proof.num_evaluations);

    let mut g_channel = DefaultVerifierChannel::<E, H>::new(
        proof.g_proof,
        proof.g_queried.queried_proofs[0].clone(),
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;
    println!("proof.num_evaluations={}", proof.num_evaluations);

    let g_verifier = FriVerifier::<B, E, DefaultVerifierChannel<E, H>, H>::new(
        &mut g_channel,
        &mut public_coin,
        proof.options.clone(),
        63//proof.g_max_degree-1, //63 (was 31) but should be 63
        // verifier_key.params.max_degree - 1,
    )?;
    println!("lincheck max_poly_degree {}", proof.g_max_degree-1);
    let g_queried_evals = proof.g_queried.queried_evals;
    //todo, are the queried position ever checked?
    println!("Sumcheck verifier indexes: {:?}", &proof.queried_positions);
    println!("Sumcheck verifier g_queried_evals: {:?}", &g_queried_evals);
    //println!("g_channel.layer_proofs.domain_size={}", g_channel.layer_proofs.dom);
    println!("g_verifier.domain_size={}", g_verifier.domain_size());
    g_verifier.verify(&mut g_channel, &g_queried_evals, &proof.queried_positions)?;
    println!("verified g");

    let mut e_channel = DefaultVerifierChannel::<E, H>::new(
        proof.e_proof,
        proof.e_queried.queried_proofs[0].clone(),
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;
    println!("proof.e_max_degree: {} ", &proof.e_max_degree);
    let e_verifier = FriVerifier::<B, E, DefaultVerifierChannel<E, H>, H>::new(
        &mut e_channel,
        &mut public_coin,
        proof.options.clone(),
        63//proof.e_max_degree-1
    )?;

    let e_queried_evals = proof.e_queried.queried_evals;
    println!("calling verify");
    Ok(e_verifier.verify(&mut e_channel, &e_queried_evals, &proof.queried_positions)?)
}
