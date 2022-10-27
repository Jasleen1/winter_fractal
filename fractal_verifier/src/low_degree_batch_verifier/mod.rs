use crate::errors::LowDegreeVerifierError;

use fractal_proofs::{FieldElement, LowDegreeBatchProof, polynom};
use fractal_utils::polynomial_utils::*;
use winter_crypto::{ElementHasher, RandomCoin};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

pub fn verify_low_degree_batch_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    proof: LowDegreeBatchProof<B, E, H>, max_degrees: Vec<usize>, public_coin: &mut RandomCoin<B,H>
) -> Result<(), LowDegreeVerifierError> {
    let mut channel = DefaultVerifierChannel::<E, H>::new(
        proof.fri_proof,
        proof.commitments,
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;

    //todo: need to be able to sample these throughout the protocol like for the batch verifier
    //todo: need to sample from the extension field?
    let mut alphas = Vec::new();
    let mut betas = Vec::new();
    for i in 0..2*max_degrees.len(){
        //todo: this doesn't mutate public coin
        alphas.push(public_coin.draw().unwrap());
        betas.push(public_coin.draw().unwrap());
    }
    let fri_verifier = FriVerifier::<B, E, DefaultVerifierChannel<E, H>, H>::new(
        &mut channel,
        public_coin,
        proof.options.clone(),
        proof.fri_max_degree,
    )?;
    //todo, are the queried position ever checked?
    fri_verifier.verify(&mut channel, &proof.composed_queried_evaluations, &proof.queried_positions)?;
    //todo: merkle branches are never verified
    verify_lower_degree_batch::<B, E, H>(proof.options.blowup_factor() * (proof.fri_max_degree+1),
    max_degrees, proof.fri_max_degree, proof.all_unpadded_queried_evaluations,
    proof.composed_queried_evaluations, proof.queried_positions.clone(), alphas, betas)?;
    Ok(())
}

fn verify_lower_degree_batch<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(eval_domain_size: usize, original_degrees: Vec<usize>, fri_max_degree: usize, 
    original_evals: Vec<Vec<E>>, final_evals: Vec<E>, positions: Vec<usize>, alphas: Vec<E>, betas: Vec<E>) -> Result<(), LowDegreeVerifierError> {
    let eval_domain_base = E::from(B::get_root_of_unity(eval_domain_size.trailing_zeros()));
    let eval_domain_pows = positions.iter().map(|&x| x as u64).collect::<Vec<u64>>();
    let eval_domain_elts = eval_domain_pows.iter().map(|&x| eval_domain_base.exp(E::PositiveInteger::from(x))).collect::<Vec<E>>();
    
    //todo: use length of queried positions here
    let mut reconstructed_evals = vec![E::ZERO;eval_domain_elts.len()];
    for pos in 0..original_degrees.len(){
        let comp_poly = get_randomized_complementary_poly::<E>(original_degrees[pos], fri_max_degree, alphas[pos], betas[pos]);
        let eval_domain_evals = polynom::eval_many(&comp_poly, &eval_domain_elts);
        for pos2 in 0..eval_domain_elts.len(){
            reconstructed_evals[pos2] += original_evals[pos][pos2] * eval_domain_evals[pos2];
        }
    }
    for (pos, _) in eval_domain_elts.iter().enumerate() {
        if reconstructed_evals[pos] != final_evals[pos] {
            println!("Position {}", pos);
            println!("reconstructed_evals = {:?}", reconstructed_evals);
            println!("Final evals = {:?}", final_evals[pos]);
            return Err(LowDegreeVerifierError::PaddingErr);
        }
    }
    Ok(())
}

#[cfg(test)]
mod test{
    use crate::{low_degree_batch_prover::LowDegreeBatchProver};
    use super::verify_low_degree_batch_proof;
    use fractal_proofs::{FieldElement, SumcheckProof};
    use fractal_prover::prover_channel::FractalProverChannel;
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{DefaultVerifierChannel, FriVerifier, FriOptions, DefaultProverChannel};
    use winter_math::StarkField;
    use winter_math::fields::f64::BaseElement;
    use winter_crypto::hashers::Rp64_256;
    use winter_math::utils;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn run_test_low_degree_proof(){
        test_low_degree_proof::<BaseElement, BaseElement, Rp64_256>();
    }

    fn test_low_degree_proof<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
        >() {
        let lde_blowup = 4;
        let num_queries = 16;
        let fri_options = FriOptions::new(lde_blowup, 4, 32);
        let max_degree:usize = 63;
        let l_field_size: usize = 4 * max_degree.next_power_of_two();
        let l_field_base = B::get_root_of_unity(l_field_size.trailing_zeros());
        let evaluation_domain = utils::get_power_series(l_field_base, l_field_size);

        //TODO: RandomCoin needs to be able to sample from the extension field, but currently can't
        let mut public_coin = RandomCoin::<B,H>::new(&[]);
        let mut channel = FractalProverChannel::<B,E,H>::new(evaluation_domain.len(), num_queries);
        let mut prover = LowDegreeBatchProver::<B, E, H>::new(&evaluation_domain, fri_options.clone());

        let max_degrees: Vec<usize> = vec![14, 63, 29, 31];
        let mut polys: Vec<Vec<B>> = Vec::new();
        for degree in max_degrees.iter(){
            let poly = nonrand_poly(*degree);
            prover.add_polynomial(&poly, *degree, &mut channel);
            polys.push(poly);
        }

        let proof = prover.generate_proof(&mut channel);
        assert!(verify_low_degree_batch_proof(proof, max_degrees, &mut public_coin).is_ok());
    }

    // a random-ish polynomial that isn't actually random at all. Instead, it uses the system clock since that doesn't require a new crate import
    fn nonrand_poly<B: StarkField>(degree: usize) -> Vec<B>{
        let mut out: Vec<B> = Vec::new();
        for _ in 0..degree+1{
            let nanos = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .subsec_nanos();
            out.push(B::from(nanos as u128));
        }
        out
    }

}