use crate::errors::SumcheckVerifierError;

use fractal_accumulator::accumulator;
use fractal_proofs::{compute_vanishing_poly, FieldElement, LayeredSumcheckProof, SumcheckProof};

use low_degree_verifier::low_degree_batch_verifier::verify_low_degree_batch_proof;
use low_degree_verifier::low_degree_verifier::verify_low_degree_proof;
use winter_crypto::{ElementHasher, RandomCoin};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

// pub struct SumcheckVerifier<E>  {
//     context: VerifierContext,
//     proof: SumcheckProof,
// }

#[cfg_attr(feature = "flame_it", flame("sumcheck_verifier"))]
pub fn verify_sumcheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    proof: SumcheckProof<B, E, H>,
    g_max_degree: usize,
    e_max_degree: usize,
    public_coin: &mut RandomCoin<B, H>,
    num_queries: usize,
) -> Result<(), SumcheckVerifierError> {
    // let mut public_coin = RandomCoin::new(&[]);
    verify_low_degree_batch_proof(
        proof.batch_proof,
        vec![g_max_degree, e_max_degree],
        public_coin,
        num_queries,
    )?;
    //verify_low_degree_proof(proof.g_proof, g_max_degree, public_coin)?;
    //verify_low_degree_proof(proof.e_proof, e_max_degree, public_coin)?;
    // FIXME: This proof verification should also check that e and g are correct wrt the Az, Bz and Cz.
    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("sumcheck_verifier"))]
pub fn verify_layered_sumcheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    queried_positions: &Vec<usize>,
    proof: LayeredSumcheckProof<B, E>,
    eval_domain_size: usize,
    summing_domain_size: usize,
    eval_domain_offset: B,
    summing_domain_offset: B,
    gamma: E,
    starting_layer: usize,
) -> Result<(), SumcheckVerifierError> {
    let summing_domain_size_u64: u64 = summing_domain_size.try_into().unwrap();
    let summing_domain_size_field = E::from(summing_domain_size_u64);
    let l_field_base = E::from(B::get_root_of_unity(
        eval_domain_size.trailing_zeros().try_into().unwrap(),
    ));
    let eta = summing_domain_offset;
    for i in 0..proof.numerator_vals.len() {
        let position_u64: u64 = queried_positions[i].try_into().unwrap();
        let x_val =
            l_field_base.exp(E::PositiveInteger::from(position_u64)) * E::from(eval_domain_offset);
        let denom_val = compute_vanishing_poly::<E>(x_val, E::from(eta), summing_domain_size);
        let lhs = ((((x_val * proof.sumcheck_g_vals[i]) + (gamma / summing_domain_size_field))
            * proof.denominator_vals[i])
            - proof.numerator_vals[i])
            / denom_val;
        if lhs != proof.sumcheck_e_vals[i] {
            println!("lhs = {:?}, e = {:?}", lhs, proof.sumcheck_e_vals[i]);
            return Err(SumcheckVerifierError::ConsistentValuesErr(i));
        }
    }
    Ok(())
}
