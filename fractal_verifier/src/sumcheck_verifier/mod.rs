use crate::errors::SumcheckVerifierError;

use fractal_proofs::{FieldElement, SumcheckProof};

use crate::low_degree_batch_verifier::verify_low_degree_batch_proof;
use crate::low_degree_verifier::verify_low_degree_proof;
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
