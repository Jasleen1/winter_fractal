use crate::{
    channel::DefaultFractalVerifierChannel, errors::RowcheckVerifierError,
    low_degree_verifier::verify_low_degree_proof,
};

use fractal_indexer::snark_keys::VerifierKey;
use fractal_proofs::{get_complementary_poly, polynom, FieldElement, RowcheckProof, TryInto};

use log::debug;
use winter_crypto::{ElementHasher, MerkleTree, RandomCoin};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

pub fn verify_rowcheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<H, B>,
    proof: RowcheckProof<B, E, H>,
    public_coin: &mut RandomCoin<B, H>,
) -> Result<(), RowcheckVerifierError> {
    verify_low_degree_proof(
        proof.s_proof,
        verifier_key.params.max_degree - 1,
        public_coin,
    )?;

    Ok(())
}
