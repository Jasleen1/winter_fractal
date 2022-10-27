use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree};
use winter_fri::VerifierChannel as VerifierChannel;
use winter_math::{FieldElement, StarkField};

pub struct FractalVerifierChannel<E: FieldElement, H: ElementHasher<BaseField = E::BaseField>> {
    
}

impl<E, H> VerifierChannel<E> for FractalVerifierChannel<E, H> 
where
    E: FieldElement,
    H: ElementHasher<BaseField = E::BaseField>,
{
    type Hasher = H;
}