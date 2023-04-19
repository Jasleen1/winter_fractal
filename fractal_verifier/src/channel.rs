use fractal_proofs::{DeserializationError, FieldElement};
use winter_crypto::{BatchMerkleProof, ElementHasher};
use winter_fri::{FriProof, VerifierChannel};

/// Provides a default implementation of the [VerifierChannel] trait.
///
/// Default verifier channel can be instantiated directly from a [FriProof] struct.
///
/// Though this implementation is primarily intended for testing purposes, it can be used in
/// production use cases as well.
pub struct DefaultFractalVerifierChannel<
    E: FieldElement,
    H: ElementHasher<BaseField = E::BaseField>,
> {
    layer_commitments: Vec<H::Digest>,
    layer_proofs: Vec<BatchMerkleProof<H>>,
    layer_queries: Vec<Vec<E>>,
    remainder: Vec<E>,
    num_partitions: usize,
}

impl<E, H> DefaultFractalVerifierChannel<E, H>
where
    E: FieldElement,
    H: ElementHasher<BaseField = E::BaseField>,
{
    /// Builds a new verifier channel from the specified [FriProof].
    ///
    /// # Errors
    /// Returns an error if the specified `proof` could not be parsed correctly.
    pub fn new(
        proof: FriProof,
        layer_commitments: Vec<H::Digest>,
        domain_size: usize,
        folding_factor: usize,
    ) -> Result<Self, DeserializationError> {
        let num_partitions = proof.num_partitions();

        let remainder = proof.parse_remainder()?;
        let (layer_queries, layer_proofs) =
            proof.parse_layers::<H, E>(domain_size, folding_factor)?;

        Ok(DefaultFractalVerifierChannel {
            layer_commitments,
            layer_proofs,
            layer_queries,
            remainder,
            num_partitions,
        })
    }
}

impl<E, H> VerifierChannel<E> for DefaultFractalVerifierChannel<E, H>
where
    E: FieldElement,
    H: ElementHasher<BaseField = E::BaseField>,
{
    type Hasher = H;

    fn read_fri_num_partitions(&self) -> usize {
        self.num_partitions
    }

    fn read_fri_layer_commitments(&mut self) -> Vec<H::Digest> {
        self.layer_commitments.drain(..).collect()
    }

    fn take_next_fri_layer_proof(&mut self) -> BatchMerkleProof<H> {
        self.layer_proofs.remove(0)
    }

    fn take_next_fri_layer_queries(&mut self) -> Vec<E> {
        self.layer_queries.remove(0)
    }

    fn take_fri_remainder(&mut self) -> Vec<E> {
        self.remainder.clone()
    }
}
