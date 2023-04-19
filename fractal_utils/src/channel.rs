use std::marker::PhantomData;

use winter_math::{FieldElement, StarkField};
use winter_crypto::{Hasher, RandomCoin};
use winter_fri::ProverChannel;

use winter_crypto::{BatchMerkleProof, ElementHasher};
use winter_fri::{VerifierChannel, FriProof};
use winter_utils::DeserializationError;

/// This file basically contains a replica of [winter_fri::DefaultProverChannel] with some extra functions for our purposes.
/// Provides a default implementation of the [ProverChannel] trait.
///
/// Though this implementation is intended primarily for testing purposes, it can be used in
/// production use cases as well.
pub struct DefaultFractalProverChannel<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub public_coin: RandomCoin<B, H>,
    pub(crate) commitments: Vec<H::Digest>,
    domain_size: usize,
    num_queries: usize,
    _field_element: PhantomData<E>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher>
    DefaultFractalProverChannel<B, E, H>
{
    pub fn set_coin(&mut self, coin: RandomCoin<B, H>) {
        self.public_coin = coin;
    }
    /// Returns a new prover channel instantiated from the specified parameters.
    ///
    /// # Panics
    /// Panics if:
    /// * `domain_size` is smaller than 8 or is not a power of two.
    /// * `num_queries` is zero.
    pub fn new(domain_size: usize, num_queries: usize, pub_inputs_bytes: Vec<u8>) -> Self {
        assert!(
            domain_size >= 8,
            "domain size must be at least 8, but was {}",
            domain_size
        );
        assert!(
            domain_size.is_power_of_two(),
            "domain size must be a power of two, but was {}",
            domain_size
        );
        assert!(
            num_queries > 0,
            "number of queries must be greater than zero"
        );
        let coin_seed = pub_inputs_bytes;
        DefaultFractalProverChannel {
            public_coin: RandomCoin::new(&coin_seed),
            commitments: Vec::new(),
            domain_size,
            num_queries,
            _field_element: PhantomData,
        }
    }

    /// Draws a set of positions at which the polynomial evaluations committed at the first FRI
    /// layer should be queried.
    ///
    /// The positions are pseudo-randomly generated based on the values the prover has written
    /// into this channel.
    ///
    /// # Panics
    /// Panics if the specified number of unique positions could not be drawn from the specified
    /// domain. Both number of queried positions and domain size are specified during
    /// construction of the channel.
    pub fn draw_query_positions(&mut self) -> Vec<usize> {
        self.public_coin
            .draw_integers(self.num_queries, self.domain_size)
            .expect("failed to draw query position")
    }

    /// Returns a list of FRI layer commitments written by the prover into this channel.
    pub fn layer_commitments(&self) -> &[H::Digest] {
        &self.commitments
    }

    pub fn commit_fractal_iop_layer(&mut self, layer_root: H::Digest) {
        self.commitments.push(layer_root);
        self.public_coin.reseed(layer_root);
    }

    pub fn draw_random_b_pt(&mut self) -> B {
        self.public_coin.draw().expect("failed to draw FRI alpha")
    }
}

impl<B, E, H> ProverChannel<E> for DefaultFractalProverChannel<B, E, H>
where
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: Hasher,
{
    type Hasher = H;

    fn commit_fri_layer(&mut self, layer_root: H::Digest) {
        self.commitments.push(layer_root);
        self.public_coin.reseed(layer_root);
    }

    fn draw_fri_alpha(&mut self) -> E {
        self.public_coin.draw().expect("failed to draw FRI alpha")
    }
}

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
