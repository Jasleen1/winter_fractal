use fractal_proofs::{LowDegreeBatchProof, MultiPoly};
use fractal_utils::polynomial_utils::MultiEval;
use std::{convert::TryInto, marker::PhantomData};
use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree, RandomCoin};
use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
use winter_math::{fft, FieldElement, StarkField};

use crate::low_degree_batch_verifier::verify_low_degree_batch_proof;

pub struct AccumulatorVerifier<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub evaluation_domain_len: usize,
    pub offset: B,
    pub evaluation_domain: Vec<B>,
    pub num_queries: usize,
    pub fri_options: FriOptions,
    pub max_degrees: Vec<usize>,
    pub public_coin: RandomCoin<B, H>,
    _e: PhantomData<E>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > AccumulatorVerifier<B, E, H>
{
    // should take pub_bytes here?
    pub fn new(
        evaluation_domain_len: usize,
        offset: B,
        evaluation_domain: Vec<B>,
        num_queries: usize,
        fri_options: FriOptions,
    ) -> Self {
        Self {
            evaluation_domain_len,
            offset,
            evaluation_domain,
            num_queries,
            fri_options,
            max_degrees: Vec::new(),
            public_coin: RandomCoin::<B, H>::new(&vec![]),
            _e: PhantomData,
        }
    }

    pub fn add_constraint(&mut self, max_degree: usize) {
        self.max_degrees.push(max_degree);
    }

    // verify batch incluion proof, update channel state
    pub fn verify_layer(
        &mut self,
        layer_commit: H::Digest,
        decommit: Vec<Vec<E>>,
        proof: BatchMerkleProof<H>,
    ) -> bool {
        let mut coin = RandomCoin::<B, H>::new(&vec![]);
        coin.reseed(layer_commit);
        let indices = coin
            .draw_integers(self.num_queries, self.evaluation_domain_len)
            .expect("failed to draw query position");
        MultiEval::<B, E, H>::batch_verify_values_and_proofs_at(
            decommit, // todo: this should be decommit once this function is fixed,
            &proof.get_root(&indices).unwrap(), //todo: is this okay
            &proof,
            indices.to_vec(),
        )
        .is_ok()
    }

    // run at the end
    pub fn verify_fri_proof(
        &mut self,
        last_layer_commit: H::Digest,
        proof: LowDegreeBatchProof<B, E, H>,
    ) -> bool {
        let mut coin = RandomCoin::<B, H>::new(&vec![]);
        coin.reseed(last_layer_commit);
        verify_low_degree_batch_proof(proof, self.max_degrees.clone(), &mut coin, self.num_queries)
            .is_ok()
    }
}

#[cfg(test)]
mod test {
    use fractal_proofs::{fields::QuadExtension, utils, BaseElement, MultiPoly};
    use fractal_utils::polynomial_utils::MultiEval;
    use std::{convert::TryInto, marker::PhantomData};
    use winter_crypto::{hashers::Blake3_256, BatchMerkleProof, ElementHasher, MerkleTree};
    use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
    use winter_math::{fft, FieldElement, StarkField};

    use super::AccumulatorVerifier;
    /*#[test]
    fn test_accumulator() {
        let lde_blowup = 4;
        let num_queries = 16;
        let fri_options = FriOptions::new(lde_blowup, 4, 32);
        let max_degree: usize = 63;
        let l_field_size: usize = 4 * max_degree.next_power_of_two();
        let l_field_base = BaseElement::get_root_of_unity(l_field_size.trailing_zeros());
        let offset = BaseElement::ONE;
        let evaluation_domain = utils::get_power_series(l_field_base, l_field_size);
        let acc = Accumulator::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>::new(evaluation_domain.len(), offset, evaluation_domain, num_queries, fri_options);
        let alphas = acc.draw_queries(20);
        assert!(alphas.len() == 20)
    }*/
}
