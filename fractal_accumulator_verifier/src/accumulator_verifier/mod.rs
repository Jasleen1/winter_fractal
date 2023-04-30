use crate::errors::AccumulatorVerifierError;
use fractal_proofs::{LowDegreeBatchProof, MultiPoly};
use fractal_utils::polynomial_utils::MultiEval;
use low_degree_verifier::low_degree_batch_verifier::verify_low_degree_batch_proof;
use std::{convert::TryInto, marker::PhantomData};
use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree, RandomCoin};
use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel, VerifierError};
use winter_math::{fft, FieldElement, StarkField}; //, FractalVerifierError};

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
    pub max_degrees_by_layer: Vec<Vec<usize>>,
    //pub public_coin: RandomCoin<B, H>,
    pub public_inputs_bytes: Vec<u8>,
    _e: PhantomData<E>,
    _h: PhantomData<H>,
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
        public_inputs_bytes: Vec<u8>,
    ) -> Self {
        Self {
            evaluation_domain_len,
            offset,
            evaluation_domain,
            num_queries,
            fri_options,
            max_degrees: Vec::new(),
            max_degrees_by_layer: Vec::new(),
            public_inputs_bytes,
            //public_coin: RandomCoin::<B, H>::new(&pub_inputs_bytes), //todo: this is unused
            _e: PhantomData,
            _h: PhantomData,
        }
    }

    #[cfg_attr(feature = "flame_it", flame("accumulator_verifier"))]
    pub fn add_constraint(&mut self, max_degree: usize, current_layer: usize) {
        self.max_degrees.push(max_degree);
        while self.max_degrees_by_layer.len() <= current_layer {
            self.max_degrees_by_layer.push(Vec::new());
        }
        self.max_degrees_by_layer[current_layer].push(max_degree);
    }

    // verify batch incluion proof, update channel state
    #[cfg_attr(feature = "flame_it", flame("accumulator_verifier"))]
    pub fn verify_layer(
        &mut self,
        layer_commit: H::Digest,
        query_seed: H::Digest,
        decommit: &Vec<Vec<E>>,
        proof: &BatchMerkleProof<H>,
    ) -> Result<(), AccumulatorVerifierError> {
        let mut coin = RandomCoin::<B, H>::new(&self.public_inputs_bytes);
        coin.reseed(query_seed);
        let indices = coin
            .draw_integers(self.num_queries, self.evaluation_domain_len)
            .expect("failed to draw query position");
        let claimed_root = proof.get_root(&indices).unwrap();
        if layer_commit != claimed_root {
            return Err(AccumulatorVerifierError::CommitMatchErr(format!(
                "Claimed root = {:?}, Layer commitment = {:?}",
                claimed_root, claimed_root
            )));
        }
        MultiEval::<B, E, H>::batch_verify_values_and_proofs_at(
            decommit,      // todo: this should be decommit once this function is fixed,
            &claimed_root, //todo: is this okay
            &proof,
            &indices.to_vec(),
        )?;
        Ok(())
    }

    // verify batch incluion proof, update channel state
    #[cfg_attr(feature = "flame_it", flame("accumulator_verifier"))]
    pub fn verify_layer_with_queries(
        &mut self,
        layer_commit: H::Digest,
        query_indices: &Vec<usize>,
        decommit: &Vec<Vec<E>>,
        proof: &BatchMerkleProof<H>,
    ) -> Result<(), AccumulatorVerifierError> {
        let claimed_root = proof.get_root(&query_indices).unwrap();
        if layer_commit != claimed_root {
            return Err(AccumulatorVerifierError::CommitMatchErr(format!(
                "Claimed root = {:?}, Layer commitment = {:?}",
                claimed_root, claimed_root
            )));
        }

        MultiEval::<B, E, H>::batch_verify_values_and_proofs_at(
            &decommit,     // todo: this should be decommit once this function is fixed,
            &claimed_root, //todo: is this okay
            &proof,
            &query_indices.to_vec(),
        )?;

        Ok(())
    }

    // verify batch incluion proof, update channel state
    #[cfg_attr(feature = "flame_it", flame("accumulator_verifier"))]
    pub fn verify_transposed_layer_with_queries(
        &mut self,
        layer_commit: H::Digest,
        query_indices: Vec<usize>,
        decommit: Vec<[E; 1]>,
        proof: &BatchMerkleProof<H>,
    ) -> Result<(), AccumulatorVerifierError> {
        let claimed_root = proof.get_root(&query_indices).unwrap();
        if layer_commit != claimed_root {
            return Err(AccumulatorVerifierError::CommitMatchErr(format!(
                "Claimed root = {:?}, Layer commitment = {:?}",
                claimed_root, claimed_root
            )));
        }

        MultiEval::<B, E, H>::batch_verify_transposed_values_and_proofs_at(
            decommit,      // todo: this should be decommit once this function is fixed,
            &claimed_root, //todo: is this okay
            &proof,
            &query_indices.to_vec(),
        )?;

        Ok(())
    }

    // run at the end
    #[cfg_attr(feature = "flame_it", flame("accumulator_verifier"))]
    pub fn verify_fri_proof(
        &mut self,
        last_layer_commit: H::Digest,
        proof: &LowDegreeBatchProof<B, E, H>,
        pub_inputs_bytes: &Vec<u8>,
    ) -> Result<(), AccumulatorVerifierError> {
        let mut coin = RandomCoin::<B, H>::new(&pub_inputs_bytes);
        coin.reseed(last_layer_commit);
        let mut max_degrees = Vec::new();
        for v in self.max_degrees_by_layer.iter() {
            max_degrees.extend(v);
        }
        let res = verify_low_degree_batch_proof(proof, max_degrees, &mut coin, self.num_queries);
        println!("res = {:?}", res);
        Ok(res?)
    }

    #[cfg_attr(feature = "flame_it", flame("accumulator_verifier"))]
    pub fn get_query_indices(
        &self,
        query_seed: H::Digest,
        pub_inputs_bytes: Vec<u8>,
    ) -> Result<Vec<usize>, AccumulatorVerifierError> {
        let mut coin = RandomCoin::<B, H>::new(&pub_inputs_bytes);
        coin.reseed(query_seed);
        let indices = coin.draw_integers(self.num_queries, self.evaluation_domain_len)?;
        Ok(indices)
    }
}

/*
pub struct FriAccumulatorVerifier<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub evaluation_domain_len: usize,
    pub offset: B,
    pub evaluation_domain: Vec<B>,
    pub num_queries: usize,
    pub fri_options: FriOptions,
    pub fri_max_degrees: Vec<usize>,
    pub public_coin: RandomCoin<B, H>,
    _e: PhantomData<E>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > FriAccumulatorVerifier<B, E, H>
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
            fri_max_degrees: Vec::new(),
            public_coin: RandomCoin::<B, H>::new(&vec![]),
            _e: PhantomData,
        }
    }


    // run at the end
    pub fn verify_fri_proof(
        &mut self,
        last_layer_commit: H::Digest,
        proof: LowDegreeBatchProof<B, E, H>,
    ) -> Result<bool, FractalVerifierError> {
        let mut coin = RandomCoin::<B, H>::new(&vec![]);
        coin.reseed(last_layer_commit);
        Ok(verify_low_degree_batch_proof(proof, self.fri_max_degrees.clone(), &mut coin, self.num_queries).is_ok())
    }
}
*/

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
