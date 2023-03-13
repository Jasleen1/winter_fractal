use fractal_proofs::{MultiPoly, LowDegreeBatchProof};
use winter_math::{fft, FieldElement, StarkField};
use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
use std::{convert::TryInto, marker::PhantomData};
use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree};
use fractal_utils::polynomial_utils::MultiEval;

use crate::{low_degree_batch_prover::LowDegreeBatchProver, channel::DefaultFractalProverChannel};

pub struct Accumulator<
B: StarkField,
E: FieldElement<BaseField = B>,
H: ElementHasher + ElementHasher<BaseField = B>,
> {
pub evaluation_domain_len: usize,
pub offset: B,
pub evaluation_domain: Vec<B>,
pub num_queries: usize,
pub fri_options: FriOptions,
pub coefficients: Vec<Vec<B>>,
pub coefficients_e: Vec<Vec<E>>,
pub max_degrees: Vec<usize>,
pub max_degrees_e: Vec<usize>,
_h: PhantomData<H>,
}

impl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> Accumulator<B, E, H>{
    pub fn new(evaluation_domain_len: usize, offset: B, evaluation_domain: Vec<B>, num_queries: usize, fri_options:FriOptions) -> Self{
        Self {
            evaluation_domain_len,
            offset,
            evaluation_domain,
            num_queries,
            fri_options,
            coefficients: Vec::new(),
            coefficients_e: Vec::new(),
            max_degrees: Vec::new(),
            max_degrees_e: Vec::new(),
            _h: PhantomData,
        }
    }
    pub fn add_polynomial(&mut self, coefficients: Vec<B>, max_degree: usize) -> (){
        self.coefficients.push(coefficients);
        self.max_degrees.push(max_degree);
    }
    pub fn add_polynomial_e(&mut self, coefficients: Vec<E>, max_degree: usize) -> (){
        self.coefficients_e.push(coefficients);
        self.max_degrees_e.push(max_degree);
    }
    pub fn commit_layer(&self) -> <H>::Digest{
        let mut multi_eval = MultiEval::<B,E,H>::new(self.coefficients.clone(), self.coefficients_e.clone(), self.evaluation_domain_len, self.offset);
        multi_eval.commit_polynomial_evaluations().unwrap();
        multi_eval.get_commitment().unwrap().clone()
    }
    pub fn draw_queries(&self, count: usize) -> Vec<E> {
        let channel_state = self.commit_layer();
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            Vec::new()
        );
        channel.commit_fractal_iop_layer(channel_state);
        let queries = (0..count).map(|_| channel.draw_fri_alpha()).collect();
        queries
    }
    pub fn create_fri_proof(&self) -> LowDegreeBatchProof<B,E,H>{
        let channel_state = self.commit_layer();
        let mut channel = &mut DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            Vec::new()
        );
        channel.public_coin.reseed(channel_state);
        let mut low_degree_prover = LowDegreeBatchProver::<B,E,H>::new(&self.evaluation_domain, self.fri_options.clone());
        for i in 0..self.max_degrees.len() {
            low_degree_prover.add_polynomial(self.coefficients.get(i).unwrap(), *self.max_degrees.get(i).unwrap(), &mut channel);
        }
        for i in 0..self.max_degrees_e.len() {
            low_degree_prover.add_polynomial_e(self.coefficients_e.get(i).unwrap(), *self.max_degrees_e.get(i).unwrap(), &mut channel);
        }

        low_degree_prover.generate_proof(&mut channel)
    }
}

#[cfg(test)]
mod test {
    use fractal_proofs::{MultiPoly, BaseElement, utils, fields::QuadExtension};
    use winter_math::{fft, FieldElement, StarkField};
    use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
    use std::{convert::TryInto, marker::PhantomData};
    use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree, hashers::Blake3_256};
    use fractal_utils::polynomial_utils::MultiEval;

    use super::Accumulator;
    #[test]
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
    }
}