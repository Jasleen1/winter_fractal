use crate::errors::AccumulatorProverError;
use fractal_proofs::{LowDegreeBatchProof, MultiPoly};
use fractal_utils::channel::DefaultFractalProverChannel;
use fractal_utils::polynomial_utils::MultiEval;
use low_degree_prover::low_degree_batch_prover::LowDegreeBatchProver;
use std::{convert::TryInto, marker::PhantomData};
use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree};
use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
use winter_math::{fft, FieldElement, StarkField};

pub struct Accumulator<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub evaluation_domain_len: usize,
    pub eval_domain_offset: B,
    pub evaluation_domain: Vec<B>,
    pub num_queries: usize,
    pub fri_options: FriOptions,
    pub coefficients: Vec<Vec<B>>,
    pub coefficients_ext: Vec<Vec<E>>,
    //pub max_degrees: Vec<usize>,
    pub max_degrees_ext: Vec<usize>,
    pub unchecked_coefficients: Vec<Vec<B>>,
    pub unchecked_coefficients_ext: Vec<Vec<E>>,
    //pub fri_coefficients: Vec<Vec<B>>,
    pub fri_coefficients_ext: Vec<Vec<E>>,
    //pub fri_max_degrees: Vec<usize>,
    pub fri_max_degrees_ext: Vec<usize>,
    pub layer_evals: Vec<MultiEval<B, E, H>>,
    pub public_inputs_bytes: Vec<u8>,
    pub max_degree: usize,
    _h: PhantomData<H>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > Accumulator<B, E, H>
{
    pub fn new(
        evaluation_domain_len: usize,
        eval_domain_offset: B,
        evaluation_domain: Vec<B>,
        num_queries: usize,
        fri_options: FriOptions,
        public_inputs_bytes: Vec<u8>,
        max_degree: usize,
    ) -> Self {
        Self {
            evaluation_domain_len,
            eval_domain_offset,
            evaluation_domain,
            num_queries,
            fri_options,
            coefficients: Vec::new(),
            coefficients_ext: Vec::new(),
            //max_degrees: Vec::new(),
            max_degrees_ext: Vec::new(),
            unchecked_coefficients: Vec::new(),
            unchecked_coefficients_ext: Vec::new(),
            //fri_coefficients: Vec::new(),
            fri_coefficients_ext: Vec::new(),
            //fri_max_degrees: Vec::new(),
            fri_max_degrees_ext: Vec::new(),
            layer_evals: Vec::new(),
            public_inputs_bytes,
            max_degree,
            _h: PhantomData,
        }
    }

    pub fn add_polynomial(&mut self, coefficients: Vec<B>, max_degree: usize) {
        let coeffs_ext: Vec<E> = coefficients.iter().map(|y| E::from(*y)).collect();
        self.coefficients_ext.push(coeffs_ext);
        self.max_degrees_ext.push(max_degree);
    }

    pub fn add_polynomial_e(&mut self, coefficients: Vec<E>, max_degree: usize) {
        self.coefficients_ext.push(coefficients);
        self.max_degrees_ext.push(max_degree);
    }

    // commit to a polynomial which does not need to be part of a degree proof
    pub fn add_unchecked_polynomial(&mut self, coefficients: Vec<B>) {
        self.unchecked_coefficients.push(coefficients);
    }

    pub fn commit_layer(&mut self) -> Result<<H>::Digest, AccumulatorProverError> {
        let mut coeffs_b = self.unchecked_coefficients.clone();
        let mut coeffs_b2 = self.coefficients.clone();
        coeffs_b.append(&mut coeffs_b2);
        let mut multi_eval = MultiEval::<B, E, H>::new(
            coeffs_b,
            self.coefficients_ext.clone(),
            self.evaluation_domain_len,
            self.eval_domain_offset,
        );
        //let mut multi_eval = MultiEval::<B,E,H>::new(self.coefficients.clone(), self.coefficients_ext.clone(), self.evaluation_domain_len, self.offset);
        //self.fri_coefficients.append(&mut self.coefficients.clone());
        //self.fri_max_degrees.append(&mut self.max_degrees.clone());
        self.fri_coefficients_ext.append(&mut self.coefficients_ext);
        self.fri_max_degrees_ext.append(&mut self.max_degrees_ext);
        self.coefficients = Vec::new();
        self.coefficients_ext = Vec::new();
        self.unchecked_coefficients = Vec::new();
        //self.max_degrees = Vec::new();
        self.max_degrees_ext = Vec::new();
        multi_eval.commit_polynomial_evaluations()?;
        let com = *multi_eval.get_commitment()?;
        self.layer_evals.push(multi_eval);
        Ok(com)
    }

    pub fn draw_query_positions(&mut self) -> Result<Vec<usize>, AccumulatorProverError> {
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            self.public_inputs_bytes.clone(), // make sure there's actually chainging between layers
        );
        let latest_eval = self
            .layer_evals
            .last()
            .ok_or(AccumulatorProverError::QueryErr(
                "You tried to query the accumulator before anything was committed".to_string(),
            ))?;
        let coin_val = latest_eval.get_commitment()?;
        channel.commit_fractal_iop_layer(*coin_val);
        let queries = channel.draw_query_positions();
        Ok(queries)
    }

    pub fn draw_queries(&mut self, count: Option<usize>) -> Result<Vec<E>, AccumulatorProverError> {
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            self.public_inputs_bytes.clone(), // make sure there's actually chainging between layers
        );
        let latest_eval = self
            .layer_evals
            .last()
            .ok_or(AccumulatorProverError::QueryErr(
                "You tried to query the accumulator before anything was committed".to_string(),
            ))?;
        let coin_val = latest_eval.get_commitment()?;
        channel.commit_fractal_iop_layer(*coin_val);
        match count {
            Some(count) => {
                let queries = (0..count).map(|_| channel.draw_fri_alpha()).collect();
                Ok(queries)
            }
            None => {
                let queries = (0..self.num_queries)
                    .map(|_| channel.draw_fri_alpha())
                    .collect();
                Ok(queries)
            }
        }
    }

    /// This function, implemented for the accumulator,
    /// expects as input the layer index, indexed starting at 1, since we
    /// numbered the layers that way.
    /// We'll subtract 1 from layer_idx to retrieve the actual index of the polynomial
    /// evals we are looking for.
    pub fn decommit_layer(
        &mut self,
        layer_idx: usize,
    ) -> Result<(Vec<Vec<E>>, BatchMerkleProof<H>), AccumulatorProverError> {
        // let mut coeffs_b = self.unchecked_coefficients.clone();
        // let mut coeffs_b2 = self.coefficients.clone();
        // coeffs_b.append(&mut coeffs_b2);
        // let mut multi_eval = MultiEval::<B, E, H>::new(
        //     coeffs_b,
        //     self.coefficients_ext.clone(),
        //     self.evaluation_domain_len,
        //     self.offset,
        // );
        // multi_eval.commit_polynomial_evaluations()?;

        let multi_eval =
            self.layer_evals
                .get(layer_idx - 1)
                .ok_or(AccumulatorProverError::DecommitErr(
                    layer_idx,
                    "Tried to access some strange position in the multi_evals".to_string(),
                ))?;
        let channel_state = multi_eval.get_commitment()?.clone();
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            self.public_inputs_bytes.clone(), // make sure there's actually chaining between layers
        );
        channel.commit_fractal_iop_layer(channel_state);
        let queries = channel.draw_query_positions();

        Ok(multi_eval.batch_get_values_and_proofs_at(&queries)?)
    }

    /// This is the same as decommit_layer but with queries.
    pub fn decommit_layer_with_queries(
        &self,
        layer_idx: usize,
        queries: &Vec<usize>,
    ) -> Result<(Vec<Vec<E>>, BatchMerkleProof<H>), AccumulatorProverError> {
        // let mut coeffs_b = self.unchecked_coefficients.clone();
        // let mut coeffs_b2 = self.coefficients.clone();
        // coeffs_b.append(&mut coeffs_b2);
        // let mut multi_eval = MultiEval::<B, E, H>::new(
        //     coeffs_b,
        //     self.coefficients_ext.clone(),
        //     self.evaluation_domain_len,
        //     self.offset,
        // );
        // multi_eval.commit_polynomial_evaluations()?;

        let multi_eval =
            self.layer_evals
                .get(layer_idx - 1)
                .ok_or(AccumulatorProverError::DecommitErr(
                    layer_idx,
                    "Tried to access some strange position in the multi_evals".to_string(),
                ))?;

        Ok(multi_eval.batch_get_values_and_proofs_at(queries)?)
    }

    /// This is the same as decommit_layer but with queries.
    pub fn decommit_layer_with_pub_input(
        &mut self,
        layer_idx: usize,
        pub_input: H::Digest,
    ) -> Result<(Vec<Vec<E>>, BatchMerkleProof<H>), AccumulatorProverError> {
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            self.public_inputs_bytes.clone(), // make sure there's actually chaining between layers
        );
        channel.commit_fractal_iop_layer(pub_input);
        let queries = channel.draw_query_positions();

        let multi_eval =
            self.layer_evals
                .get(layer_idx - 1)
                .ok_or(AccumulatorProverError::DecommitErr(
                    layer_idx,
                    "Tried to access some strange position in the multi_evals".to_string(),
                ))?;

        Ok(multi_eval.batch_get_values_and_proofs_at(&queries)?)
    }

    // could be named something like "finish"
    pub fn create_fri_proof(
        &mut self,
    ) -> Result<LowDegreeBatchProof<B, E, H>, AccumulatorProverError> {
        // let channel_state = self.commit_layer()?;

        let multi_eval = self
            .layer_evals
            .last()
            .ok_or(AccumulatorProverError::QueryErr(
                "You tried to query the accumulator before anything was committed".to_string(),
            ))?;
        let channel_state = *multi_eval.get_commitment()?;
        let mut channel = &mut DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            self.public_inputs_bytes.clone(),
        );

        channel.public_coin.reseed(channel_state);
        let mut low_degree_prover = LowDegreeBatchProver::<B, E, H>::new(
            &self.evaluation_domain,
            self.fri_options.clone(),
            self.max_degree,
        );

        for i in 0..self.fri_max_degrees_ext.len() {
            //println!("prover adding max_degree_ext {}", self.fri_max_degrees_ext.get(i).unwrap());
            low_degree_prover.add_polynomial_e(
                self.fri_coefficients_ext.get(i).unwrap(),
                *self.fri_max_degrees_ext.get(i).unwrap(),
                &mut channel,
            );
        }

        Ok(low_degree_prover.generate_proof(&mut channel))
    }

    /// This function takes a one-indexed layer_idx and returns the hash for that layer
    pub fn get_layer_commitment(
        &self,
        layer_idx: usize,
    ) -> Result<H::Digest, AccumulatorProverError> {
        let layer =
            self.layer_evals
                .get(layer_idx - 1)
                .ok_or(AccumulatorProverError::DecommitErr(
                    layer_idx,
                    "You tried to get a layer that doesn't exist yet.".to_string(),
                ))?;
        Ok(*layer.get_commitment()?)
    }
}

/*
pub struct FriAccumulator<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub evaluation_domain_len: usize,
    pub offset: B,
    pub evaluation_domain: Vec<B>,
    pub num_queries: usize,
    pub fri_options: FriOptions,
    pub fri_coefficients: Vec<Vec<B>>,
    pub fri_coefficients_ext: Vec<Vec<E>>,
    pub fri_max_degrees: Vec<usize>,
    pub fri_max_degrees_ext: Vec<usize>,
    _h: PhantomData<H>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > FriAccumulator<B, E, H>
{
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
            fri_coefficients: Vec::new(),
            fri_coefficients_ext: Vec::new(),
            fri_max_degrees: Vec::new(),
            fri_max_degrees_ext: Vec::new(),
            _h: PhantomData,
        }
    }

    pub fn add_polynomial(&mut self, coefficients: Vec<B>, max_degree: usize) {
        self.fri_coefficients.push(coefficients);
        self.fri_max_degrees.push(max_degree);
    }

    pub fn add_polynomial_e(&mut self, coefficients: Vec<E>, max_degree: usize) {
        self.fri_coefficients_ext.push(coefficients);
        self.fri_max_degrees_ext.push(max_degree);
    }


    // pub fn commit_layer(&mut self) -> Result<<H>::Digest, AccumulatorError> {
    //     let mut coeffs_b = self.unchecked_coefficients.clone();
    //     let mut coeffs_b2 = self.coefficients.clone();
    //     coeffs_b.append(&mut coeffs_b2);
    //     let mut multi_eval = MultiEval::<B, E, H>::new(
    //         coeffs_b,
    //         self.coefficients_ext.clone(),
    //         self.evaluation_domain_len,
    //         self.offset,
    //     );
    //     //let mut multi_eval = MultiEval::<B,E,H>::new(self.coefficients.clone(), self.coefficients_ext.clone(), self.evaluation_domain_len, self.offset);
    //     self.fri_coefficients.append(&mut self.coefficients.clone());
    //     self.fri_max_degrees.append(&mut self.max_degrees.clone());
    //     // // self.coefficients = Vec::new();
    //     // self.max_degrees = Vec::new();
    //     multi_eval.commit_polynomial_evaluations()?;
    //     Ok(multi_eval.get_commitment()?.clone())
    // }

    // pub fn draw_queries(&mut self, count: usize) -> Result<Vec<E>, AccumulatorError> {
    //     let channel_state = self.commit_layer()?;
    //     let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
    //         self.evaluation_domain_len,
    //         self.num_queries,
    //         Vec::new(), // make sure there's actually chainging between layers
    //     );
    //     channel.commit_fractal_iop_layer(channel_state);
    //     let queries = (0..count).map(|_| channel.draw_fri_alpha()).collect();
    //     Ok(queries)
    // }

    // pub fn decommit_layer(&mut self) -> Result<(Vec<Vec<E>>, BatchMerkleProof<H>), AccumulatorError> {
    //     //let mut multi_eval = MultiEval::<B,E,H>::new(self.coefficients.clone(), self.coefficients_ext.clone(), self.evaluation_domain_len, self.offset);
    //     let mut coeffs_b = self.unchecked_coefficients.clone();
    //     let mut coeffs_b2 = self.coefficients.clone();
    //     coeffs_b.append(&mut coeffs_b2);
    //     let mut multi_eval = MultiEval::<B, E, H>::new(
    //         coeffs_b,
    //         self.coefficients_ext.clone(),
    //         self.evaluation_domain_len,
    //         self.offset,
    //     );
    //     multi_eval.commit_polynomial_evaluations()?;

    //     let channel_state = multi_eval.get_commitment()?.clone();
    //     let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
    //         self.evaluation_domain_len,
    //         self.num_queries,
    //         Vec::new(), // make sure there's actually chaining between layers
    //     );
    //     channel.commit_fractal_iop_layer(channel_state);
    //     let queries = channel.draw_query_positions();
    //     println!("queries: {:?}", &queries);
    //     Ok(multi_eval.batch_get_values_and_proofs_at(queries)?)
    // }

    // could be named something like "finish"
    pub fn create_fri_proof(&mut self, pub_input_bytes: Vec<u8>) -> Result<LowDegreeBatchProof<B, E, H>, AccumulatorError> {

        let mut channel = &mut DefaultFractalProverChannel::<B, E, H>::new(
            self.evaluation_domain_len,
            self.num_queries,
            pub_input_bytes,
        );
        let mut low_degree_prover =
            LowDegreeBatchProver::<B, E, H>::new(&self.evaluation_domain, self.fri_options.clone());
        for i in 0..self.fri_max_degrees.len() {
            low_degree_prover.add_polynomial(self.fri_coefficients.get(i).unwrap(), *self.fri_max_degrees.get(i).unwrap(), &mut channel);

        }
        for i in 0..self.fri_max_degrees_ext.len() {
            low_degree_prover.add_polynomial_e(self.fri_coefficients_ext.get(i).unwrap(), *self.fri_max_degrees_ext.get(i).unwrap(), &mut channel);
        }

        Ok(low_degree_prover.generate_proof(&mut channel))
    }
}
*/

#[cfg(test)]
mod test {
    use fractal_proofs::{fields::QuadExtension, BaseElement, MultiPoly};
    use fractal_utils::polynomial_utils::MultiEval;
    use std::{convert::TryInto, marker::PhantomData, thread::AccessError};
    use winter_crypto::{hashers::Blake3_256, BatchMerkleProof, ElementHasher, MerkleTree};
    use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
    use winter_math::{fft, FieldElement, StarkField};

    use crate::errors::AccumulatorProverError;

    use super::Accumulator;
    #[test]
    fn test_accumulator() -> Result<(), AccumulatorProverError> {
        let lde_blowup = 4;
        let num_queries = 16;
        let fri_options = FriOptions::new(lde_blowup, 4, 32);
        let max_degree: usize = 63;
        let l_field_size: usize = 4 * max_degree.next_power_of_two();
        let l_field_base = BaseElement::get_root_of_unity(l_field_size.trailing_zeros());
        let offset = BaseElement::ONE;
        let evaluation_domain = winter_math::get_power_series(l_field_base, l_field_size);
        let mut acc =
            Accumulator::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>::new(
                evaluation_domain.len(),
                offset,
                evaluation_domain,
                num_queries,
                fri_options,
                vec![],
                max_degree,
            );
        acc.commit_layer()?;
        let alphas = acc.draw_queries(Some(20))?;
        assert!(alphas.len() == 20);
        Ok(())
    }
}
