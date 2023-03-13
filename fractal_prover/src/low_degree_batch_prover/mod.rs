use std::{convert::TryInto, marker::PhantomData, ops::Add};

use fractal_indexer::hash_values;
use fractal_utils::polynomial_utils::*;
use log::debug;
use winter_crypto::{BatchMerkleProof, ElementHasher, Hasher, MerkleTree};
use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
use winter_math::{fft, FieldElement, StarkField};
use winter_utils::transpose_slice;

use fractal_proofs::{
    polynom::{self, eval},
    LowDegreeBatchProof, OracleQueries,
};
use crate::channel::DefaultFractalProverChannel;


//This should be able to accumulate polynomials over time and prove at the end
pub struct LowDegreeBatchProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
> {
    // The (ongoing) random linear combination of polynomials
    randomized_sum: Vec<E>,
    // the original polynomials we're creating a batch proof over
    constituant_polynomials: Vec<Vec<E>>,
    evaluation_domain: Vec<E>,
    // The maximum allowable degrees of the polynomials you're proving about
    // Each can be any degree up to (evaluation_domain.len / fri_options.blowup_factor)-1
    max_degrees: Vec<usize>,
    // Maximum degree of the "padded" polynomial used in FRI
    // (Derived automatically by doing the opposite of how eval_domain size is derived in the winterfell fri verifier)
    fri_max_degree: usize,
    fri_options: FriOptions,
    _h: PhantomData<H>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    LowDegreeBatchProver<B, E, H>
{
    pub fn new(evaluation_domain: &Vec<B>, fri_options: FriOptions) -> Self {
        let evaluation_domain_e = evaluation_domain.iter().map(|y| E::from(*y)).collect();
        let fri_max_degree = evaluation_domain.len() / fri_options.blowup_factor() - 1;
        LowDegreeBatchProver {
            randomized_sum: Vec::new(),
            constituant_polynomials: Vec::new(),
            evaluation_domain: evaluation_domain_e,
            max_degrees: Vec::new(),
            fri_max_degree,
            fri_options,
            _h: PhantomData,
        }
    }

    //the channel is updated each time a polynomial is added so that this composes with other proofs
    pub fn add_polynomial(
        &mut self,
        polynomial_coeffs: &Vec<B>,
        max_degree: usize,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
    ) {
        let polynomial_coeffs_e: Vec<E> = polynomial_coeffs.iter().map(|y| E::from(*y)).collect();
        let alpha = channel.draw_fri_alpha();
        let beta = channel.draw_fri_alpha();
        println!("alpha: {:?}", &alpha);
        println!("beta: {:?}", &beta);
        let comp_coeffs =
            get_randomized_complementary_poly::<E>(max_degree, self.fri_max_degree, alpha, beta);

        let randomized_padded_coeffs = polynom::mul(&polynomial_coeffs_e, &comp_coeffs);
        self.randomized_sum = polynom::add(&self.randomized_sum, &randomized_padded_coeffs);
        self.max_degrees.push(max_degree);
        self.constituant_polynomials.push(polynomial_coeffs_e);
    }

    pub fn add_polynomial_e(
        &mut self,
        polynomial_coeffs: &Vec<E>,
        max_degree: usize,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
    ) {
        let alpha = channel.draw_fri_alpha();
        let beta = channel.draw_fri_alpha();
        println!("alpha: {:?}", &alpha);
        println!("beta: {:?}", &beta);
        let comp_coeffs =
            get_randomized_complementary_poly::<E>(max_degree, self.fri_max_degree, alpha, beta);

        let randomized_padded_coeffs = polynom::mul(&polynomial_coeffs, &comp_coeffs);
        self.randomized_sum = polynom::add(&self.randomized_sum, &randomized_padded_coeffs);
        self.max_degrees.push(max_degree);
        self.constituant_polynomials.push(polynomial_coeffs.clone());
    }

    pub fn generate_proof(
        &self,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
        //queried_positions: Vec<usize>,
    ) -> LowDegreeBatchProof<B, E, H> {
        // variable containing the result of evaluating each consitiuant polynomial on the set of queried eval points
        let mut all_unpadded_queried_evaluations: Vec<Vec<E>> = Vec::new();
        let mut trees: Vec<MerkleTree<H>> = Vec::new();
        let mut tree_roots: Vec<<H as Hasher>::Digest> = Vec::new();
        for poly in self.constituant_polynomials.iter() {
            let unpadded_evaluations = polynom::eval_many(&poly, &self.evaluation_domain);
            let transposed_evaluations = transpose_slice(&unpadded_evaluations);
            let hashed_evaluations = hash_values::<H, E, 1>(&transposed_evaluations);
            let tree = MerkleTree::<H>::new(hashed_evaluations).unwrap();
            tree_roots.push(*tree.root());
            trees.push(tree);
        }

        // commit to all evaluations before drawing queries
        for root in tree_roots.iter(){
            channel.commit_fri_layer(*root);
        }
        let queried_positions = channel.draw_query_positions();
        let commitment_idx = channel.layer_commitments().len();
        let eval_domain_queried = queried_positions
            .iter()
            .map(|&pos| self.evaluation_domain[pos])
            .collect::<Vec<_>>();

        for poly in self.constituant_polynomials.iter() {
            let unpadded_queried_evaluations = polynom::eval_many(&poly, &eval_domain_queried);
            all_unpadded_queried_evaluations.push(unpadded_queried_evaluations);
        }
        let tree_proofs = trees.iter().map(|tree| tree.prove_batch(&queried_positions).unwrap()).collect();
        let composed_evals: Vec<E> =
            polynom::eval_many(&self.randomized_sum, &self.evaluation_domain);
        let mut fri_prover = winter_fri::FriProver::<B, E, DefaultFractalProverChannel<B, E, H>, H>::new(
            self.fri_options.clone(),
        );
        fri_prover.build_layers(channel, composed_evals.clone());
        let fri_proof = fri_prover.build_proof(&queried_positions);
        // use only the commitments that we just added
        let commitments = channel.layer_commitments()[commitment_idx..].to_vec();
        let composed_queried_evaluations = queried_positions
            .iter()
            .map(|&p| composed_evals[p])
            .collect::<Vec<_>>();

        LowDegreeBatchProof {
            options: self.fri_options.clone(),
            num_evaluations: self.evaluation_domain.len(),
            queried_positions: queried_positions.to_vec(),
            all_unpadded_queried_evaluations,
            composed_queried_evaluations,
            commitments,
            tree_roots,
            tree_proofs,
            fri_proof: fri_proof,
            max_degrees: self.max_degrees.clone(),
            fri_max_degree: self.fri_max_degree,
        }
    }
}
