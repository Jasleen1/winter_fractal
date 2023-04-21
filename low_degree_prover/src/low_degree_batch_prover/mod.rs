use std::{convert::TryInto, marker::PhantomData, ops::Add};

use fractal_utils::polynomial_utils::*;
use log::debug;
use winter_crypto::{BatchMerkleProof, ElementHasher, Hasher, MerkleTree};
use winter_fri::utils::hash_values;
use winter_fri::{DefaultProverChannel, FriOptions, ProverChannel};
use winter_math::{fft, FieldElement, StarkField};
use winter_utils::transpose_slice;

use fractal_proofs::{
    polynom::{self, eval},
    LowDegreeBatchProof, OracleQueries,
};
use fractal_utils::channel::DefaultFractalProverChannel;

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
    /// Creates a new low degree batch prover
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

    /// Adds a polynomial to the low degree batch prover.
    #[cfg_attr(feature = "flame_it", flame("low_degree_prover"))]
    pub fn add_polynomial(
        &mut self,
        polynomial_coeffs: &Vec<B>,
        max_degree: usize,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
    ) {
        //the channel is updated each time a polynomial is added so that this composes with other proofs

        let polynomial_coeffs_e: Vec<E> = polynomial_coeffs.iter().map(|y| E::from(*y)).collect();
        self.add_polynomial_e(&polynomial_coeffs_e, max_degree, channel);
    }

    /// Adds a polynomial with coefficients in the extension field.
    #[cfg_attr(feature = "flame_it", flame("low_degree_prover"))]
    pub fn add_polynomial_e(
        &mut self,
        polynomial_coeffs: &Vec<E>,
        max_degree: usize,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
    ) {
        let alpha = channel.draw_fri_alpha();
        let beta = channel.draw_fri_alpha();

        let comp_coeffs =
            get_randomized_complementary_poly::<E>(max_degree, self.fri_max_degree, alpha, beta);

        // easy multiplication, don't use fft_mul here
        let randomized_padded_coeffs = polynom::mul(&polynomial_coeffs, &comp_coeffs);
        self.randomized_sum = polynom::add(&self.randomized_sum, &randomized_padded_coeffs);
        self.max_degrees.push(max_degree);
        self.constituant_polynomials.push(polynomial_coeffs.clone());
    }

    /// Helper function to zip the evaluations so that each element of the output is of the
    /// form [poly_1(e), ..., poly_n(e)] i.e. evaluations of all the polynomials are included
    /// in the same array.
    #[cfg_attr(feature = "flame_it", flame("low_degree_prover"))]
    fn zip_evals(separate_evals: Vec<Vec<E>>, evaluation_domain_len: usize) -> Vec<Vec<E>> {
        let mut zipped_evals = vec![Vec::<E>::new(); evaluation_domain_len];
        for (_, eval) in separate_evals.iter().enumerate() {
            for (loc, &val) in eval.iter().enumerate() {
                zipped_evals[loc].push(val);
            }
        }
        zipped_evals
    }

    #[cfg_attr(feature = "flame_it", flame("low_degree_prover"))]
    /// Generates a low degree proof for a batch of values.
    pub fn generate_proof(
        &self,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
    ) -> LowDegreeBatchProof<B, E, H> {
        // variable containing the result of evaluating each consitiuant polynomial on the set of queried eval points
        let mut all_unpadded_queried_evaluations: Vec<Vec<E>> = Vec::new();

        let mut all_unpadded_evaluations = vec![];
        let eval_domain_size = self.evaluation_domain.len();
        let eval_domain_twiddles: Vec<B> = fft::get_twiddles(eval_domain_size);

        flame::start("loop1");
        for poly in self.constituant_polynomials.iter() {
            let mut unpadded_evals = poly.clone();
            pad_with_zeroes(&mut unpadded_evals, eval_domain_size);
            fft::evaluate_poly(&mut unpadded_evals, &eval_domain_twiddles);
            all_unpadded_evaluations.push(unpadded_evals);
        }
        flame::end("loop1");

        flame::start("make tree");
        let zipped_evals = Self::zip_evals(
            all_unpadded_evaluations.clone(),
            self.evaluation_domain.len(),
        );
        let eval_hashes = zipped_evals
            .iter()
            .map(|evals| H::hash_elements(evals))
            .collect::<Vec<_>>();
        let tree = MerkleTree::<H>::new(eval_hashes).unwrap();
        let tree_root = *tree.root();
        flame::end("make tree");

        flame::start("commit_fri_layer");
        channel.commit_fri_layer(tree_root);
        flame::end("commit_fri_layer");

        let queried_positions = channel.draw_query_positions();

        flame::start("tree_proof");
        let tree_proof = tree.prove_batch(&queried_positions).unwrap();
        flame::end("tree_proof");

        let commitment_idx = channel.layer_commitments().len();

        for evals in all_unpadded_evaluations {
            let unpadded_queried_evaluations = queried_positions
                .iter()
                .map(|&pos| evals[pos])
                .collect::<Vec<_>>();
            all_unpadded_queried_evaluations.push(unpadded_queried_evaluations);
        }

        flame::start("composed_evals");
        let composed_evals: Vec<E> =
            polynom::eval_many(&self.randomized_sum, &self.evaluation_domain);
        flame::end("composed_evals");

        let mut fri_prover =
            winter_fri::FriProver::<B, E, DefaultFractalProverChannel<B, E, H>, H>::new(
                self.fri_options.clone(),
            );
        flame::start("fri_proof");
        fri_prover.build_layers(channel, composed_evals.clone());
        let fri_proof = fri_prover.build_proof(&queried_positions);
        flame::end("fri_proof");
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
            tree_root,
            tree_proof,
            fri_proof: fri_proof,
            max_degrees: self.max_degrees.clone(),
            fri_max_degree: self.fri_max_degree,
        }
    }
}
