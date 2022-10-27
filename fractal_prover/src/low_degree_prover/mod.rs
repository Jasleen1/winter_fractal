use std::{convert::TryInto, marker::PhantomData};

use fractal_utils::{polynomial_utils::*, FractalOptions};
use crate::prover_channel::FractalProverChannel;
use winter_crypto::{ElementHasher, Hasher, MerkleTree};
use winter_fri::{FriOptions};
use winter_math::{fft, FieldElement, StarkField};
use winter_utils::{transpose_slice};
use fractal_indexer::hash_values;



use fractal_proofs::{OracleQueries, LowDegreeProof, polynom::{self, eval}};

pub struct LowDegreeProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
> {
    // The polynomial we're proving has a low degree
    polynomial_coeffs: Vec<E>,
    // evaluations of the polynomial over the Evaluation Domain
    polynomial_evals: Vec<E>,
    evaluation_domain: Vec<E>,
    // The maximum allowable degree of the polynomial you're proving about
    // This can be any degree up to (evaluation_domain.len / fri_options.blowup_factor)-1
    max_degree: usize,
    // Maximum degree of the "padded" polynomial used in FRI
    // (Derived automatically by doing the opposite of how eval_domain size is derived in the winterfell fri verifier)
    fri_max_degree: usize,
    fri_options: FriOptions,
    _h: PhantomData<H>
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>,>
    LowDegreeProver<B, E, H>
{
    pub fn from_polynomial(
        polynomial: &Vec<B>,
        evaluation_domain: &Vec<B>,
        max_degree: usize,
        fri_options: FriOptions,
    ) -> Self {
        let polynomial_evals = polynom::eval_many(&polynomial, &evaluation_domain).iter().map(|x| E::from(*x)).collect();
        let polynomial_e = polynomial.iter().map(|c| E::from(*c)).collect();
        let fri_max_degree = evaluation_domain.len() / fri_options.blowup_factor() -1;
        assert!(polynom::degree_of(&polynomial) <= max_degree);
        let evaluation_domain_e = evaluation_domain.iter().map(|y| E::from(*y)).collect();
        LowDegreeProver {
            polynomial_coeffs: polynomial_e,
            polynomial_evals,
            evaluation_domain: evaluation_domain_e,
            max_degree,
            fri_max_degree,
            fri_options,
            _h: PhantomData
        }
    }

    pub fn from_evals(
        polynomial_evals: Vec<E>,
        evaluation_domain: &Vec<E>,
        max_degree: usize,
        fri_options: FriOptions,
    ) -> Self {
        assert_eq!(polynomial_evals.len(), evaluation_domain.len());
        let polynomial_coeffs = polynom::interpolate(&evaluation_domain, &polynomial_evals, true);
        assert!(polynom::degree_of(&polynomial_coeffs) <= max_degree);
        let fri_max_degree = evaluation_domain.len() / fri_options.blowup_factor() -1;
        LowDegreeProver {
            polynomial_coeffs,
            polynomial_evals,
            evaluation_domain: evaluation_domain.clone(),
            max_degree,
            fri_max_degree,
            fri_options,
            _h: PhantomData
        }
    }

    pub fn generate_proof(&self, channel: &mut FractalProverChannel<B, E, H>) -> LowDegreeProof<B, E, H> {
        let queried_positions = channel.get_query_positions();
        // let commitment_idx = channel.commitments.0.len();
        let unpadded_queried_evaluations = queried_positions
            .iter()
            .map(|&p| self.polynomial_evals[p])
            .collect::<Vec<_>>();

        let transposed_evaluations = transpose_slice(&self.polynomial_evals);
        let hashed_evaluations = hash_values::<H, E, 1>(&transposed_evaluations);
        let tree = MerkleTree::<H>::new(hashed_evaluations).unwrap();
        let tree_root = *tree.root();
        let tree_proof = tree.prove_batch(&queried_positions).unwrap();

        let comp_coeffs = get_complementary_poly::<E>(self.max_degree, self.fri_max_degree);
        let padded_coeffs = polynom::mul(&self.polynomial_coeffs, &comp_coeffs);
        let padded_evals: Vec<E> = polynom::eval_many(&padded_coeffs, &self.evaluation_domain);

        let mut fri_prover =
            winter_fri::FriProver::<B, E, FractalProverChannel<B, E, H>, H>::new(self.fri_options.clone());
        fri_prover.build_layers(channel, padded_evals.clone());
        let fri_proof = fri_prover.build_proof(&queried_positions);
        // use only the commitments that we just added
        // let commitments = channel.commitments[commitment_idx..].to_vec();
        let padded_queried_evaluations = queried_positions
            .iter()
            .map(|&p| padded_evals[p])
            .collect::<Vec<_>>();

        LowDegreeProof {
            options: self.fri_options.clone(),
            num_evaluations: self.polynomial_evals.len(),
            queried_positions: queried_positions.to_vec(),
            unpadded_queried_evaluations,
            padded_queried_evaluations,
            // commitments,
            tree_root,
            tree_proof,
            fri_proof: fri_proof,
            max_degree: self.max_degree,
            fri_max_degree: self.fri_max_degree,
        }
    }
}