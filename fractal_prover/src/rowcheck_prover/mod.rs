use std::{convert::TryInto, marker::PhantomData};

use fractal_indexer::hash_values;
use fractal_utils::polynomial_utils::*;
use fractal_proofs::{RowcheckProof, polynom};

use winter_crypto::{ElementHasher, Hasher, MerkleTree};
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use crate::{errors::ProverError, prover_channel::FractalProverChannel};

pub struct RowcheckProver<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    f_az_coeffs: Vec<B>,
    f_bz_coeffs: Vec<B>,
    f_cz_coeffs: Vec<B>,
    degree_fs: usize,
    size_subgroup_h: usize,
    evaluation_domain: Vec<B>,
    fri_options: FriOptions,
    num_queries: usize,
    max_degree: usize,
    eta: B,
    _h: PhantomData<H>,
    _e: PhantomData<E>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    RowcheckProver<B, E, H>
{
    pub fn new(
        f_az_coeffs: Vec<B>,
        f_bz_coeffs: Vec<B>,
        f_cz_coeffs: Vec<B>,
        degree_fs: usize,
        size_subgroup_h: usize,
        evaluation_domain: Vec<B>,
        fri_options: FriOptions,
        num_queries: usize,
        max_degree: usize,
        eta: B,
    ) -> Self {
        RowcheckProver {
            f_az_coeffs,
            f_bz_coeffs,
            f_cz_coeffs,
            degree_fs,
            size_subgroup_h,
            evaluation_domain,
            fri_options,
            num_queries,
            max_degree,
            eta,
            _h: PhantomData,
            _e: PhantomData,
        }
    }

    pub fn generate_proof(&self, channel: &mut  FractalProverChannel<B, E, H>) -> Result<RowcheckProof<B, E, H>, ProverError> {
        let mut denom_poly = vec![B::ZERO; self.size_subgroup_h-1];
        denom_poly.push(B::ONE);
        let h_size_32: u32 = self.size_subgroup_h.try_into().unwrap();
        let eta_pow = B::PositiveInteger::from(h_size_32);
        denom_poly[0] = self.eta.exp(eta_pow).neg();
        let s_coeffs = polynom::div(
            &polynom::sub(&polynom::mul(&self.f_az_coeffs, &self.f_bz_coeffs), &self.f_cz_coeffs),
            &denom_poly
        );   
        let old_s_evals_b: Vec<B> = polynom::eval_many(s_coeffs.clone().as_slice(), self.evaluation_domain.clone().as_slice());// Vec::new();
        let old_s_evals: Vec<E> = old_s_evals_b.into_iter().map(|x: B| {E::from(x)}).collect();
        let transposed_evaluations = transpose_slice(&old_s_evals);
        let hashed_evaluations = hash_values::<H, E, 1>(&transposed_evaluations);
        let s_tree = MerkleTree::<H>::new(hashed_evaluations)?;
        
        let s_comp_coeffs = get_complementary_poly::<B>(polynom::degree_of(&s_coeffs), self.max_degree - 1);
        let new_s = polynom::mul(&s_coeffs, &s_comp_coeffs);

        let s_evals_b: Vec<B> = polynom::eval_many(new_s.clone().as_slice(), self.evaluation_domain.clone().as_slice());// Vec::new();
        let s_evals: Vec<E> = s_evals_b.into_iter().map(|x: B| {E::from(x)}).collect();
        

        let mut fri_prover =
            winter_fri::FriProver::<B, E, FractalProverChannel<B, E, H>, H>::new(self.fri_options.clone());

        let query_positions = channel.get_query_positions();
        let queried_positions = query_positions.clone();

        let s_eval_root = *s_tree.root();
        let s_original_evals = query_positions
            .iter()
            .map(|&p| old_s_evals[p])
            .collect::<Vec<_>>();
        
        let s_original_proof = s_tree.prove_batch(&queried_positions)?;

        // Build proofs for the polynomial g
        fri_prover.build_layers(channel, s_evals.clone());
        let s_proof = fri_prover.build_proof(&query_positions);
        let s_queried_evals = query_positions
            .iter()
            .map(|&p| s_evals[p])
            .collect::<Vec<_>>();
        // let s_commitments = channel.layer_commitments().to_vec();
        Ok(RowcheckProof {
            options: self.fri_options.clone(),
            num_evaluations: self.evaluation_domain.len(),
            queried_positions,
            s_eval_root,
            s_original_evals,
            s_original_proof,
            s_proof,
            s_queried_evals,
            // s_commitments,
            s_max_degree: self.size_subgroup_h - 2,
        })
    }
}
