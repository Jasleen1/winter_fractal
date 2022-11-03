use std::{convert::TryInto, marker::PhantomData};

use fractal_indexer::hash_values;
use fractal_proofs::{polynom, RowcheckProof};
use fractal_utils::polynomial_utils::*;

use winter_crypto::{ElementHasher, Hasher, MerkleTree};
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use crate::{
    channel::DefaultFractalProverChannel, errors::ProverError, low_degree_prover::LowDegreeProver,
};

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

    /// The rowcheck proof generation function. Takes as input a channel and returns either an error or a rowcheck proof.
    pub fn generate_proof(
        &self,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
        initial_query_positions: Vec<usize>,
    ) -> Result<RowcheckProof<B, E, H>, ProverError> {
        // The rowcheck is supposed to prove whether f_az * f_bz - f_cz = 0 on all of H.
        // Which means that the polynomial f_az * f_bz - f_cz must be divisible by the
        // vanishing polynomial for H.
        // Since the degree of f_az and f_bz is each |H| - 1, the degree of the polynomial
        // s = (f_az * f_bz - f_cz) / vanishing_H is upper bounded by |H| - 2.
        let example_point = self.evaluation_domain[initial_query_positions[0]];

        

        // Generate coefficients for vanishing_polynomial(H)
        let denom_poly = get_vanishing_poly(self.eta, self.size_subgroup_h);
        let numerator = polynom::sub(
            &polynom::mul(&self.f_az_coeffs, &self.f_bz_coeffs),
            &self.f_cz_coeffs,
        );



        // Generate the polynomial s = (f_az * f_bz - f_cz) / vanishing_H
        let s_coeffs = polynom::div(
            &polynom::sub(
                &polynom::mul(&self.f_az_coeffs, &self.f_bz_coeffs),
                &self.f_cz_coeffs,
            ),
            &denom_poly,
        );
        println!("LHS = {:?}", polynom::mul(&s_coeffs, &denom_poly));
        println!("RHS = {:?}", numerator);
        println!("Equality = {}", polynom::mul(&s_coeffs, &denom_poly) == numerator);
        let f_az0 = polynom::eval(&self.f_az_coeffs.clone(), example_point);
        let f_bz0 = polynom::eval(&self.f_bz_coeffs.clone(), example_point);
        let f_cz0 = polynom::eval(&self.f_cz_coeffs.clone(), example_point);
        let denom_poly_0 = polynom::eval(&denom_poly.clone(), example_point);

        println!("Evals inside rowcheck = {:?}, {:?}, {:?}", f_az0, f_bz0, f_cz0);
        
        let s_poly_0 = polynom::eval(&s_coeffs.clone(), self.evaluation_domain[initial_query_positions[0]]);
        
        println!("Evals inside rowcheck prover s_computed_0 = {:?}", (f_az0 * f_bz0 - f_cz0)/denom_poly_0);

        println!("Evals inside rowcheck prover denom, s_coeffs = {:?}, {:?}", denom_poly_0, s_poly_0);

        // Build proofs for the polynomial s
        let s_prover = LowDegreeProver::<B, E, H>::from_polynomial(
            &s_coeffs,
            &self.evaluation_domain,
            self.size_subgroup_h - 1,
            self.fri_options.clone(),
        );

        let s_proof = s_prover.generate_proof(channel, initial_query_positions);
        // let old_s_evals_b: Vec<B> = polynom::eval_many(
        //     s_coeffs.clone().as_slice(),
        //     self.evaluation_domain.clone().as_slice(),
        // );
        // let old_s_evals: Vec<E> = old_s_evals_b.into_iter().map(|x: B| E::from(x)).collect();
        // let transposed_evaluations = transpose_slice(&old_s_evals);
        // let hashed_evaluations = hash_values::<H, E, 1>(&transposed_evaluations);
        // let s_tree = MerkleTree::<H>::new(hashed_evaluations)?;

        // let s_comp_coeffs =
        //     get_complementary_poly::<B>(polynom::degree_of(&s_coeffs), self.max_degree - 1);
        // let new_s = polynom::mul(&s_coeffs, &s_comp_coeffs);

        // let s_evals_b: Vec<B> = polynom::eval_many(
        //     new_s.clone().as_slice(),
        //     self.evaluation_domain.clone().as_slice(),
        // );
        // let s_evals: Vec<E> = s_evals_b.into_iter().map(|x: B| E::from(x)).collect();

        // let mut fri_prover =
        //     winter_fri::FriProver::<B, E, DefaultFractalProverChannel<B, E, H>, H>::new(
        //         self.fri_options.clone(),
        //     );

        // let query_positions = channel.draw_query_positions();
        // let queried_positions = query_positions.clone();

        // let s_eval_root = *s_tree.root();
        // let s_original_evals = query_positions
        //     .iter()
        //     .map(|&p| old_s_evals[p])
        //     .collect::<Vec<_>>();

        // let s_original_proof = s_tree.prove_batch(&queried_positions)?;
        // let commitment_idx = channel.layer_commitments().len();
        // // Build proofs for the polynomial g
        // fri_prover.build_layers(channel, s_evals.clone());
        // let s_proof = fri_prover.build_proof(&query_positions);
        // let s_queried_evals = query_positions
        //     .iter()
        //     .map(|&p| s_evals[p])
        //     .collect::<Vec<_>>();
        // let s_commitments = channel.layer_commitments()[commitment_idx..].to_vec(); //channel.layer_commitments().to_vec();
        Ok(RowcheckProof {
            options: self.fri_options.clone(),
            num_evaluations: self.evaluation_domain.len(),
            s_proof,
            s_max_degree: self.size_subgroup_h - 1,
        })
    }
}
