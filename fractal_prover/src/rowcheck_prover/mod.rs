#![allow(dead_code,unused_imports)]
use std::{convert::TryInto, marker::PhantomData};

use fractal_indexer::{hash_values, index::IndexParams};
use fractal_proofs::{polynom, RowcheckProof};
use fractal_utils::{
    channel::DefaultFractalProverChannel, polynomial_utils::*, FractalProverOptions,
};

use winter_crypto::{ElementHasher, Hasher, MerkleTree};
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use fractal_accumulator::accumulator::Accumulator;
use low_degree_prover::low_degree_prover::LowDegreeProver;

use crate::{errors::ProverError, LayeredSubProver};

pub struct RowcheckProver<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    f_az_coeffs: Vec<B>,
    f_bz_coeffs: Vec<B>,
    f_cz_coeffs: Vec<B>,
    // size_subgroup_h: usize,
    // fractal_options: FractalProverOptions<B>,
    _h: PhantomData<H>,
    _e: PhantomData<E>,
    current_layer: usize,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>
    RowcheckProver<B, E, H>
{
    /// Generates a new prover for Fractal's Rowcheck operation.
    pub fn new(f_az_coeffs: Vec<B>, f_bz_coeffs: Vec<B>, f_cz_coeffs: Vec<B>) -> Self {
        RowcheckProver {
            f_az_coeffs,
            f_bz_coeffs,
            f_cz_coeffs,
            _h: PhantomData,
            _e: PhantomData,
            current_layer: 0,
        }
    }

    /// The rowcheck proof generation function. Takes as input a channel and returns either an error or a rowcheck proof.
    #[cfg_attr(feature = "flame_it", flame("rowcheck_prover"))]
    pub fn generate_proof(
        &self,
        channel: &mut DefaultFractalProverChannel<B, E, H>,
        options: &FractalProverOptions<B>,
    ) -> Result<RowcheckProof<B, E, H>, ProverError> {
        // The rowcheck is supposed to prove whether f_az * f_bz - f_cz = 0 on all of H.
        // Which means that the polynomial f_az * f_bz - f_cz must be divisible by the
        // vanishing polynomial for H.
        // Since the degree of f_az and f_bz is each |H| - 1, the degree of the polynomial
        // s = (f_az * f_bz - f_cz) / vanishing_H is upper bounded by |H| - 2.

        // Generate the polynomial s = (f_az * f_bz - f_cz) / vanishing_H
        let mut s_coeffs = polynom::sub(
            &fft_mul(&self.f_az_coeffs, &self.f_bz_coeffs),
            &self.f_cz_coeffs,
        );
        divide_by_vanishing_in_place(&mut s_coeffs, options.eta, options.h_domain.len());

        // Build proofs for the polynomial s
        let s_prover = LowDegreeProver::<B, E, H>::from_polynomial(
            &s_coeffs,
            &options.evaluation_domain,
            options.size_subgroup_h - 1,
            options.fri_options.clone(),
        );

        let s_proof = s_prover.generate_proof(channel);

        Ok(RowcheckProof {
            options: options.fri_options.clone(),
            num_evaluations: options.evaluation_domain.len(),
            s_proof,
            s_max_degree: options.size_subgroup_h - 1,
        })
    }
    #[cfg_attr(feature = "flame_it", flame("rowcheck_prover"))]
    fn rowcheck_layer_one(
        &self,
        accumulator: &mut Accumulator<B, E, H>,
        options: &FractalProverOptions<B>,
    ) {
        // The rowcheck is supposed to prove whether f_az * f_bz - f_cz = 0 on all of H.
        // Which means that the polynomial f_az * f_bz - f_cz must be divisible by the
        // vanishing polynomial for H.
        // Since the degree of f_az and f_bz is each |H| - 1, the degree of the polynomial
        // s = (f_az * f_bz - f_cz) / vanishing_H is upper bounded by |H| - 2.

        // Generate the polynomial s = (f_az * f_bz - f_cz) / vanishing_H
        let mut s_coeffs = polynom::sub(
            &fft_mul(&self.f_az_coeffs, &self.f_bz_coeffs),
            &self.f_cz_coeffs,
        );
        divide_by_vanishing_in_place(&mut s_coeffs, options.eta, options.h_domain.len());

        accumulator.add_polynomial(s_coeffs, options.size_subgroup_h - 2);
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > LayeredSubProver<B, E, H> for RowcheckProver<B, E, H>
{
    fn run_next_layer(
        &mut self,
        _query: E,
        accumulator: &mut Accumulator<B, E, H>,
        options: &FractalProverOptions<B>,
    ) -> Result<(), ProverError> {
        if self.current_layer == 0 {
            self.rowcheck_layer_one(accumulator, options);
            self.current_layer += 1;
        }
        Ok(())
    }

    fn get_num_layers(&self) -> usize {
        1
    }

    fn get_current_layer(&self) -> usize {
        self.current_layer
    }

    fn get_max_degree_constraint(
        num_input_variables: usize,
        _num_non_zero: usize,
        num_constraints: usize,
    ) -> usize {
        let h_domain_len = std::cmp::max(num_input_variables, num_constraints);
        (h_domain_len - 2).next_power_of_two()
    }
}
