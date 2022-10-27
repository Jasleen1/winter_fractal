// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use fractal_proofs::ByteWriter;
use winter_math::StarkField;
use winter_air::{
    proof::{Commitments, Queries}, ConstraintCompositionCoefficients,
};
use core::marker::PhantomData;
use winter_crypto::{ElementHasher, RandomCoin};
use winter_fri::{self, FriProof, ProverChannel};
use winter_math::FieldElement;
use winter_utils::{collections::Vec, Serializable};
use fractal_indexer::VerifierKey;

#[cfg(feature = "concurrent")]
use utils::iterators::*;

use fractal_utils::FractalOptions;

#[derive(Debug, Copy, Clone)]
pub struct FractalContext<B, H> 
where
    B: StarkField,
    H: ElementHasher<BaseField = B>,
{    
    pub(crate) index_commitments: VerifierKey<H, B>, 
}

impl<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> Serializable for FractalContext<B, H> {
    /// Serializes `self` and writes the resulting bytes into the `target`.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        self.index_commitments.write_into(target);
    }
}

// TODO: Make a FractalCommitments type.


// TYPES AND INTERFACES
// ================================================================================================

pub struct FractalProverChannel<B, E, H>
where
    B: StarkField, 
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
{
    // Next TODO: Add context here to be able to see how many queries and eval points etc.
    context: FractalContext<B, H>,
    pub public_coin: RandomCoin<B, H>,
    pub commitments: Commitments,
    pow_nonce: u64,
    _field_element: PhantomData<E>,
}

// PROVER CHANNEL IMPLEMENTATION
// ================================================================================================

impl<B, E, H> FractalProverChannel<B, E, H>
where
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
{
    // CONSTRUCTOR
    // --------------------------------------------------------------------------------------------
    /// Creates a new prover channel for the specified `air` and public inputs.
    pub fn new(context: FractalContext<B, H>, pub_inputs_bytes: Vec<u8>) -> Self {

        // build a seed for the public coin; the initial seed is the hash of public inputs but as the protocol progresses, the coin will be reseeded with the info sent to
        // the verifier
        let mut coin_seed = pub_inputs_bytes;
        context.write_into(&mut coin_seed);

        FractalProverChannel {
            context,
            public_coin: RandomCoin::new(&coin_seed),
            commitments: Commitments::default(),
            pow_nonce: 0,
            _field_element: PhantomData,
        }
    }

    // COMMITMENT METHODS
    // --------------------------------------------------------------------------------------------


    /// Commits the prover to the evaluations of the polynomials f_z, f_az, f_bz.
    pub fn commit_initial_witness_evals(&mut self, initial_witness_root: H::Digest) {
        self.commitments.add::<H>(&initial_witness_root);
        self.public_coin.reseed(initial_witness_root);
    }

    // Commits the prover to the evaluations of the polynomials generated during 
    // lincheck and rowcheck
    pub fn commit_initial_polynomials(&mut self, initial_poly_com_root: H::Digest) {
        self.commitments.add::<H>(&initial_poly_com_root);
        self.public_coin.reseed(initial_poly_com_root);
    }

    // Commits to the evaluations of the low degree polynomial 
    pub fn commit_low_degree_poly(&mut self, low_degree_poly_com_root: H::Digest) {
        self.commitments.add::<H>(&low_degree_poly_com_root);
        self.public_coin.reseed(low_degree_poly_com_root);
    }
    

    // PUBLIC COIN METHODS
    // --------------------------------------------------------------------------------------------

    /// Returns a set of coefficients for TODO.
    ///
    /// The coefficients are drawn from the public coin uniformly at random.
    pub fn get_composition_coeffs(&mut self) -> ConstraintCompositionCoefficients<E> {
        // self.air
        //     .get_constraint_composition_coefficients(&mut self.public_coin)
        //     .expect("failed to draw composition coefficients")
        unimplemented!()
        // TODO we need some version of this to get the coeffs of the composition polynomials for 
        // an R1CS instance, depending on the fractal_options
    }

    /// Returns an out-of-domain point drawn uniformly at random from the public coin. 
    /// This is used in lincheck to get the beta value.
    /// Beta should be selected from Field \ h_domain, and we know h_domain is a multiplicative
    /// coset, so we just check that before returning.
    pub fn get_lincheck_beta(&mut self) -> E {
        let mut drawn_pt: E = self.public_coin.draw().expect("failed to draw OOD point");
        let eta = E::from(self.context.index_commitments.params.eta);
        let h_domain_size = self.context.index_commitments.params.num_input_variables;
        let h_domain_size_exp = E::PositiveInteger::from(h_domain_size as u64);
        while  (drawn_pt.div(eta)).exp(h_domain_size_exp) == E::ONE {
            self.public_coin.reseed(H::hash(&drawn_pt.to_bytes()));
            drawn_pt = self.public_coin.draw().expect("failed to draw OOD point");
        }
        drawn_pt
    }

    /// Returns an out-of-domain point drawn uniformly at random from the public coin. 
    /// This is used in lincheck to get the alpha value.
    /// Alpha should be selected from Field \ h_domain, and we know h_domain is a multiplicative
    /// coset, so we just check that before returning.
    pub fn get_lincheck_alpha(&mut self) -> E {
        let mut drawn_pt: E = self.public_coin.draw().expect("failed to draw OOD point");
        let eta = E::from(self.context.index_commitments.params.eta);
        let h_domain_size = self.context.index_commitments.params.num_input_variables;
        let h_domain_size_exp = E::PositiveInteger::from(h_domain_size as u64);
        while  (drawn_pt.div(eta)).exp(h_domain_size_exp) == E::ONE {
            self.public_coin.reseed(H::hash(&drawn_pt.to_bytes()));
            drawn_pt = self.public_coin.draw().expect("failed to draw OOD point");
        }
        drawn_pt
    }


    /// Returns a set of positions in the LDE domain against which the evaluations of trace and
    /// constraint composition polynomials should be queried.
    ///
    /// The positions are drawn from the public coin uniformly at random.
    pub fn get_query_positions(&mut self) -> Vec<usize> {
        let num_queries = self.context.index_commitments.params.num_queries;
        let eval_domain_size = self.context.index_commitments.params.blowup_factor * self.context.index_commitments.params.max_degree;
        self.public_coin
            .draw_integers(num_queries, eval_domain_size)
            .expect("failed to draw query position")
    }


    // /// Determines a nonce, which when hashed with the current seed of the public coin results
    // /// in a new seed with the number of leading zeros equal to the grinding_factor specified
    // /// in the proof options.
    // pub fn grind_query_seed(&mut self) {
    //     let grinding_factor = self.context.options().grinding_factor();

    //     #[cfg(not(feature = "concurrent"))]
    //     let nonce = (1..u64::MAX)
    //         .find(|&nonce| self.public_coin.check_leading_zeros(nonce) >= grinding_factor)
    //         .expect("nonce not found");

    //     #[cfg(feature = "concurrent")]
    //     let nonce = (1..u64::MAX)
    //         .into_par_iter()
    //         .find_any(|&nonce| self.public_coin.check_leading_zeros(nonce) >= grinding_factor)
    //         .expect("nonce not found");

    //     self.pow_nonce = nonce;
    //     self.public_coin.reseed_with_int(nonce);
    // }

    // PROOF BUILDER
    // --------------------------------------------------------------------------------------------
    /// Builds a proof from the previously committed values as well as values passed into
    /// this method.
    pub fn build_proof(
        self,
        trace_queries: Vec<Queries>,
        constraint_queries: Queries,
        fri_proof: FriProof,
    )  {
        // We may not need this function at all
        unimplemented!()
    }
}

// FRI PROVER CHANNEL IMPLEMENTATION
// ================================================================================================

impl<B, E, H> ProverChannel<E> for FractalProverChannel<B, E, H>
where
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
{
    type Hasher = H;

    /// Commits the prover to a FRI layer.
    fn commit_fri_layer(&mut self, layer_root: H::Digest) {
        self.commitments.add::<H>(&layer_root);
        self.public_coin.reseed(layer_root);
    }

    /// Returns a new alpha drawn from the public coin.
    fn draw_fri_alpha(&mut self) -> E {
        self.public_coin.draw().expect("failed to draw FRI alpha")
    }
}
