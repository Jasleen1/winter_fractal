use std::{marker::PhantomData, sync::Arc};

use fractal_indexer::{index::IndexParams, snark_keys::*};
use fractal_proofs::{
    fft, polynom, FractalProof, FractalProverOptions, InitialPolyProof, IopData,
    LayeredFractalProof, LayeredLincheckProof, LayeredRowcheckProof, LincheckProof,
    LowDegreeBatchProof, MultiEval, MultiPoly, TopLevelProof, TryInto,
};
use models::r1cs::Matrix;
use winter_fri::DefaultProverChannel;

use winter_crypto::{BatchMerkleProof, ElementHasher, Hasher, MerkleTree, RandomCoin};
use winter_fri::{FriOptions, ProverChannel};
use winter_math::{FieldElement, StarkField};
use winter_utils::transpose_slice;

use fractal_accumulator::accumulator::Accumulator;
use fractal_utils::channel::DefaultFractalProverChannel;

use crate::{
    errors::ProverError,
    LayeredProver, LayeredSubProver, FRACTAL_LAYERS, prover::FractalProver, CommonRSData,
};

pub struct CommonPolyVec<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    // The common polynomials in coefficient form
    pub common_polys: Vec<Vec<E>>,
}

impl<B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,> CommonRSData<B, E, H> for CommonPolyVec<B, E, H> {
    // Nothing to do here for now. 
}


pub struct LinkedR1CSAIRProver<

>