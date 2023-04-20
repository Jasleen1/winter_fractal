mod tests;

pub use std::convert::TryInto;
use std::{fmt::Display, marker::PhantomData, usize};

pub use fractal_utils::{errors::MatrixError, matrix_utils::*, polynomial_utils::*, *};
use winter_crypto::{BatchMerkleProof, Hasher};
use winter_fri::{DefaultProverChannel, FriOptions, FriProof};
pub use winter_math::{fft, fields::f128::BaseElement, FieldElement, StarkField, *};
pub use winter_utils::{
    ByteReader, ByteWriter, Deserializable, DeserializationError, Serializable,
};

pub struct InitialPolyProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub commitment: H::Digest,
    pub evals: Vec<Vec<B>>,
    pub proof: BatchMerkleProof<H>,
    pub _phantom: PhantomData<E>,
}
pub struct FractalProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub initial_poly_proof: InitialPolyProof<B, E, H>,
    pub rowcheck_proof: RowcheckProof<B, E, H>,
    pub lincheck_a: LincheckProof<B, E, H>,
    pub lincheck_b: LincheckProof<B, E, H>,
    pub lincheck_c: LincheckProof<B, E, H>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> Serializable
    for FractalProof<B, E, H>
{
    /// Serializes `self` and writes the resulting bytes into the `target` writer.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        self.rowcheck_proof.write_into(target);
        self.lincheck_a.write_into(target);
        self.lincheck_b.write_into(target);
        self.lincheck_c.write_into(target);
    }
}

pub struct RowcheckProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub options: FriOptions,
    pub num_evaluations: usize,
    pub s_proof: LowDegreeProof<B, E, H>,
    pub s_max_degree: usize,
}

// pub struct LayeredRowcheckProof<B: StarkField, E: FieldElement<BaseField = B>> {
//     pub num_evaluations: usize,
//     pub rowcheck_poly_evals: Vec<Vec<E>>,
//     pub s_max_degree: usize,
// }

impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> Serializable
    for RowcheckProof<B, E, H>
{
    /// Serializes `self` and writes the resulting bytes into the `target` writer.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        target.write_u8(self.num_evaluations as u8);
        self.s_proof.write_into(target);
        target.write_u8(self.s_max_degree as u8);
    }
}

pub struct SumcheckProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub options: FriOptions,
    pub num_evaluations: usize,
    // Question: is it ok to use the same queried positions for both
    // g and e of different degrees?
    pub queried_positions: Vec<usize>,
    //pub g_proof: LowDegreeProof<B, E, H>,
    pub g_max_degree: usize,
    //pub e_proof: LowDegreeProof<B, E, H>,
    pub e_max_degree: usize,
    pub batch_proof: LowDegreeBatchProof<B, E, H>,
}

// TODO: FIX once interface is stable
impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> Serializable
    for SumcheckProof<B, E, H>
{
    /// Serializes `self` and writes the resulting bytes into the `target` writer.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        target.write_u8(self.num_evaluations as u8);
        target.write_u8(self.queried_positions.len() as u8);
        for pos in 0..self.queried_positions.len() {
            target.write_u8(self.queried_positions[pos] as u8);
        }
        //self.g_proof.write_into(target);
        // self.g_queried.write_into(target);
        target.write_u8(self.g_max_degree as u8);

        //self.e_proof.write_into(target);
        // self.e_queried.write_into(target);
        target.write_u8(self.e_max_degree as u8);
    }
}

pub struct LincheckProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub options: FriOptions,
    pub num_evaluations: usize,
    pub alpha: B,
    pub beta: B,
    pub t_alpha_commitment: H::Digest,
    pub t_alpha_queried: OracleQueries<B, E, H>,
    pub products_sumcheck_proof: SumcheckProof<B, E, H>,
    pub gamma: B,
    pub row_queried: OracleQueries<B, E, H>,
    pub col_queried: OracleQueries<B, E, H>,
    pub val_queried: OracleQueries<B, E, H>,
    pub matrix_sumcheck_proof: SumcheckProof<B, E, H>,
    pub _e: PhantomData<E>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> Serializable
    for LincheckProof<B, E, H>
{
    /// Serializes `self` and writes the resulting bytes into the `target` writer.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        target.write_u8(self.num_evaluations as u8);
        self.alpha.write_into(target);
        self.beta.write_into(target);
        self.t_alpha_commitment.write_into(target);
        self.t_alpha_queried.write_into(target);
        self.products_sumcheck_proof.write_into(target);
        self.gamma.write_into(target);
        self.row_queried.write_into(target);
        self.col_queried.write_into(target);
        self.val_queried.write_into(target);
        self.matrix_sumcheck_proof.write_into(target);
    }
}

pub trait LayeredProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {}

pub struct LayeredLincheckProof<B: StarkField, E: FieldElement<BaseField = B>> {
    pub row_vals: Vec<E>,
    pub col_vals: Vec<E>,
    pub val_vals: Vec<E>,
    pub f_z_vals: Vec<E>,
    pub f_mz_vals: Vec<E>,
    pub t_alpha_vals: Vec<E>,
    pub product_sumcheck_vals: Vec<(E, E)>,
    pub matrix_sumcheck_vals: Vec<(E, E)>,
    pub alpha: E,
    pub beta: E,
    pub gamma: E,
}

// impl<B: StarkField, E: FieldElement<BaseField = B>> Serializable
//     for LayeredLincheckProof<B, E>
// {
//     /// Serializes `self` and writes the resulting bytes into the `target` writer.
//     fn write_into<W: ByteWriter>(&self, target: &mut W) {
//         self.row_vals.write_into(target);
//         self.col_vals.write_into(target);
//         self.val_vals.write_into(target);
//         self.f_z_vals.write_into(target);
//         self.f_mz_vals.write_into(target);
//         self.t_alpha_vals.write_into(target);
//         self.product_sumcheck_vals.write_into(target);
//         self.matrix_sumcheck_vals.write_into(target);
//         self.alpha.write_into(target);
//         self.beta.write_into(target);
//         self.gamma.write_into(target);
//     }
// }

pub struct LayeredRowcheckProof<B: StarkField, E: FieldElement<BaseField = B>> {
    pub f_z_vals: Vec<E>,
    pub f_az_vals: Vec<E>,
    pub f_bz_vals: Vec<E>,
    pub f_cz_vals: Vec<E>,
    pub s_vals: Vec<E>,
}

pub struct LayeredSumcheckProof<B: StarkField, E: FieldElement<BaseField = B>> {
    pub numerator_vals: Vec<E>,
    pub denominator_vals: Vec<E>,
    pub sumcheck_g_vals: Vec<E>,
    pub sumcheck_e_vals: Vec<E>,
}

/*pub struct LayeredFractalProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub preprocessing_decommits_a: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub preprocessing_decommits_b: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub preprocessing_decommits_c: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub layer_commitments: [H::Digest; 3],
    pub gammas: [E; 3],
    pub layer_decommits: [(Vec<Vec<E>>, BatchMerkleProof<H>); 3],
    pub low_degree_proof: LowDegreeBatchProof<B, E, H>,
}*/

pub struct LayeredFractalProof<B: StarkField, E: FieldElement<BaseField = B>> {
    pub rowcheck: LayeredRowcheckProof<B,E>,
    pub lincheck_a: LayeredLincheckProof<B,E>,
    pub lincheck_b: LayeredLincheckProof<B,E>,
    pub lincheck_c: LayeredLincheckProof<B,E>
}

/*impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> LayeredProof<B, E, H>
    for LayeredFractalProof<B, E, H>
{
}*/

#[derive(Debug, Clone)]
pub struct QueriedPositions {
    pub queried_positions: Vec<usize>,
}
pub struct OracleQueries<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub queried_evals: Vec<E>,
    pub queried_proofs: Vec<Vec<H::Digest>>,
}

// FIXME: change this to return a Result and throw an error if qeuried_evals.len() != queried_proofs.len()
impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> OracleQueries<B, E, H> {
    pub fn new(queried_evals: Vec<E>, queried_proofs: Vec<Vec<H::Digest>>) -> Self {
        OracleQueries {
            queried_evals,
            queried_proofs,
        }
    }
}

impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> Serializable
    for OracleQueries<B, E, H>
{
    /// Serializes `self` and writes the resulting bytes into the `target` writer.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        self.queried_evals.write_into(target);
        self.queried_proofs.write_into(target);
    }
}

pub struct LowDegreeProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub options: FriOptions,
    pub num_evaluations: usize,
    pub queried_positions: Vec<usize>,
    pub unpadded_queried_evaluations: Vec<E>,
    pub padded_queried_evaluations: Vec<E>,
    pub commitments: Vec<<H>::Digest>,
    pub tree_root: H::Digest,
    pub tree_proof: BatchMerkleProof<H>,
    pub fri_proof: FriProof,
    pub max_degree: usize,
    pub fri_max_degree: usize,
}
// TODO: fix once interface is finalized (should this just be a serde macro?)
impl<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> Serializable
    for LowDegreeProof<B, E, H>
{
    /// Serializes `self` and writes the resulting bytes into the `target` writer.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        target.write_u8(self.num_evaluations as u8);
        target.write_u8(self.queried_positions.len() as u8);
        for pos in 0..self.queried_positions.len() {
            target.write_u8(self.queried_positions[pos] as u8);
        }
        self.fri_proof.write_into(target);
        //self.queried.write_into(target);
        target.write_u8(self.max_degree as u8);
    }
}
pub struct LowDegreeBatchProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher> {
    pub options: FriOptions,
    pub num_evaluations: usize,
    pub queried_positions: Vec<usize>,
    pub all_unpadded_queried_evaluations: Vec<Vec<E>>,
    pub composed_queried_evaluations: Vec<E>,
    pub commitments: Vec<<H>::Digest>,
    pub tree_roots: Vec<H::Digest>,
    pub tree_proofs: Vec<BatchMerkleProof<H>>,
    pub fri_proof: FriProof,
    pub max_degrees: Vec<usize>,
    pub fri_max_degree: usize,
}

// identifies structs with only the field elements used in an IOP. No hashes / decommitment proofs
pub trait IopData<B: StarkField, E: FieldElement<BaseField = B>> {}
impl<B: StarkField, E: FieldElement<BaseField = B>> IopData<B, E>
    for LayeredRowcheckProof<B, E>{}
impl<B: StarkField, E: FieldElement<BaseField = B>> IopData<B, E>
    for LayeredLincheckProof<B, E>{}
impl<B: StarkField, E: FieldElement<BaseField = B>> IopData<B, E>
    for LayeredSumcheckProof<B, E>{}
impl<B: StarkField, E: FieldElement<BaseField = B>> IopData<B, E>
    for LayeredFractalProof<B, E>{}

pub struct TopLevelProof<B: StarkField, E: FieldElement<BaseField = B>, H: Hasher>{
    pub preprocessing_decommitment: (Vec<Vec<E>>, BatchMerkleProof<H>),
    pub layer_commitments: Vec<H::Digest>,
    pub layer_decommitments: Vec<(Vec<Vec<E>>, BatchMerkleProof<H>)>,
    pub unverified_misc: Vec<E>,
    pub low_degree_proof: LowDegreeBatchProof<B, E, H>
}