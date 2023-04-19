use displaydoc::Display;
use thiserror::Error;
use fractal_utils::errors::FractalUtilError;
use winter_crypto::MerkleTreeError;
/// Represents a generic error type
#[derive(Debug, Display, Error, PartialEq)]
pub enum AccumulatorProverError {
    /// If the accumulator's decommit leads to an error
    DecommitErr(usize, String),
    /// Merkle tree error within the accumulator
    MerkleTreeErr(MerkleTreeError),
    /// Util Error
    FractalUtilErr(FractalUtilError),
    /// If the caller tries to operate on an accumulator which doesn't yet have commitments.
    QueryErr(String),
}
impl From<MerkleTreeError> for AccumulatorProverError {
    fn from(e: MerkleTreeError) -> AccumulatorProverError {
        AccumulatorProverError::MerkleTreeErr(e)
    }
}

impl From<FractalUtilError> for AccumulatorProverError {
    fn from(e: FractalUtilError) -> AccumulatorProverError {
        AccumulatorProverError::FractalUtilErr(e)
    }
}