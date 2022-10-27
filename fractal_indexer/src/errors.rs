//! A list of error types which are produced during an execution of the indexing protocol

use displaydoc::Display;
use models::errors::R1CSError;
use thiserror::Error;
use winter_crypto::MerkleTreeError;

/// Represents a generic error type
#[derive(Debug, Display, Error)]
pub enum IndexerError {
    /// Error produced by the prover
    R1CS(R1CSError),
    /// If the Merkle Tree leads to an error
    MerkleTreeErr(MerkleTreeError),
}

impl From<R1CSError> for IndexerError {
    fn from(e: R1CSError) -> IndexerError {
        IndexerError::R1CS(e)
    }
}

impl From<MerkleTreeError> for IndexerError {
    fn from(e: MerkleTreeError) -> IndexerError {
        IndexerError::MerkleTreeErr(e)
    }
}
