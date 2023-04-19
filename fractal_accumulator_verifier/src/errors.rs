use fractal_utils::errors::FractalUtilError;
use low_degree_verifier::errors::LowDegreeVerifierError;
use thiserror::Error;
use winter_crypto::{MerkleTreeError, RandomCoinError};

/// Represents a generic error type
#[derive(Debug, Error, PartialEq)]
pub enum AccumulatorVerifierError {
    /// If the accumulator's decommit leads to an error
    QueryVerificationErr(String),
    /// Merkle tree error within the accumulator
    MerkleTreeErr(MerkleTreeError),
    /// Util Error
    FractalUtilErr(FractalUtilError),
    /// If the caller tries to operate on an accumulator which doesn't yet have commitments.
    QueryErr(String),
    /// Claimed root and layer commit don't match
    CommitMatchErr(String),
    /// Low degree verification error
    LowDegreeVerifierErr(LowDegreeVerifierError),
    /// Random coin error
    RandomCoinErr(RandomCoinError),
}

impl From<MerkleTreeError> for AccumulatorVerifierError {
    fn from(e: MerkleTreeError) -> AccumulatorVerifierError {
        AccumulatorVerifierError::MerkleTreeErr(e)
    }
}

impl From<FractalUtilError> for AccumulatorVerifierError {
    fn from(e: FractalUtilError) -> AccumulatorVerifierError {
        AccumulatorVerifierError::FractalUtilErr(e)
    }
}

impl From<LowDegreeVerifierError> for AccumulatorVerifierError {
    fn from(e: LowDegreeVerifierError) -> AccumulatorVerifierError {
        AccumulatorVerifierError::LowDegreeVerifierErr(e)
    }
}

impl From<RandomCoinError> for AccumulatorVerifierError {
    fn from(e: RandomCoinError) -> AccumulatorVerifierError {
        AccumulatorVerifierError::RandomCoinErr(e)
    }
}

impl std::fmt::Display for AccumulatorVerifierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            AccumulatorVerifierError::QueryVerificationErr(err) => {
                writeln!(f, "Accumulator verification error from queries: {}", err)
            }
            AccumulatorVerifierError::MerkleTreeErr(err) => {
                writeln!(f, "Fractal accumulator verifier Merkle Tree error: {}", err)
            }
            AccumulatorVerifierError::FractalUtilErr(err) => {
                writeln!(f, "Fractal accumulator util error: {}", err)
            }
            AccumulatorVerifierError::QueryErr(err) => {
                writeln!(f, "Problem with query in the accumulator: {}", err)
            }
            AccumulatorVerifierError::CommitMatchErr(err) => {
                writeln!(
                    f,
                    "The commitment input here didn't match that derived from the proof: {}",
                    err
                )
            }
            AccumulatorVerifierError::LowDegreeVerifierErr(err) => {
                writeln!(f, "The low degree proof didn't verify: {}", err)
            }
            AccumulatorVerifierError::RandomCoinErr(err) => {
                writeln!(f, "Problem with the random coin: {}", err)
            }
        }
    }
}