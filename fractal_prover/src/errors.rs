//! A list of error types which are produced during an execution of the protocol

use core::fmt;

use displaydoc::Display;
use fractal_proofs::errors::FractalUtilError;
use models::errors::R1CSError;
use thiserror::Error;
use winter_crypto::MerkleTreeError;

#[derive(Debug, Error, PartialEq)]
pub enum ProverError {
    LincheckErr(LincheckError),
    R1CSErr(R1CSError),
    InvalidMatrixName(String),
    MerkleTreeErr(MerkleTreeError),
    MultiPolyErr(String),
    FractalUtilErr(FractalUtilError),
    AccumulatorErr(AccumulatorProverError),
    ProverKeyNoneErr(),
}

impl From<LincheckError> for ProverError {
    fn from(e: LincheckError) -> ProverError {
        ProverError::LincheckErr(e)
    }
}

impl From<R1CSError> for ProverError {
    fn from(e: R1CSError) -> ProverError {
        ProverError::R1CSErr(e)
    }
}

impl From<MerkleTreeError> for ProverError {
    fn from(e: MerkleTreeError) -> ProverError {
        ProverError::MerkleTreeErr(e)
    }
}

impl From<FractalUtilError> for ProverError {
    fn from(e: FractalUtilError) -> ProverError {
        ProverError::FractalUtilErr(e)
    }
}

impl From<AccumulatorProverError> for ProverError {
    fn from(e: AccumulatorProverError) -> ProverError {
        ProverError::AccumulatorErr(e)
    }
}

/// Represents a generic error type
#[derive(Debug, PartialEq, Error)]
pub enum LincheckError {
    /// If the Merkle Tree leads to an error
    MerkleTreeErr(MerkleTreeError),
    /// If you tried to compute gamma without having set alpha or t_alpha
    GammaCompErr(String),
}

impl fmt::Display for LincheckError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MerkleTreeErr(err) => {
                write!(f, "Encountered an error in Lincheck: {:?}", err,)
            }
            Self::GammaCompErr(err) => {
                write!(
                    f,
                    "Encountered an error in Lincheck, you tried to compute gamma: {:?}",
                    err,
                )
            }
        }
    }
}

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

impl From<MerkleTreeError> for LincheckError {
    fn from(e: MerkleTreeError) -> LincheckError {
        LincheckError::MerkleTreeErr(e)
    }
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

impl fmt::Display for ProverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMatrixName(matrix_name) => {
                write!(f, "Invalid matrix name for multiplying {}", matrix_name)
            }
            Self::LincheckErr(err) => {
                write!(f, "Encountered an error in Lincheck: {:?}", err,)
            }
            Self::R1CSErr(err) => {
                write!(
                    f,
                    "Encountered an R1CS error in the fractal prover: {:?}",
                    err,
                )
            }
            Self::MerkleTreeErr(err) => {
                write!(
                    f,
                    "Encountered a Merkle Tree error in the fractal prover: {:?}",
                    err,
                )
            }
            Self::FractalUtilErr(err) => {
                write!(
                    f,
                    "Encountered an error using utils in the fractal prover: {:?}",
                    err,
                )
            }
            Self::MultiPolyErr(err) => {
                write!(
                    f,
                    "Encountered an error while trying to deal with commitments of multiple polynomials the fractal prover: {:?}",
                    err,
                )
            }
            Self::AccumulatorErr(err) => {
                write!(
                    f,
                    "Encountered an error in the accumulator somewhere in the fractal prover: {:?}",
                    err,
                )
            }
            Self::ProverKeyNoneErr() => {
                write!(
                    f,
                    "Encountered an error in the proof generation: you tried to unwrap a None ProverKey"
                )
            }
        }
    }
}
