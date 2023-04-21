//! A list of error types which are produced during an execution of the protocol

use core::fmt;

use displaydoc::Display;
use fractal_accumulator::errors::AccumulatorProverError;
use fractal_proofs::errors::FractalUtilError;
use models::errors::R1CSError;
use thiserror::Error;
use winter_crypto::MerkleTreeError;


/// The errors for a Fractal Prover
#[derive(Debug, Error, PartialEq)]
pub enum ProverError {
    /// Error handling for errors in a [`crate::lincheck_prover::LincheckProver`]
    LincheckErr(LincheckError),
    /// Error handling for R1CS data structure related errors
    R1CSErr(R1CSError),
    /// Used in testing sometimes, if the matrix name provided is not valid.
    InvalidMatrixName(String),
    /// Error related to Merkle Tree operations
    MerkleTreeErr(MerkleTreeError),
    /// Error related to the [`fractal_utils::polynomial_utils::MultiEval`] structs
    MultiPolyErr(String),
    /// Other errors related to [`fractal_utils`]
    FractalUtilErr(FractalUtilError),
    /// Errors related to the [`fractal_accumulator`] crate.
    AccumulatorErr(AccumulatorProverError),
    /// In some cases, a prover key for a struct my be an option and may not be set. 
    /// Logically speaking it shouldn't be accessed in such a situation. 
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

/// Represents a generic error type for lincheck-related operations.
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

impl From<MerkleTreeError> for LincheckError {
    fn from(e: MerkleTreeError) -> LincheckError {
        LincheckError::MerkleTreeErr(e)
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
