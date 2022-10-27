//! A list of error types which are produced during an execution of the protocol

use core::fmt;

use displaydoc::Display;
use models::errors::R1CSError;
use thiserror::Error;
use winter_crypto::MerkleTreeError;

#[derive(Debug, Error)]
pub enum ProverError {
    LincheckErr(LincheckError),
    R1CSErr(R1CSError),
    InvalidMatrixName(String),
    MerkleTreeErr(MerkleTreeError),
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

/// Represents a generic error type
#[derive(Debug, Display, Error)]
pub enum LincheckError {
    /// If the Merkle Tree leads to an error
    MerkleTreeErr(MerkleTreeError),
}

impl From<MerkleTreeError> for LincheckError {
    fn from(e: MerkleTreeError) -> LincheckError {
        LincheckError::MerkleTreeErr(e)
    }
}

// impl fmt::Display for LincheckError {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         match self {
//             Self::MerkleTreeErr(err) => {
//                 write!(f, "Encountered an error in Lincheck: {:?}", err,)
//             }
//         }
//     }
// }

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
        }
    }
}
