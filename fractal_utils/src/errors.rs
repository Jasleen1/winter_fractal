//! A list of error types which are produced during an execution of the protocol

use displaydoc::Display;
use thiserror::Error;

/// Represents a generic error type
#[derive(Debug, Display, Error)]
pub enum FractalUtilError {
    /// Error produced by the prover
    MATRIX(MatrixError),
}

impl From<MatrixError> for FractalUtilError {
    fn from(e: MatrixError) -> FractalUtilError {
        FractalUtilError::MATRIX(e)
    }
}

/// Represents errors in instantiating R1CS types
#[derive(Debug, Display, Error)]
pub enum MatrixError {
    /// Matrix should consist of a vector of equal length vectors. Not the case here.
    InvalidMatrix(String),
    /// All matrices in R1CS should have equal dimensions
    MatrixSizeMismatch(String, String),
    /// Number of cols in the first matrix should equal the number of rows in the second.
    MatrixMultiplicationSizeErr(String, String),
}
