//! Error types associated with representation models like R1CS.

use displaydoc::Display;
use thiserror::Error;

/// Represents errors in instantiating R1CS types
#[derive(Debug, Display, Error)]
pub enum R1CSError {
    /// Matrix should consist of a vector of equal length vectors. Not the case here.
    InvalidMatrix(String),
    /// All matrices in R1CS should have equal dimensions
    MatrixSizeMismatch(String, String),
}

/// Represents errors in instantiating input wire value vectors
#[derive(Debug, Display, Error)]
pub enum InputWireError {
    /// Generic error.
    GenericError(String),
}
