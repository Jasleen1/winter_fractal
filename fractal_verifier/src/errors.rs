// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under both the MIT license found in the
// LICENSE-MIT file in the root directory of this source tree and the Apache
// License, Version 2.0 found in the LICENSE-APACHE file in the root directory
// of this source tree.

//! Errors for various data structure operations.
use fractal_proofs::DeserializationError;
use fractal_sumcheck::errors::SumcheckVerifierError;
use winter_fri::VerifierError;

#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub enum LincheckVerifierError {
    /// Error propagation
    UnsoundProduct(SumcheckVerifierError),
    /// Error propagation
    UnsoundMatrix(SumcheckVerifierError),
}

impl From<SumcheckVerifierError> for LincheckVerifierError {
    fn from(error: SumcheckVerifierError) -> Self {
        Self::UnsoundProduct(error)
    }
}

impl std::fmt::Display for LincheckVerifierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            LincheckVerifierError::UnsoundProduct(err) => {
                writeln!(f, "Lincheck error: unsound product: {}", err)
            }
            LincheckVerifierError::UnsoundMatrix(err) => {
                writeln!(f, "Lincheck error: unsound matrix: {}", err)
            }
        }
    }
}

#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub enum RowcheckVerifierError {
    /// Error propagation
    DeserializationErr(winter_utils::DeserializationError),
    /// Error propagation
    MerkleTreeErr(winter_crypto::MerkleTreeError),
    /// Error propagation
    SmallPolyAdjustmentErr(),
    /// Error propagation
    FriVerifierErr(winter_fri::VerifierError),
}

impl From<winter_utils::DeserializationError> for RowcheckVerifierError {
    fn from(error: winter_utils::DeserializationError) -> Self {
        Self::DeserializationErr(error)
    }
}
impl From<winter_crypto::MerkleTreeError> for RowcheckVerifierError {
    fn from(error: winter_crypto::MerkleTreeError) -> Self {
        Self::MerkleTreeErr(error)
    }
}

impl From<winter_fri::VerifierError> for RowcheckVerifierError {
    fn from(error: winter_fri::VerifierError) -> Self {
        Self::FriVerifierErr(error)
    }
}

impl std::fmt::Display for RowcheckVerifierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            RowcheckVerifierError::DeserializationErr(err) => {
                writeln!(f, "Rowcheck deserialization error: {}", err)
            }
            RowcheckVerifierError::MerkleTreeErr(err) => {
                writeln!(f, "Rowcheck Merkle Tree error: {}", err)
            }
            RowcheckVerifierError::SmallPolyAdjustmentErr() => {
                writeln!(f, "Rowcheck Small Poly Adjustment error")
            }
            RowcheckVerifierError::FriVerifierErr(err) => {
                writeln!(f, "Rowcheck Fri error: {}", err)
            }
        }
    }
}

#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub enum FractalVerifierError {
    /// Error propagation
    LincheckVerifierErr(LincheckVerifierError),
    /// Error propagation
    RowcheckVerifierErr(RowcheckVerifierError),
}

impl From<LincheckVerifierError> for FractalVerifierError {
    fn from(error: LincheckVerifierError) -> Self {
        Self::LincheckVerifierErr(error)
    }
}

impl From<RowcheckVerifierError> for FractalVerifierError {
    fn from(error: RowcheckVerifierError) -> Self {
        Self::RowcheckVerifierErr(error)
    }
}

impl std::fmt::Display for FractalVerifierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            FractalVerifierError::LincheckVerifierErr(err) => {
                writeln!(f, "Lincheck error: {}", err)
            }
            FractalVerifierError::RowcheckVerifierErr(err) => {
                writeln!(f, "Rowcheck error: {}", err)
            }
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum LowDegreeVerifierError {
    /// Error propagation
    FriVerifierErr(VerifierError),
    /// Error propagation
    DeserializationErr(DeserializationError),
    PaddingErr,
}

impl From<VerifierError> for LowDegreeVerifierError {
    fn from(error: VerifierError) -> Self {
        Self::FriVerifierErr(error)
    }
}

impl From<DeserializationError> for LowDegreeVerifierError {
    fn from(error: DeserializationError) -> Self {
        Self::DeserializationErr(error)
    }
}

impl std::fmt::Display for LowDegreeVerifierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            LowDegreeVerifierError::FriVerifierErr(err) => {
                writeln!(f, "FRI Verifier Error: {}", err)
            }
            LowDegreeVerifierError::DeserializationErr(err) => {
                writeln!(f, "Winterfell Utils Deserialization Error: {}", err)
            }
            LowDegreeVerifierError::PaddingErr => {
                writeln!(f, "Complimentary Polynomial Check Failed")
            }
        }
    }
}