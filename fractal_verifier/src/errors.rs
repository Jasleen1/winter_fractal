// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under both the MIT license found in the
// LICENSE-MIT file in the root directory of this source tree and the Apache
// License, Version 2.0 found in the LICENSE-APACHE file in the root directory
// of this source tree.

//! Errors for various data structure operations.
use fractal_accumulator::errors::AccumulatorProverError;
use fractal_accumulator_verifier::errors::AccumulatorVerifierError;
use fractal_proofs::{errors::FractalUtilError, DeserializationError};
use fractal_prover::errors::{LincheckError, ProverError};
use low_degree_verifier::errors::LowDegreeVerifierError;
use thiserror::Error;
use winter_crypto::{MerkleTreeError, RandomCoinError};
use winter_fri::VerifierError;

#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub enum LincheckVerifierError {
    /// Error propagation
    UnsoundProduct(SumcheckVerifierError),
    /// Error propagation
    UnsoundMatrix(SumcheckVerifierError),
    /// Error propagation
    AccumulatorVerifierErr(AccumulatorVerifierError),
}

impl From<SumcheckVerifierError> for LincheckVerifierError {
    fn from(error: SumcheckVerifierError) -> Self {
        Self::UnsoundProduct(error)
    }
}

impl From<AccumulatorVerifierError> for LincheckVerifierError {
    fn from(error: AccumulatorVerifierError) -> Self {
        Self::AccumulatorVerifierErr(error)
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
            LincheckVerifierError::AccumulatorVerifierErr(err) => {
                writeln!(f, "Accumulator Verifer error: {}", err)
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
    /// Error propagation
    LowDegreeVerifierErr(LowDegreeVerifierError),
    /// This error is thrown if the computed value of the polynomial s in rowcheck
    /// does not match the value that is sent from the prover
    ComputedValueMismatchErr(String),
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

impl From<LowDegreeVerifierError> for RowcheckVerifierError {
    fn from(error: LowDegreeVerifierError) -> Self {
        Self::LowDegreeVerifierErr(error)
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
            RowcheckVerifierError::LowDegreeVerifierErr(err) => {
                writeln!(f, "Rowcheck Low Degree Verifier error: {}", err)
            }
            RowcheckVerifierError::ComputedValueMismatchErr(err) => {
                writeln!(f, "Rowcheck error in checking computed values: {}", err)
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
    /// Error propagation
    FractalUtilErr(FractalUtilError),
    /// Error propagation
    AccumulatorVerifierErr(AccumulatorVerifierError),
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

impl From<FractalUtilError> for FractalVerifierError {
    fn from(error: FractalUtilError) -> Self {
        Self::FractalUtilErr(error)
    }
}

impl From<AccumulatorVerifierError> for FractalVerifierError {
    fn from(error: AccumulatorVerifierError) -> Self {
        Self::AccumulatorVerifierErr(error)
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
            FractalVerifierError::FractalUtilErr(err) => {
                writeln!(f, "Fractal utils error: {}", err)
            }
            FractalVerifierError::AccumulatorVerifierErr(err) => {
                writeln!(f, "Accumulator Verifer error: {}", err)
            }
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum SumcheckVerifierError {
    /// Error propagation
    FriVerifierErr(LowDegreeVerifierError),
    /// Error propagation
    DeserializationErr(DeserializationError),
    /// The e polynomial does not match up with the g_polynomial as needed
    ConsistentValuesErr(usize),
}

impl From<LowDegreeVerifierError> for SumcheckVerifierError {
    fn from(error: LowDegreeVerifierError) -> Self {
        Self::FriVerifierErr(error)
    }
}

impl From<DeserializationError> for SumcheckVerifierError {
    fn from(error: DeserializationError) -> Self {
        Self::DeserializationErr(error)
    }
}

impl std::fmt::Display for SumcheckVerifierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            SumcheckVerifierError::FriVerifierErr(err) => {
                writeln!(f, "FRI Verifier Error: {}", err)
            }
            SumcheckVerifierError::DeserializationErr(err) => {
                writeln!(f, "Winterfell Utils Deserialization Error: {}", err)
            }
            SumcheckVerifierError::ConsistentValuesErr(err) => {
                writeln!(
                    f,
                    "Sumcheck verifier's equality test failed at position: {}",
                    err
                )
            }
        }
    }
}

#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub enum TestingError {
    ProverErr(ProverError),
    VerifierErr(FractalVerifierError),
    LincheckProverErr(LincheckError),
    RowcheckVerifierErr(RowcheckVerifierError),
    LincheckVerifierErr(LincheckVerifierError),
    AccumulatorProverErr(AccumulatorProverError),
    MerkleTreeErr(MerkleTreeError),
    AccumulatorVerifierErr(AccumulatorVerifierError),
}

impl From<FractalVerifierError> for TestingError {
    fn from(err: FractalVerifierError) -> Self {
        TestingError::VerifierErr(err)
    }
}

impl From<ProverError> for TestingError {
    fn from(err: ProverError) -> Self {
        TestingError::ProverErr(err)
    }
}

impl From<LincheckVerifierError> for TestingError {
    fn from(err: LincheckVerifierError) -> Self {
        TestingError::LincheckVerifierErr(err)
    }
}

impl From<RowcheckVerifierError> for TestingError {
    fn from(err: RowcheckVerifierError) -> Self {
        TestingError::RowcheckVerifierErr(err)
    }
}

impl From<LincheckError> for TestingError {
    fn from(err: LincheckError) -> Self {
        TestingError::LincheckProverErr(err)
    }
}

impl From<AccumulatorProverError> for TestingError {
    fn from(err: AccumulatorProverError) -> Self {
        TestingError::AccumulatorProverErr(err)
    }
}

impl From<AccumulatorVerifierError> for TestingError {
    fn from(err: AccumulatorVerifierError) -> Self {
        TestingError::AccumulatorVerifierErr(err)
    }
}

impl From<MerkleTreeError> for TestingError {
    fn from(err: MerkleTreeError) -> Self {
        TestingError::MerkleTreeErr(err)
    }
}

impl std::fmt::Display for TestingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            TestingError::ProverErr(err) => {
                writeln!(f, "Fractal prover error: {}", err)
            }
            TestingError::VerifierErr(err) => {
                writeln!(f, "Fractal verifier error: {}", err)
            }
            TestingError::RowcheckVerifierErr(err) => {
                writeln!(f, "Fractal rowcheck verifier error: {}", err)
            }
            TestingError::LincheckVerifierErr(err) => {
                writeln!(f, "Fractal lincheck verifier error: {}", err)
            }
            TestingError::AccumulatorProverErr(err) => {
                writeln!(f, "Fractal accumulator error: {}", err)
            }
            TestingError::AccumulatorVerifierErr(err) => {
                writeln!(f, "Fractal accumulator error: {}", err)
            }
            TestingError::MerkleTreeErr(err) => {
                writeln!(f, "Fractal accumulator error: {}", err)
            }
            TestingError::LincheckProverErr(err) => {
                writeln!(f, "Fractal lincheck prover error: {}", err)
            }
        }
    }
}
