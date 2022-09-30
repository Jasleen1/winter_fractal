// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under both the MIT license found in the
// LICENSE-MIT file in the root directory of this source tree and the Apache
// License, Version 2.0 found in the LICENSE-APACHE file in the root directory
// of this source tree.

//! Errors for various data structure operations.
//use winter_fri::VerifierError;
use winter_utils::DeserializationError;
use low_degree::errors::LowDegreeVerifierError;

#[derive(Debug, PartialEq)]
pub enum SumcheckVerifierError {
    /// Error propagation
    FriVerifierErr(LowDegreeVerifierError),
    /// Error propagation
    DeserializationErr(DeserializationError),
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
        }
    }
}