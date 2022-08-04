// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under both the MIT license found in the
// LICENSE-MIT file in the root directory of this source tree and the Apache
// License, Version 2.0 found in the LICENSE-APACHE file in the root directory
// of this source tree.

//! Errors for various data structure operations.
use winter_fri::VerifierError;
use winter_utils::DeserializationError;

/// Symbolizes a AkdError, thrown by the akd.
#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub enum SumcheckError {
    /// Error propagation
    FriVerifierErr(VerifierError),
    /// Error propagation
    DeserializationErr(DeserializationError),
}

impl From<VerifierError> for SumcheckError {
    fn from(error: VerifierError) -> Self {
        Self::FriVerifierErr(error)
    }
}

impl From<DeserializationError> for SumcheckError {
    fn from(error: DeserializationError) -> Self {
        Self::DeserializationErr(error)
    }
}

impl std::fmt::Display for SumcheckError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            SumcheckError::FriVerifierErr(err) => {
                writeln!(f, "FRI Verifier Error: {}", err)
            }
            SumcheckError::DeserializationErr(err) => {
                writeln!(f, "Winterfell Utils Deserialization Error: {}", err)
            }
        }
    }
}