use winter_fri::VerifierError;
use winter_utils::DeserializationError;

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