// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use winter_crypto::Hasher;
use winter_utils::{
    collections::Vec, ByteReader, ByteWriter, Deserializable, DeserializationError, Serializable,
    SliceReader,
};


#[derive(Debug, Clone, Default, Eq, PartialEq)]
pub struct FractalCommitments(Vec<u8>);

impl FractalCommitments {

    pub fn new<H: Hasher>(
        lincheck_a_roots: Vec<H::Digest>,
        lincheck_b_roots: Vec<H::Digest>,
        lincheck_c_roots: Vec<H::Digest>,
        rowcheck_roots: Vec<H::Digest>,
    ) -> Self {
        let mut bytes = Vec::new();
        bytes.write(lincheck_a_roots);
        bytes.write(lincheck_b_roots);
        bytes.write(lincheck_c_roots);
        bytes.write(rowcheck_roots);
        // bytes.write(fri_roots);
        FractalCommitments(bytes)
    }

    pub fn add<H: Hasher>(&mut self, com: &H::Digest) {
        com.write_into(&mut self.0);
    }

    #[allow(clippy::type_complexity)]
    pub fn parse<H: Hasher>(
        self,
        num_fri_layers: usize,
    ) -> Result<(Vec<H::Digest>, Vec<H::Digest>, Vec<H::Digest>, Vec<H::Digest>), DeserializationError> {
        let mut reader = SliceReader::new(&self.0);

    }   


}