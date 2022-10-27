use fractal_math::StarkField;
use winter_fri::{FriOptions, ProverChannel};
use winter_utils::{Serializable, ByteWriter};

pub mod errors;
pub mod matrix_utils;
pub mod polynomial_utils;

#[cfg(test)]
mod tests;
pub type SmallFieldElement17 = fractal_math::smallprimefield::BaseElement<17, 3, 4>;
pub type SmallFieldElement13 = fractal_math::smallprimefield::BaseElement<13, 2, 2>;

pub static BLOWUP_FACTOR: usize = 8;
pub static FOLDING_FACTOR: usize = 4;

#[derive(Clone)]
pub struct FractalOptions<B: StarkField> {
    pub degree_fs: usize,
    pub size_subgroup_h: usize,
    pub size_subgroup_k: usize,
    // K domain in paper
    pub summing_domain: Vec<B>,
    // L domain in paper
    pub evaluation_domain: Vec<B>,
    // H domain in paper
    pub h_domain: Vec<B>,
    pub eta: B,
    pub eta_k: B,
    pub fri_options: FriOptions,
    pub num_queries: usize,
}

impl<B: StarkField> Serializable for FractalOptions<B> {
    /// Serializes `self` and writes the resulting bytes into the `target`.
    fn write_into<W: ByteWriter>(&self, target: &mut W) {
        target.write_u16(self.degree_fs as u16); 
        target.write_u16(self.size_subgroup_h as u16);
        target.write_u16(self.size_subgroup_k as u16);
        target.write_u16(self.evaluation_domain.len() as u16);
        target.write_u16(self.num_queries as u16);
        self.eta.write_into(target);
        self.eta_k.write_into(target);
    }
}

// impl Deserializable for FractalOptions<E> {
//     /// Reads proof context from the specified `source` and returns the result.
//     ///
//     /// # Errors
//     /// Returns an error of a valid Context struct could not be read from the specified `source`.
//     fn read_from<R: ByteReader>(source: &mut R) -> Result<Self, DeserializationError> {
//         // read and validate trace layout info
//         let trace_layout = TraceLayout::read_from(source)?;

//         // read and validate trace length (which was stored as a power of two)
//         let degree_fs = source.read_u16()? as usize;
//         let size_subgroup_h = source.read_u16()? as usize;
//         let size_subgroup_k = source.read_u16()? as usize;
//         let eval_domain_len = source.read_u16()? as usize;
//         let evaluation_domain = compute

//         // read trace metadata
//         let num_meta_bytes = source.read_u16()? as usize;
//         let trace_meta = if num_meta_bytes != 0 {
//             source.read_u8_vec(num_meta_bytes)?
//         } else {
//             vec![]
//         };

//         // read and validate field modulus bytes
//         let num_modulus_bytes = source.read_u8()? as usize;
//         if num_modulus_bytes == 0 {
//             return Err(DeserializationError::InvalidValue(
//                 "field modulus cannot be an empty value".to_string(),
//             ));
//         }
//         let field_modulus_bytes = source.read_u8_vec(num_modulus_bytes)?;

//         // read options
//         let options = ProofOptions::read_from(source)?;

//         Ok(Context {
//             trace_layout,
//             trace_length,
//             trace_meta,
//             field_modulus_bytes,
//             options,
//         })
//     }
// }