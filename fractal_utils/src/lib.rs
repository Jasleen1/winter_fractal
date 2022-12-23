pub mod errors;
pub mod matrix_utils;
pub mod polynomial_utils;

#[cfg(test)]
mod tests;
pub type SmallFieldElement17 = fractal_math::smallprimefield::BaseElement<17, 3, 4>;
pub type SmallFieldElement13 = fractal_math::smallprimefield::BaseElement<13, 2, 2>;

pub static BLOWUP_FACTOR: usize = 8;
pub static FOLDING_FACTOR: usize = 4;
