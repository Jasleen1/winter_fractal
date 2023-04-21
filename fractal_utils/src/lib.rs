use fractal_math::StarkField;
use winter_fri::FriOptions;

pub mod channel;
pub mod errors;
pub mod matrix_utils;
pub mod polynomial_utils;

#[cfg(test)]
mod tests;
pub type SmallFieldElement17 = fractal_math::smallprimefield::BaseElement<17, 3, 4>;
pub type SmallFieldElement13 = fractal_math::smallprimefield::BaseElement<13, 2, 2>;

pub static BLOWUP_FACTOR: usize = 4;
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


#[derive(Clone)]
pub struct FractalProverOptions<B: StarkField> {
    pub degree_fs: usize,
    pub size_subgroup_h: usize,
    pub size_subgroup_k: usize,
    // K domain in paper
    pub summing_domain: Vec<B>,
    // L domain in paper
    pub evaluation_domain: Vec<B>,
    // H domain in paper
    pub h_domain: Vec<B>,
    pub h_domain_twiddles: Vec<B>,
    pub h_domain_inv_twiddles: Vec<B>,
    pub k_domain_twiddles: Vec<B>,
    pub k_domain_inv_twiddles: Vec<B>,
    pub l_domain_twiddles: Vec<B>,
    pub l_domain_inv_twiddles: Vec<B>,
    pub eta: B,
    pub eta_k: B,
    pub fri_options: FriOptions,
    pub num_queries: usize,
}

#[derive(Clone)]
pub struct FractalVerifierOptions<B: StarkField> {
    pub degree_fs: usize,
    // H domain in paper
    pub size_subgroup_h: usize,
    // K domain in paper
    pub size_subgroup_k: usize,
    // L domain in paper
    pub size_subgroup_l: usize,
    pub eta: B,
    pub eta_k: B,
    pub fri_options: FriOptions,
    pub num_queries: usize,
}