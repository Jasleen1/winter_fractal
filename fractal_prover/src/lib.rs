use log;
use winter_fri::FriOptions;
use winter_math::StarkField;
mod errors;
mod lincheck_prover;
pub mod low_degree_batch_prover;
pub mod low_degree_prover;
pub mod prover;
mod rowcheck_prover;
pub mod sumcheck_prover;
#[cfg(test)]
mod tests;

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
