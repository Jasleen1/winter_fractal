use log;
use winter_fri::FriOptions;
use winter_math::StarkField;
pub mod channel;
mod errors;
mod lincheck_prover;
pub mod low_degree_batch_prover;
pub mod low_degree_prover;
pub mod prover;
mod rowcheck_prover;
pub mod sumcheck_prover;
pub mod accumulator;
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

//multiple proofs can be run in parallel starting from the same transcript state
//Lincheck is 2 layers, Rowcheck is 1. Before all this, need to commit to f_az, f_bz, etc

//take some commitable information, add it to a vec of such things (accumulator)
//fn add_to_next_layer (&mut self, channel, mut someVec encoding/IOP data?) -> Result<()>

//mutate the channel with previously-committed information
//fn finalize_layer(&mut self, &mut channel) -> Result<()>

//instead of blackbox lincheck, have lincheck layer 1, layer 2, run FRI
// get passed a type which implements layeredproof

//do beta sampling outside of lincheck?
//then pass the same beta into all three linchecks
//t_alpha(x)
//instead of creating three different merkle trees for committing to t_alpha, put all the evaluations in one tree
//i'th leaf of the tree contains i'th evaulation for each t_alpha_evals in each of the three linchecks
//check the paths outside of the lincheck
//additionally commit to the two polynomials in rationalSumcheck

//non-interaction IOP protocol trait -> include associated type: num IOP layers
//do a for loop: 0..num_layers, each loop you do one sample and make a new synced channel
//layeredproof divorrced from FRI

//Layerd proof object-thing finalizes at each layer to sync the channels. At the end run the BatchFRIProver
//-you probably need to pass it around to get all the polynomials and degrees you need for the end part
//maybe pass a cloned channel into each then smoosh them together

/*
pub trait LayeredProof{
    fn add_proof_layer(&mut self, channel) -> 
}

pub trait LayeredProver<
B: StarkField,
E: FieldElement<BaseField = B>
>{
    pub type LayeredProof;
    pub type LayeredProverState;
    pub type DeferredAccumulator; //multipoly
    pub fn run_next_layer(self, query: E, accumulator: DeferredAccumulator, proof: &mut LayeredProof, state: &mut LayeredProverState) -> Result<(), E>;

}

let lincheck_prover_a = LincheckProver::<B, E, H>::new(
    //alpha,
    &matrix_index,
    prod_m_z_coeffs.to_vec(),
    z_coeffs.to_vec(),
    &self.options,
);
// ditto for b and c

let mut commitment_vec: Vec<H> = vec![];
let mut deferred_batch_low_degree_prover = LowDegreeBatchProver::new();
// make a data type to handle both low_degree_batch_prover and the mult_poly, 
// pass this in mutably into run_next_layer instead
//multipoly commits to eval_domain evals
for i in 0..lincheck.layers.len(){
    let layer_i_multi_eval = MultiEval::new();
    let query = channel.draw_value_depending_on_lincheck_layer(i);
    // could input query into the ProverState
    // proving function decides if the query is valid. If not, use it to seed another channel and draw until good
    let (proof_state_a: LincheckProver::ProverState, proof_option: Option<LincheckProof>) = lincheck_prover_a.run_next_layer(query, layer_i_multi_eval, deferred_batch_low_degree_prover);
    // same for b and c
    // some kind of evaluate and commit function needed below VVV
    let next_com = layer_i_multi_eval.commit_polynomial_evaluations();
    commitment_vec.push(next_com);
    channel.reseed(next_com);
}


impl LayeredProver for LincheckProver{
    type LayeredProverState = LincheckProverState;
}

*/