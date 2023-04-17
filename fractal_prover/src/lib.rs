use std::thread::AccessError;

use accumulator::Accumulator;
use channel::DefaultFractalProverChannel;
use errors::{AccumulatorProverError, ProverError};
use fractal_proofs::{FieldElement, LayeredProof, LowDegreeBatchProof, IopData, TopLevelProof};
use log;
use winter_crypto::ElementHasher;
use winter_fri::{FriOptions, ProverChannel};
use winter_math::StarkField;
pub mod accumulator;
pub mod channel;
pub mod errors;
pub mod lincheck_prover;
pub mod low_degree_batch_prover;
pub mod low_degree_prover;
pub mod prover;
pub mod rowcheck_prover;
pub mod sumcheck_prover;
#[cfg(test)]
mod tests;

pub const FRACTAL_LAYERS: usize = 3;

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

pub trait LayeredSubProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
>
{
    //pub type LayeredProof;
    //pub type LayeredProverState;
    //pub type DeferredAccumulator; //multipoly
    //pub fn run_next_layer(self, query: E, accumulator: DeferredAccumulator, proof: &mut LayeredProof, state: &mut LayeredProverState) -> Result<(), E>;
    // just one proof at the end, no need for multiple
    // You already need to keep track of each prover state, might as well have provers be stateful
    fn run_next_layer(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
    ) -> Result<(), ProverError>;

    /// Gets the id of the current layer, count starts at zero
    fn get_current_layer(&self) -> usize;

    /// Gets the total number of layers for this layered prover
    fn get_num_layers(&self) -> usize;

    /// Gets options for fractal proofs
    fn get_fractal_options(&self) -> FractalOptions<B>;

    // fn generate_proof(
    //     &mut self,
    //     public_input_bytes: Vec<u8>,
    // ) -> Result<LowDegreeBatchProof<B, E, H>, ProverError> {
    //     let options = self.get_fractal_options();
    //     let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
    //         options.evaluation_domain.len(),
    //         options.num_queries,
    //         public_input_bytes,
    //     );
    //     let mut acc = Accumulator::<B, E, H>::new(
    //         options.evaluation_domain.len(),
    //         B::ONE,
    //         options.evaluation_domain.clone(),
    //         options.num_queries,
    //         options.fri_options.clone(),
    //     );
    //     for i in 0..self.get_num_layers() {
    //         let query = channel.draw_fri_alpha();
    //         self.run_next_layer(query, &mut acc);
    //         acc.commit_layer(); //todo: do something with this
    //     }
    //     Ok(acc.create_fri_proof()?)
    // }
}

pub trait LayeredProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    D: IopData<B, E>,
>: LayeredSubProver<B, E, H>
{
    fn generate_proof(&mut self, public_input_bytes: Vec<u8>) -> Result<TopLevelProof<B,E,H>, ProverError>;
    // {
    //     let options = self.get_fractal_options();
    //     let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
    //         options.evaluation_domain.len(),
    //         options.num_queries,
    //         public_input_bytes,
    //     );
    //     let mut acc = Accumulator::<B, E, H>::new(
    //         options.evaluation_domain.len(),
    //         B::ONE,
    //         options.evaluation_domain.clone(),
    //         options.num_queries,
    //         options.fri_options.clone(),
    //     );
    //     for i in 0..self.get_num_layers() {
    //         let query = channel.draw_fri_alpha();
    //         self.run_next_layer(query, &mut acc);
    //         acc.commit_layer(); //todo: do something with this
    //     }
    //     Ok(acc.create_fri_proof()?)
    // }
}

/*
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
