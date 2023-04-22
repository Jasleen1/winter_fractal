use std::thread::AccessError;

use errors::ProverError;
use fractal_accumulator::{accumulator::Accumulator, errors::AccumulatorProverError};
use fractal_indexer::{snark_keys::ProverKey, index::IndexParams};
use fractal_proofs::{
    FieldElement, FractalProverOptions, IopData, LayeredProof, LowDegreeBatchProof, TopLevelProof,
};
use fractal_utils::channel::DefaultFractalProverChannel;
use log;
use winter_crypto::ElementHasher;
use winter_fri::{FriOptions, ProverChannel};
use winter_math::StarkField;
pub mod errors;
pub mod lincheck_prover;
pub mod prover;
pub mod rowcheck_prover;
pub mod sumcheck_prover;

#[cfg(feature = "flame_it")]
extern crate flame;
#[cfg(feature = "flame_it")]
#[macro_use]
extern crate flamer;

#[cfg(test)]
mod tests;

/// This constant stores the number of layers in the Fractal IOP
pub const FRACTAL_LAYERS: usize = 3;

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

/// This trait is meant to be _like_ a layered IOP prover but we don't want it to do any commitments.
/// This is why we called it the LayeredSubProver, since we will be implementing it in subroutines of an actual IOP
/// prover to maintain a semblance of modularity.
/// This trait includes subroutines associated with a layered IOP.
pub trait LayeredSubProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
>
{
    /// Run the next layer of this IOP prover
    fn run_next_layer(
        &mut self,
        query: E,
        accumulator: &mut Accumulator<B, E, H>,
        options: &FractalProverOptions<B>,
    ) -> Result<(), ProverError>;

    /// Gets the id of the current layer, count starts at zero
    fn get_current_layer(&self) -> usize;

    /// Gets the total number of layers for this layered prover
    fn get_num_layers(&self) -> usize;

    fn get_max_degree_constraint(num_input_variables: usize, num_non_zero: usize, num_constraints: usize) -> usize;
}

/// This is a trait for a layered IOP prover which also implements the trait
/// [`LayeredSubProver`]. The main additional function is the actual proof generation,
/// which takes place in the [`LayeredProver::generate_proof`] function and returns a
/// proof of type [`TopLevelProof`].
pub trait LayeredProver<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    D: IopData<B, E>,
>: LayeredSubProver<B, E, H>
{
    /// Generate proof for a [`LayeredProver`]
    fn generate_proof(
        &mut self,
        prover_key: &Option<ProverKey<B, E, H>>,
        public_input_bytes: Vec<u8>,
        options: &FractalProverOptions<B>,
    ) -> Result<TopLevelProof<B, E, H>, ProverError>;
    // BELOW IS A SAMPLE IMPLEMENTATION OF THIS FUNCTION
    // This function, however, needs a special-purpose implementation,
    // depending on the specific IOP.
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
