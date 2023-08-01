// Copyright (c) Jasleen Malvai, Tom Jurek, Don Beaver.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

// This is the canonical orchestrator for running R1CS proof systems.
// It orchestrates proof generation and verification.
// It can run the prover and the verifier independently
// (useful for benchmarking, testing, developing).
// It supports a choice among prover implementations
// (provers who do provide the honest-prover proof artifact
// but may have optimizations such as batching, or maybe
// eventually dishonest provers).

use std::cmp::max;
use std::time::Instant;

use fractal_proofs::{fft, FractalProverOptions, TopLevelProof};
use fractal_prover::LayeredProver;
use fractal_prover::{prover::FractalProver, LayeredSubProver};
use fractal_prover::batched_lincheck_full_prover::BatchedFractalProver;
use fractal_utils::FractalOptions;
use fractal_verifier::verifier::verify_layered_fractal_proof_from_top as plain_verify_fractal_top;
use fractal_verifier::verifier_with_batched_lincheck::verify_layered_fractal_proof_from_top as batched_verify_fractal_top;
use models::r1cs::Matrix;
use winter_fri::FriOptions;

use structopt::StructOpt;

use fractal_indexer::{
    index::{build_index_domains, Index, IndexParams},
    indexed_matrix::index_matrix,
    snark_keys::*,
};

use models::jsnark_arith_parser::JsnarkArithReaderParser;
use models::jsnark_wire_parser::JsnarkWireReaderParser;
use reports::reporter::generate_flame_report;

use winter_crypto::hashers::{Blake3_256, Rp64_256};
use winter_crypto::ElementHasher;

use winter_math::fields::f64::BaseElement;
use winter_math::FieldElement;
use winter_math::StarkField;

#[cfg(feature = "flame_it")]
extern crate flame;

#[cfg(feature = "flame_it")]
extern crate flamer;

// At some point, someone wanted to try this.
//orchestrate_r1cs_example::<BaseElement, QuadExtension<BaseElement>, Rp64_256, 1>(

// Verbose-controlled output.
macro_rules! println_if {
    ($verbose:expr, $($x:tt)*) => { if $verbose { println!($($x)*) } }
}

#[cfg_attr(feature = "flame_it", flamer::flame("main"))]
fn main() {
    let options = OrchestratorOptions::from_args();
    println_if!(options.verbose,
        "Arith file {}, wire value file {}",
        options.arith_file, options.wires_file
    );

    let orchestrator = ProofSystemOrchestrator::<BaseElement, BaseElement, Blake3_256<BaseElement>, 1>::new(
        options.arith_file.clone(), options.wires_file.clone(), options.batched, options.verbose
    );
    orchestrator.orchestrate();

    #[cfg(feature = "flame_it")]
    {
        let filename = std::path::Path::new(&options.arith_file).file_stem().unwrap().to_str().unwrap();
        generate_flame_report(None, format!("r1cs:{filename}").as_str());
    }
}

// Orchestrates a proof system (P, V).
//
// TODO: pluggable Prover with a trait, rather than manual choice.
//
// Implementation note: there are two (P,V) implementations, selected by 'batched'.
// Each has its own (P,V) and the proof artifacts are not interchangable.
// There are three places where the code bifurcates: see import of
// plain_verify_fractal_top; see get_max_degree_constraint_batched;
// see Proving trait.

struct ProofSystemOrchestrator<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
 >{
    arith_file: String,
    wire_file: String,
    batched: bool,  // use plain or batched system, P and V both affected
    verbose: bool,

    _phantom_b: std::marker::PhantomData<B>,
    _phantom_e: std::marker::PhantomData<E>,
    _phantom_h: std::marker::PhantomData<H>,
}

impl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize> ProofSystemOrchestrator::<B, E, H, N> {

    fn new(arith_file: String, wire_file: String, batched: bool, verbose: bool) -> Self {
        Self {
            arith_file,
            wire_file,
            batched,
            verbose,
            _phantom_b: core::marker::PhantomData,
            _phantom_e: core::marker::PhantomData,
            _phantom_h: core::marker::PhantomData
        }
    }

    #[cfg_attr(feature = "flame_it", flamer::flame)]
    pub fn orchestrate(&self) {
        let (prover_key, verifier_key, fractal_options, wires, prover_options) = self.prepare();

        //orchestrate_r1cs_example_run::<B,E,H,N>(prover_key, verifier_key, fractal_options, wires, prover_options);
        self.prove_and_verify(prover_key, verifier_key, &fractal_options, &wires, prover_options);
    }

    pub fn prepare(&self) -> (ProverKey<B,E,H>, VerifierKey<B,H>, FractalOptions<B>, Vec<B>, FractalProverOptions<B>) {

        let now = Instant::now();

        let mut arith_parser = JsnarkArithReaderParser::<B>::new().unwrap();
        arith_parser.parse_arith_file(&self.arith_file, self.verbose);
        let mut r1cs = arith_parser.r1cs_instance;
        println_if!(self.verbose,
            "---------------------\nArith File parsed in {} ms",
            now.elapsed().as_millis()
        );

        let now = Instant::now();
        let mut wires_parser = JsnarkWireReaderParser::<B>::new().unwrap();
        wires_parser.parse_wire_file(&self.wire_file, self.verbose);
        let wires = wires_parser.wires;
        println_if!(self.verbose, "wire count = {}", wires.len());
        println_if !(self.verbose,
            "---------------------\nWire File parsed in {} ms",
            now.elapsed().as_millis()
        );

        // 0. Compute num_non_zero by counting max(number of non-zero elts across A, B, C).

        // let num_input_variables = r1cs.clone().num_cols();
        // let num_constraints = r1cs.clone().num_rows();
        // let num_non_zero = max(max(r1cs.A.l0_norm(), r1cs.B.l0_norm()), r1cs.C.l0_norm());
        // 1. Index this R1CS
        let num_input_variables = r1cs.num_cols().next_power_of_two();
        let num_non_zero = r1cs.max_num_nonzero().next_power_of_two();
        let num_constraints =
            max(max(r1cs.A.num_rows(), r1cs.B.num_rows()), r1cs.C.num_rows()).next_power_of_two();

        // Dependent on strategy (plain vs batched)
        let max_degree = match self.batched {
            false => FractalProver::<B, E, H>::get_max_degree_constraint(
                num_input_variables, num_non_zero, num_constraints),
            true => FractalProver::<B, E, H>::get_max_degree_constraint_batched(
                num_input_variables, num_non_zero, num_constraints),
        };
        // TODO: make the calculation of eta automated
        let eta = B::GENERATOR.exp(B::PositiveInteger::from(2 * B::TWO_ADICITY));
        let eta_k = B::GENERATOR.exp(B::PositiveInteger::from(1337 * B::TWO_ADICITY));
        // if num_non_zero <= num_vars {
        //     num_non_zero = num_non_zero * 2;
        // }
        let index_params = IndexParams::<B> {
            num_input_variables,
            num_constraints,
            num_non_zero,
            max_degree,
            eta,
            eta_k,
        };

        let degree_fs = r1cs.num_cols();

        let index_domains = build_index_domains::<B>(index_params.clone());
        println_if!(self.verbose, "built index domains");
        let indexed_a = index_matrix::<B>(&mut r1cs.A, &index_domains);
        r1cs.A = Matrix::new("dummy A", Vec::<Vec<B>>::new()).unwrap();
        println_if!(self.verbose, "indexed matrix a");
        let indexed_b = index_matrix::<B>(&mut r1cs.B, &index_domains);
        r1cs.B = Matrix::new("dummy B", Vec::<Vec<B>>::new()).unwrap();
        println_if!(self.verbose, "indexed matrix b");
        let indexed_c = index_matrix::<B>(&mut r1cs.C, &index_domains);
        r1cs.C = Matrix::new("dummy C", Vec::<Vec<B>>::new()).unwrap();
        println_if!(self.verbose, "indexed matrices");

        // This is the index i.e. the pre-processed data for this r1cs
        let preproc_index = Index::new(index_params.clone(), indexed_a, indexed_b, indexed_c);

        // TODO: the IndexDomains should already guarantee powers of two, so why add extraneous bit or use next_power_of_two?

        let size_subgroup_h = index_domains.h_field.len().next_power_of_two();
        let size_subgroup_k = index_domains.k_field.len().next_power_of_two();

        let evaluation_domain =
            winter_math::get_power_series(index_domains.l_field_base, index_domains.l_field_len);

        let summing_domain = index_domains.k_field;

        let h_domain = index_domains.h_field;
        let lde_blowup = 4;
        let num_queries = 16;
        let fri_options = FriOptions::new(lde_blowup, 4, 32);
        //println!("h_domain: {:?}, summing_domain: {:?}, evaluation_domain: {:?}", &h_domain, &summing_domain, &evaluation_domain);
        let fractal_options: FractalOptions<B> = FractalOptions::<B> {
            degree_fs,
            size_subgroup_h,
            size_subgroup_k,
            summing_domain: summing_domain.clone(),
            evaluation_domain: evaluation_domain.clone(),
            h_domain: h_domain.clone(),
            eta,
            eta_k,
            fri_options: fri_options.clone(),
            num_queries,
        };

        let h_domain_twiddles = fft::get_twiddles(size_subgroup_h);
        let h_domain_inv_twiddles = fft::get_inv_twiddles(size_subgroup_h);
        let k_domain_twiddles = fft::get_twiddles(size_subgroup_k);
        let k_domain_inv_twiddles = fft::get_inv_twiddles(size_subgroup_k);
        let l_domain_twiddles = fft::get_twiddles(evaluation_domain.len());
        let l_domain_inv_twiddles = fft::get_inv_twiddles(evaluation_domain.len());
        let prover_options: FractalProverOptions<B> = FractalProverOptions::<B> {
            degree_fs,
            size_subgroup_h,
            size_subgroup_k,
            summing_domain,
            evaluation_domain,
            h_domain,
            h_domain_twiddles,
            h_domain_inv_twiddles,
            k_domain_twiddles,
            k_domain_inv_twiddles,
            l_domain_twiddles,
            l_domain_inv_twiddles,
            eta,
            eta_k,
            fri_options: fri_options.clone(),
            num_queries,
        };

        let (prover_key, verifier_key) =
            generate_prover_and_verifier_keys::<B, E, H>(preproc_index, &fractal_options).unwrap();
        println_if!(self.verbose, "Prover and verifier keys generated");

        // Return:
        (prover_key, verifier_key, fractal_options, wires, prover_options)
    }

    pub fn prove(
        &self,
        pub_inputs_bytes: &Vec<u8>,
        prover_key: ProverKey<B,E,H>,
        wires: &Vec<B>,
        prover_options: FractalProverOptions<B>,
    ) -> TopLevelProof<B, E, H> {

        match self.batched {
            false => {
                let claimant: &dyn Proving<B, E, H, N> = &FractalProverPlainImpl::new();
                claimant.issue_proof(pub_inputs_bytes, prover_key, wires, prover_options)
            }

            true => {
                let claimant: &dyn Proving<B, E, H, N> = &FractalProverBatchedImpl::new();
                claimant.issue_proof(pub_inputs_bytes, prover_key, wires, prover_options)
            }
        }
    }

    pub fn verify(
        &self,
        proof: TopLevelProof<B, E, H>,
        pub_inputs_bytes: &Vec<u8>,
        verifier_key: &VerifierKey<B,H>,
        fractal_options: &FractalOptions<B>,
    ) {

        let now = Instant::now();
        // Choose from different implementations of
        // verify_layered_fractal_proof_from_top(&verifier_key, &proof, &pub_inputs_bytes, &fractal_options).unwrap(),
        match self.batched {
            false => plain_verify_fractal_top(&verifier_key, &proof, &pub_inputs_bytes, &fractal_options).unwrap(),
            true => batched_verify_fractal_top(&verifier_key, &proof, &pub_inputs_bytes, &fractal_options).unwrap(),
        }

        println!(
            "---------------------\nProof verified in {} ms",
            now.elapsed().as_millis()
        );
        println!("Finished!");
    }

    pub fn prove_and_verify(
        &self,
        prover_key: ProverKey<B,E,H>,
        verifier_key: VerifierKey<B,H>,
        fractal_options: &FractalOptions<B>,
        wires: &Vec<B>,
        prover_options: FractalProverOptions<B>,
    ) {
        // TODO: rationalize this, it seems a vestige of fibonacci examples.
        let pub_inputs_bytes = vec![0u8, 1u8, 2u8];
        //let pub_inputs_bytes = vec![];

        let proof = self.prove(&pub_inputs_bytes, prover_key, wires, prover_options);
        self.verify(proof, &pub_inputs_bytes, &verifier_key, &fractal_options);
    }
}

// TODO: push this trait closer to provers and verifiers instead of using it here
// as a hacky wrapper.
pub trait Proving<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize> {

    fn issue_proof(&self, 
        pub_inputs_bytes: &Vec<u8>,
        prover_key: ProverKey<B,E,H>,
        wires: &Vec<B>,
        prover_options: FractalProverOptions<B>,
    ) -> TopLevelProof<B, E, H>;
}

// Fractal prover, simplest implementation.  The only meaningful line is "let prover = FractalProver".
struct FractalProverPlainImpl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>> {
    _phantom_b: std::marker::PhantomData<B>,
    _phantom_e: std::marker::PhantomData<E>,
    _phantom_h: std::marker::PhantomData<H>,
}
impl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>> FractalProverPlainImpl<B, E, H> {
    
    fn new() -> Self {
        Self {
            _phantom_b: core::marker::PhantomData,
            _phantom_e: core::marker::PhantomData,
            _phantom_h: core::marker::PhantomData
        }
    }
}

impl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize> Proving<B, E, H, N> for FractalProverPlainImpl<B, E, H> {

    fn issue_proof(
        &self,
        pub_inputs_bytes: &Vec<u8>,
        prover_key: ProverKey<B,E,H>,
        wires: &Vec<B>,
        prover_options: FractalProverOptions<B>,
    ) -> TopLevelProof<B, E, H> {
        let mut prover =
            FractalProver::<B, E, H>::new(prover_key.into(), vec![], wires.clone(), pub_inputs_bytes.clone());
        let now = Instant::now();
        let proof = prover
            .generate_proof(&None, pub_inputs_bytes.clone(), &prover_options)
            .unwrap();
        println!(
            "---------------------\nProof generated (fractal) in {} ms",
            now.elapsed().as_millis()
        );
        proof
    }
}

// Fractal prover whose implementation uses batching.  The only meaningful line is "let prover = BatchedFractalProver".
struct FractalProverBatchedImpl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>> {
    _phantom_b: std::marker::PhantomData<B>,
    _phantom_e: std::marker::PhantomData<E>,
    _phantom_h: std::marker::PhantomData<H>,
}
impl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>> FractalProverBatchedImpl<B, E, H> {

    fn new() -> Self {
        Self {
            _phantom_b: core::marker::PhantomData,
            _phantom_e: core::marker::PhantomData,
            _phantom_h: core::marker::PhantomData
        }
    }
}

impl<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize> Proving<B, E, H, N> for FractalProverBatchedImpl<B, E, H> {

    fn issue_proof(
        &self,
        pub_inputs_bytes: &Vec<u8>,
        prover_key: ProverKey<B,E,H>,
        wires: &Vec<B>,
        prover_options: FractalProverOptions<B>,
    ) -> TopLevelProof<B, E, H> {
        let mut prover =
            BatchedFractalProver::<B, E, H>::new(prover_key.into(), vec![], wires.clone(), pub_inputs_bytes.clone());
        let now = Instant::now();
        let proof = prover
            .generate_proof(&None, pub_inputs_bytes.clone(), &prover_options)
            .unwrap();
        println!(
            "---------------------\nProof generated (batched fractal) in {} ms",
            now.elapsed().as_millis()
        );
        proof
    }
}


#[derive(StructOpt, Debug)]
#[structopt(name = "jsnark-parser", about = "Jsnark file parsing")]
struct OrchestratorOptions {
    /// Jsnark .arith file to parse.
    #[structopt(
        short = "a",
        long = "arith_file",
        //default_value = "fractal_examples/jsnark_outputs/fibonacciexample_17.arith"
        default_value = "fractal_examples/jsnark_outputs/fftexample_5.arith"
    )]
    arith_file: String,

    /// Jsnark .in or .wires file to parse.
    #[structopt(
        short = "w",
        long = "wire_file",
        //default_value = "fractal_examples/jsnark_outputs/fibonacciexample_17.wires"
        default_value = "fractal_examples/jsnark_outputs/fftexample_5.wires"
    )]
    wires_file: String,

    /// Elect (poly)batching implementation of (P,V)
    #[structopt(short = "b", long = "batched")]
    batched: bool,

    /// Verbose logging and reporting.
    #[structopt(short = "v", long = "verbose")]
    verbose: bool,
}
