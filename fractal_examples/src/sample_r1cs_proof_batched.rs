// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

#![allow(dead_code,unused_imports)]

use core::num;
use std::cmp::max;
use std::time::Instant;

use fractal_proofs::{fft, FractalProverOptions, Serializable};
use fractal_prover::batched_lincheck_full_prover::BatchedFractalProver;
use fractal_prover::LayeredProver;
use fractal_prover::{prover::FractalProver, LayeredSubProver};
use fractal_utils::FractalOptions;
use fractal_verifier::verifier_with_batched_lincheck::verify_layered_fractal_proof_from_top;
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
use winter_math::fields::QuadExtension;
use winter_math::FieldElement;
use winter_math::StarkField;

#[cfg(feature = "flame_it")]
extern crate flame;
#[cfg(feature = "flame_it")]
#[macro_use]
extern crate flamer;

#[cfg_attr(feature = "flame_it", flame("main"))]
fn main() {
    let mut options = ExampleOptions::from_args();
    options.verbose = false;
    if options.verbose {
        println!(
            "Arith file {}, wire value  file {}",
            options.arith_file, options.wires_file
        );
    }

    // call orchestrate_r1cs_example
    //orchestrate_r1cs_example::<BaseElement, QuadExtension<BaseElement>, Rp64_256, 1>(
    orchestrate_r1cs_example::<BaseElement, BaseElement, Blake3_256<BaseElement>, 1>(
        &options.arith_file,
        &options.wires_file,
        options.verbose,
    );

    #[cfg(feature = "flame_it")]
    {
        let filename = std::path::Path::new(&options.arith_file).file_stem().unwrap().to_str().unwrap();
        generate_flame_report(None, format!("r1cs_batched:{filename}").as_str());
    }
}

#[cfg_attr(feature = "flame_it", flame)]
pub(crate) fn orchestrate_r1cs_example<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
>(
    arith_file: &str,
    wire_file: &str,
    verbose: bool,
) {
    let now = Instant::now();
    let mut arith_parser = JsnarkArithReaderParser::<B>::new().unwrap();
    arith_parser.parse_arith_file(&arith_file, verbose);
    println!("Parsed the arith file");
    let mut r1cs = arith_parser.r1cs_instance;
    println!(
        "---------------------\nArith File parsed in {} ms",
        now.elapsed().as_millis()
    );

    let now = Instant::now();
    let mut wires_parser = JsnarkWireReaderParser::<B>::new().unwrap();
    println!("Got the wire parser");
    wires_parser.parse_wire_file(&wire_file, verbose);
    let wires = wires_parser.wires;
    println!("Len wires = {}", wires.len());
    println!(
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
    let max_degree = FractalProver::<B, E, H>::get_max_degree_constraint_batched(
        num_input_variables,
        num_non_zero,
        num_constraints,
    );
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
    println!("built index domains");
    let indexed_a = index_matrix::<B>(&mut r1cs.A, &index_domains);
    println!("ret again {}", 4);
    r1cs.A = Matrix::new("dummy A", Vec::<Vec<B>>::new()).unwrap();
    println!("indexed matrix a");
    let indexed_b = index_matrix::<B>(&mut r1cs.B, &index_domains);
    r1cs.B = Matrix::new("dummy B", Vec::<Vec<B>>::new()).unwrap();
    println!("indexed matrix b");
    let indexed_c = index_matrix::<B>(&mut r1cs.C, &index_domains);
    r1cs.C = Matrix::new("dummy C", Vec::<Vec<B>>::new()).unwrap();
    println!("indexed matrices");
    // This is the index i.e. the pre-processed data for this r1cs
    let index = Index::new(index_params.clone(), indexed_a, indexed_b, indexed_c);

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
    let options: FractalOptions<B> = FractalOptions::<B> {
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
        generate_prover_and_verifier_keys::<B, E, H>(index, &options).unwrap();
    println!("Prover and verifier keys generated");
    let pub_inputs_bytes = vec![0u8, 1u8, 2u8];
    //let pub_inputs_bytes = vec![];
    let mut prover = BatchedFractalProver::<B, E, H>::new(
        prover_key.into(),
        vec![],
        wires,
        pub_inputs_bytes.clone(),
    );
    let now = Instant::now();
    let proof = prover
        .generate_proof(&None, pub_inputs_bytes.clone(), &prover_options)
        .unwrap();
    println!(
        "---------------------\nProof generated in {} ms",
        now.elapsed().as_millis()
    );

    let now = Instant::now();
    verify_layered_fractal_proof_from_top(&verifier_key, &proof, &pub_inputs_bytes, &options)
        .unwrap();
    println!(
        "---------------------\nProof verified in {} ms",
        now.elapsed().as_millis()
    );
    println!("Success!");

    // println!(
    //     "Verified: {:?}",
    //     fractal_verifier::verifier::verify_fractal_proof::<B, E, H>(
    //         verifier_key,
    //         proof,
    //         pub_inputs_bytes,
    //         options
    //     )
    // );
}

#[derive(StructOpt, Debug)]
#[structopt(name = "jsnark-parser", about = "Jsnark file parsing")]
struct ExampleOptions {
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

    /// Verbose logging and reporting.
    #[structopt(short = "v", long = "verbose")]
    verbose: bool,
}
