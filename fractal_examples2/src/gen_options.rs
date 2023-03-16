// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use core::num;
use std::cmp::max;

use fractal_indexer::index::get_max_degree;
use fractal_proofs::FriOptions;
use fractal_prover::FractalOptions;

use fractal_indexer::{
    index::{build_index_domains, Index, IndexParams},
    indexed_matrix::index_matrix,
    snark_keys::*,
};

use models::jsnark_arith_parser::JsnarkArithReaderParser;
use models::jsnark_wire_parser::JsnarkWireReaderParser;

use winter_crypto::hashers::{Rp64_256, Blake3_256};
use winter_crypto::ElementHasher;

use winter_math::fields::f64::BaseElement;
use winter_math::fields::QuadExtension;
use winter_math::utils;
use winter_math::FieldElement;
use winter_math::StarkField;

pub fn get_example_setup<
B: StarkField,
E: FieldElement<BaseField = B>,
H: ElementHasher + ElementHasher<BaseField = B>,
>() -> (FractalOptions<B>, ProverKey<H,B>, VerifierKey<H,B>) {
    let mut options = ExampleOptions::from_args();
    options.verbose = true;
    if options.verbose {
        println!(
            "Arith file {}, wire value file {}",
            options.arith_file, options.wires_file
        );
    }

    // call orchestrate_r1cs_example
    files_to_setup_outputs::<B, E, H, 1>(
    //orchestrate_r1cs_example::<BaseElement, BaseElement, Blake3_256<BaseElement>, 1>(
        &options.arith_file,
        &options.wires_file,
        options.verbose,
    )
}

fn files_to_setup_outputs<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
>(
    arith_file: &str,
    wire_file: &str,
    verbose: bool,
) -> (FractalOptions<B>, ProverKey<H,B>, VerifierKey<H,B>){
    let mut arith_parser = JsnarkArithReaderParser::<B>::new().unwrap();
    arith_parser.parse_arith_file(&arith_file, verbose);
    let r1cs = arith_parser.clone_r1cs();

    let mut wires_parser = JsnarkWireReaderParser::<B>::new().unwrap();
    wires_parser.parse_wire_file(&wire_file, verbose);
    let wires = wires_parser.wires;
    println!("Len wires = {}", wires.len());
    // 0. Compute num_non_zero by counting max(number of non-zero elts across A, B, C).

    // let num_input_variables = r1cs.clone().num_cols();
    // let num_constraints = r1cs.clone().num_rows();
    // let num_non_zero = max(max(r1cs.A.l0_norm(), r1cs.B.l0_norm()), r1cs.C.l0_norm());
    // 1. Index this R1CS
    let num_input_variables = r1cs.num_cols().next_power_of_two();
    let num_non_zero = r1cs.max_num_nonzero().next_power_of_two();
    let num_constraints =
        max(max(r1cs.A.l0_norm(), r1cs.B.l0_norm()), r1cs.C.l0_norm()).next_power_of_two();
    let max_degree = get_max_degree(num_input_variables, num_non_zero, num_constraints);
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

    let index_domains = build_index_domains::<B>(index_params.clone());
    println!("build index domains");
    let indexed_a = index_matrix::<B>(&r1cs.A, &index_domains);
    let indexed_b = index_matrix::<B>(&r1cs.B, &index_domains);
    let indexed_c = index_matrix::<B>(&r1cs.C, &index_domains);
    println!("indexed matries");
    // This is the index i.e. the pre-processed data for this r1cs
    let index = Index::new(index_params.clone(), indexed_a, indexed_b, indexed_c);

    let (prover_key, verifier_key) = generate_prover_and_verifier_keys::<H, B, N>(index).unwrap();

    // TODO: the IndexDomains should already guarantee powers of two, so why add extraneous bit or use next_power_of_two?

    let degree_fs = r1cs.num_cols();
    let size_subgroup_h = index_domains.h_field.len().next_power_of_two();
    let size_subgroup_k = index_domains.k_field.len().next_power_of_two();

    let evaluation_domain =
        utils::get_power_series(index_domains.l_field_base, index_domains.l_field_len);

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
        summing_domain,
        evaluation_domain,
        h_domain,
        eta,
        eta_k,
        fri_options,
        num_queries,
    };
    (options, prover_key, verifier_key)
}

#[derive(Debug)]
struct ExampleOptions {
    arith_file: String,
    wires_file: String,
    verbose: bool,
}
impl ExampleOptions{
    fn from_args() -> Self{
        ExampleOptions { arith_file: "fractal_examples2/jsnark_outputs/sample.arith".to_string(), wires_file: "fractal_examples2/jsnark_outputs/sample.wires".to_string(), verbose: true }
    } 
}