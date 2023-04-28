// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use models::jsnark_arith_parser::JsnarkArithReaderParser;
use models::jsnark_wire_parser::JsnarkWireReaderParser;
use models::utils::print_vec;
use winter_math::fields::f128::BaseElement;

fn main() {
    let options = ExampleOptions::from_args();
    let verbose = options.verbose;

    if verbose {
        println!("Parse files {} {}", options.arith_file, options.wires_file);
    }

    let mut arith_file_parser = JsnarkArithReaderParser::<BaseElement>::new().unwrap();
    arith_file_parser.parse_arith_file(&options.arith_file, verbose);
    let r1cs_instance = arith_file_parser.r1cs_instance;

    let mut wire_file_parser = JsnarkWireReaderParser::<BaseElement>::new().unwrap();
    wire_file_parser.parse_wire_file(&options.wires_file, verbose);
    let wires = wire_file_parser.wires;

    if verbose {
        // r1cs_instance.debug_print_bits_horizontal();
        // r1cs_instance.debug_print_symbolic();
        print_vec(&wires);
        println!("");
    }
}

#[derive(Debug)]
struct ExampleOptions {
    arith_file: String,
    wires_file: String,
    verbose: bool,
}
impl ExampleOptions {
    fn from_args() -> Self {
        ExampleOptions {
            arith_file: "fractal_examples2/jsnark_outputs/sample.arith".to_string(),
            wires_file: "fractal_examples2/jsnark_outputs/sample.wires".to_string(),
            verbose: true,
        }
    }
}
