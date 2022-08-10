use math::fields::f128::BaseElement;

use crate::arith_parser_example::reading_arith;

#[test]
fn test_arith_parser() {
    let r1cs = reading_arith::<BaseElement>("src/sample.arith", true);
    println!("r1cs dimensions: {:?}", r1cs.A.dims);
}
