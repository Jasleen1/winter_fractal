use std::vec;

use crate::{index::*, *};
use indexed_matrix::IndexedMatrix;
use math::{fields::f128::BaseElement, FieldElement, StarkField};
use models::r1cs::Matrix;
use models::{errors::R1CSError, r1cs::*};

type SmallFieldElement17 = math::fields::smallprimefield::BaseElement<17, 3, 4>;

#[test]
fn test_indexing() {
    let m1 = make_all_ones_matrix_f128("A", 2, 2);
    let matrix_a = m1.unwrap();
    let m2 = make_all_ones_matrix_f128("B", 2, 2);
    let matrix_b = m2.unwrap();
    let m3 = make_all_ones_matrix_f128("C", 2, 2);
    let matrix_c = m3.unwrap();

    // QUESTION: For all uses of clone(), below, is there a better way to pass the value? Perhaps as
    // an immutable reference?
    let r1cs_instance_result = R1CS::new(matrix_a, matrix_b, matrix_c);
    let r1cs_instance = r1cs_instance_result.unwrap();
    let params = IndexParams::<BaseElement> {
        num_input_variables: 2,
        num_constraints: 2,
        num_non_zero: 4,
        max_degree: get_max_degree(2, 2, 4),
        eta: BaseElement::ONE
    };
    let domains = build_index_domains(params.clone());
    let indexed_a = IndexedMatrix::new(&r1cs_instance.A, &domains);
    let indexed_b = IndexedMatrix::new(&r1cs_instance.B, &domains);
    let indexed_c = IndexedMatrix::new(&r1cs_instance.C, &domains);
    let index = Index::new(params, indexed_a, indexed_b, indexed_c);
    println!("Index is {:?}", index);
}

#[test]
fn test_domain_building_17() {
    let params = IndexParams::<SmallFieldElement17> {
        num_input_variables: 2,
        num_constraints: 2,
        num_non_zero: 4,
        max_degree: get_max_degree(2, 2, 4),
        eta: SmallFieldElement17::ONE,
    };
    let domains = build_primefield_index_domains(params.clone());
    let i_field_base = domains.i_field_base;
    let k_field_base = domains.k_field_base;
    let h_field_base = domains.h_field_base;
    let l_field_base = domains.l_field_base;
    assert_eq!(
        i_field_base,
        SmallFieldElement17::new(16),
        "Bad base for i_field"
    );
    assert_eq!(
        k_field_base,
        SmallFieldElement17::new(13),
        "Bad base for k_field"
    );
    assert_eq!(
        h_field_base,
        SmallFieldElement17::new(16),
        "Bad base for h_field"
    );
    // Temp using l_field as 2 * ..., so really expect 9.
    //assert_eq!(l_field_base, SmallFieldElement17::new(3), "Bad base for l_field");
    assert_eq!(
        l_field_base,
        SmallFieldElement17::new(9),
        "Bad base for l_field"
    );
}

#[test]
fn test_getting_roots_17() {
    let test_root_16 = SmallFieldElement17::get_root_of_unity(4);
    assert_eq!(test_root_16, SmallFieldElement17::new(3));
    let test_root_8 = SmallFieldElement17::get_root_of_unity(3);
    assert_eq!(test_root_8, SmallFieldElement17::new(9));
    let test_root_2 = SmallFieldElement17::get_root_of_unity(1);
    assert_eq!(test_root_2, SmallFieldElement17::new(16));
}

#[test]
fn test_single_indexed_matrix_17() {
    let m1 = make_all_ones_matrix_f17("A", 2, 2);
    let matrix_a = m1.unwrap();
    let params = IndexParams::<SmallFieldElement17> {
        num_input_variables: 2,
        num_constraints: 2,
        num_non_zero: 4,
        max_degree: get_max_degree(2, 2, 4),
        eta: SmallFieldElement17::ONE,
    };
    let domains = build_index_domains(params.clone());
    println!("Domains {:?}", domains);
    let indexed_a = IndexedMatrix::new(&matrix_a, &domains);
    println!("Indexed a is {:?}", indexed_a);

    let row_poly = indexed_a.row_poly;
    let col_poly = indexed_a.col_poly;
    let expected_row_poly = vec![0, 0, 1, 0];
    let expected_col_poly = vec![0, 11, 0, 7];
    println!("Row computed: {:?}", row_poly);
    println!("Col computed: {:?}", col_poly);
    for i in 0..4 {
        assert_eq!(row_poly[i], SmallFieldElement17::new(expected_row_poly[i]));
        assert_eq!(col_poly[i], SmallFieldElement17::new(expected_col_poly[i]));
    }
}

#[test]
fn test_indexing_f17() {
    let m1 = make_all_ones_matrix_f17("A", 2, 2);
    let matrix_a = m1.unwrap();
    let m2 = make_all_ones_matrix_f17("B", 2, 2);
    let matrix_b = m2.unwrap();
    let m3 = make_all_ones_matrix_f17("C", 2, 2);
    let matrix_c = m3.unwrap();

    // QUESTION: For all uses of clone(), below, is there a better way to pass the value? Perhaps as
    // an immutable reference?
    let r1cs_instance_result = R1CS::new(matrix_a, matrix_b, matrix_c);
    let r1cs_instance = r1cs_instance_result.unwrap();
    let params = IndexParams::<SmallFieldElement17> {
        num_input_variables: 2,
        num_constraints: 2,
        num_non_zero: 4,
        max_degree: get_max_degree(2, 2, 4),
        eta: SmallFieldElement17::ONE,
    };
    let domains = build_primefield_index_domains(params.clone());
    let indexed_a = IndexedMatrix::new(&r1cs_instance.A, &domains);
    let indexed_b = IndexedMatrix::new(&r1cs_instance.B, &domains);
    let indexed_c = IndexedMatrix::new(&r1cs_instance.C, &domains);
    let index = Index::new(params, indexed_a, indexed_b, indexed_c);
    println!("Index is {:?}", index);
}

/// ***************  HELPERS *************** \\\
fn make_all_ones_matrix_f128(
    matrix_name: &str,
    rows: usize,
    cols: usize,
) -> Result<Matrix<BaseElement>, R1CSError> {
    let mut mat = Vec::new();
    let ones_row = vec![BaseElement::ONE; cols];
    for _i in 0..rows {
        mat.push(ones_row.clone());
    }
    Matrix::new(matrix_name, mat)
}

fn make_all_ones_matrix_f17(
    matrix_name: &str,
    rows: usize,
    cols: usize,
) -> Result<Matrix<SmallFieldElement17>, R1CSError> {
    let mut mat = Vec::new();
    let ones_row = vec![SmallFieldElement17::ONE; cols];
    for _i in 0..rows {
        mat.push(ones_row.clone());
    }
    Matrix::new(matrix_name, mat)
}
