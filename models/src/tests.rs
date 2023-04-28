use fractal_math::smallprimefield;
use winter_math::{
    fields::f128::{self, BaseElement},
    FieldElement,
};

use crate::{
    errors::R1CSError,
    r1cs::{Matrix, R1CS},
};

type SmallFieldElement17 = smallprimefield::BaseElement<17, 3, 4>;

#[test]
fn test_construct_matrix_f128() {
    let m: Result<Matrix<f128::BaseElement>, R1CSError> = make_all_ones_matrix_f128("dummy", 1, 1);
    let matrix = m.unwrap();

    let (r, c) = matrix.dims;
    assert!(r == 1);
    assert!(c == 1);
    for i in 0..1 {
        for j in 0..1 {
            assert!(matrix.mat[i][&j] == f128::BaseElement::ONE);
        }
    }
    assert!(matrix.name == "dummy");
}

#[test]
fn test_construct_matrix_f17() {
    let m: Result<Matrix<SmallFieldElement17>, R1CSError> = make_all_ones_matrix_f17("dummy", 1, 1);
    let matrix = m.unwrap();

    let (r, c) = matrix.dims;
    assert!(r == 1);
    assert!(c == 1);
    for i in 0..1 {
        for j in 0..1 {
            assert!(matrix.mat[i][&j] == SmallFieldElement17::ONE);
        }
    }
    assert!(matrix.name == "dummy");
}

#[test]
fn test_construct_r1cs() {
    let m1 = make_all_ones_matrix_f128("A", 1, 1);
    let matrix_a = m1.unwrap();
    let m2 = make_all_ones_matrix_f128("B", 1, 1);
    let matrix_b = m2.unwrap();
    let m3 = make_all_ones_matrix_f128("C", 1, 1);
    let matrix_c = m3.unwrap();

    let r1cs_instance = R1CS::new(matrix_a, matrix_b, matrix_c);
    assert!(r1cs_instance.is_ok());
    let unwrapped = r1cs_instance.unwrap();
    assert!(
        unwrapped.num_rows() == 1,
        "Matrix should have had 1 row, instead got {}",
        unwrapped.num_rows()
    );
    assert!(
        unwrapped.num_cols() == 1,
        "Matrix should have had 1 col, instead got {}",
        unwrapped.num_cols()
    );
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
