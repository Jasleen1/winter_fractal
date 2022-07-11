use crate::{errors::MatrixError, matrix_utils::*, SmallFieldElement17};
use fractal_math::{FieldElement, StarkField};

#[test]
fn test_matrix_star() {
    let original_matrix = make_all_ones_matrix_f17("test", 2, 2).unwrap();
    let field_base = SmallFieldElement17::get_root_of_unity(2);
    unsafe {
        let field = SmallFieldElement17::get_power_series(field_base, 2);
        let matrix_star = original_matrix.get_matrix_star(field.clone()).unwrap();
        let expected = vec![vec![2, 9], vec![2, 9]];
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(
                    matrix_star.get_value(i, j),
                    SmallFieldElement17::new(expected[i][j])
                );
            }
        }
    }
}

fn make_all_ones_matrix_f17(
    matrix_name: &str,
    rows: usize,
    cols: usize,
) -> Result<Matrix<SmallFieldElement17>, MatrixError> {
    let mut mat = Vec::new();
    let ones_row = vec![SmallFieldElement17::ONE; cols];
    for _i in 0..rows {
        mat.push(ones_row.clone());
    }
    Matrix::new(matrix_name, mat)
}
