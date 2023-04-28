use crate::errors::*;
use crate::polynomial_utils;
use fractal_math::FieldElement;
use std::convert::TryInto;

// TODO: Add error checking and throwing

pub(crate) type MatrixDimensions = (usize, usize);
#[derive(Clone, Debug)]
pub(crate) struct Matrix<E: FieldElement> {
    pub name: String,
    pub mat: Vec<Vec<E>>,
    pub dims: MatrixDimensions,
}

pub(crate) fn valid_matrix<E: FieldElement>(
    name: &str,
    matrix: Vec<Vec<E>>,
) -> Result<Matrix<E>, MatrixError> {
    let rows = matrix.len();
    if rows == 0 {
        let dims = (0, 0);
        return Ok(Matrix {
            name: String::from(name),
            mat: matrix,
            dims: dims,
        });
    } else {
        let cols = matrix[0].len();
        for i in 0..rows {
            if matrix[i].len() != cols {
                return Err(MatrixError::InvalidMatrix(String::from(name)));
            }
        }
        let dims = (rows, cols);
        return Ok(Matrix {
            name: String::from(name),
            mat: matrix,
            dims,
        });
    }
}

impl<E: FieldElement> Matrix<E> {
    pub(crate) fn new(name: &str, matrix: Vec<Vec<E>>) -> Result<Self, MatrixError> {
        let valid = valid_matrix(name, matrix);
        match valid {
            Ok(m) => Ok(m),
            Err(e) => Err(e),
        }
    }

    pub(crate) fn get_total_size(&self) -> usize {
        let rows = self.dims.0;
        let cols = self.dims.1;
        let total_size = rows + cols;
        total_size
    }

    pub(crate) fn get_rows(&self) -> usize {
        self.dims.0
    }

    pub(crate) fn get_cols(&self) -> usize {
        self.dims.1
    }

    // TODO: Should throw an exception if the row and col are
    // too large
    pub(crate) fn get_value(&self, row: usize, col: usize) -> E {
        self.mat[row][col]
    }

    pub(crate) fn get_transpose(&self, transpose_name: &str) -> Self {
        let mut new_mat: Vec<Vec<E>> = Vec::new();
        let new_cols = self.get_rows();
        let new_rows = self.get_cols();
        for i in 0..new_rows {
            let mut new_mat_row = Vec::new();
            for j in 0..new_cols {
                new_mat_row.push(self.get_value(j, i));
            }
            new_mat.push(new_mat_row);
        }
        Matrix {
            name: String::from(transpose_name),
            mat: new_mat,
            dims: (new_rows, new_cols),
        }
    }

    pub(crate) fn get_matrix_star(&self, index_field: Vec<E>) -> Result<Self, MatrixError> {
        let mut mat_star = Vec::new();
        let dom_size = index_field.len().try_into().unwrap();
        for i in 0..self.get_rows() {
            let mut new_row: Vec<E> = Vec::new();
            for j in 0..self.get_cols() {
                let new_elt = self.mat[j][i]
                    * polynomial_utils::compute_derivative_on_single_val(index_field[j], dom_size);
                new_row.push(new_elt);
            }
            mat_star.push(new_row);
        }
        let mut new_name = self.name.clone();
        new_name.push_str("Star");
        Self::new(&new_name, mat_star)
    }
}
