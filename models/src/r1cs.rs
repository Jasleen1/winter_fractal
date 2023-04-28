use rustc_hash::FxHashMap;

use winter_math::StarkField;

use crate::errors::*;
use crate::utils::{print_vec, print_vec_bits};

pub type MatrixDimensions = (usize, usize);
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<E: StarkField> {
    pub name: String,
    pub mat: Vec<FxHashMap<usize, E>>,
    pub dims: MatrixDimensions,
}

pub fn valid_matrix<E: StarkField>(
    name: &str,
    matrix: Vec<Vec<E>>,
) -> Result<Matrix<E>, R1CSError> {
    let rows = matrix.len();
    if rows == 0 {
        let dims = (0, 0);
        return Ok(Matrix {
            name: String::from(name),
            mat: compress_matrix(matrix),
            dims: dims,
        });
    } else {
        let cols = matrix[0].len();
        for i in 0..rows {
            if matrix[i].len() != cols {
                return Err(R1CSError::InvalidMatrix(String::from(name)));
            }
        }
        let dims = (rows, cols);
        return Ok(Matrix {
            name: String::from(name),
            mat: compress_matrix(matrix),
            dims,
        });
    }
}

fn compress_matrix<E: StarkField>(longer_matrix: Vec<Vec<E>>) -> Vec<FxHashMap<usize, E>> {
    let mut out_matrix = Vec::<FxHashMap::<usize, E>>::new();
    for row in longer_matrix.iter() {
        let mut compressed_row = FxHashMap::<usize, E>::default();
        for (loc, elt)  in row.iter().enumerate() {
            if *elt != E::ZERO {
                compressed_row.insert(loc, *elt);
            }
        }
        out_matrix.push(compressed_row);
    }
    out_matrix
}

impl<E: StarkField> Matrix<E> {
    pub fn new(name: &str, matrix: Vec<Vec<E>>) -> Result<Self, R1CSError> {
        Ok(valid_matrix(name, matrix)?)
    }

    pub fn num_rows(&self) -> usize {
        self.dims.0
    }

    pub fn num_cols(&self) -> usize {
        self.dims.1
    }

    pub fn get_total_size(&self) -> usize {
        let rows = self.dims.0;
        let cols = self.dims.1;
        let total_size = rows + cols;
        return total_size;
    }

    // L0 norm, number of nonzero elements.
    pub fn l0_norm(&self) -> usize {
        let l0_norm = self.mat.iter().fold(0, |a, row| {
            a + row.len()
        });
        l0_norm
    }

    pub fn dot(&self, vec: &Vec<E>) -> Vec<E> {
        // self.mat
        //     .iter()
        //     .map(|a| {
        //         a.iter()
        //             .zip(vec.iter())
        //             .map(|(x, y)| x.mul(*y))
        //             .fold(E::ZERO, |sum, i| sum.add(i))
        //     })
        //     .collect()
        self.mat
            .iter()
            .map(|a| {
                a.iter()
                    .map(|(&loc, val)| val.mul(vec[loc]))
                    .fold(E::ZERO, |sum, i| sum.add(i))
            })
            .collect()
    }

    pub fn define_cols(&mut self, num_cols: usize) {
        assert!(
            self.dims.1 <= num_cols,
            "Attempted to reduce number of columns."
        );
        self.dims.1 = num_cols;
        // for row in &mut self.mat {
        //     row.resize(num_cols, E::ZERO);
        // }
    }

    pub fn define_rows(&mut self, num_rows: usize) {
        assert!(
            self.dims.0 <= num_rows,
            "Attempted to reduce number of rows."
        );
        let zero_row = FxHashMap::<usize, E>::default();
        let num_to_pad = num_rows - self.dims.0;
        for _ in 0..num_to_pad {
            self.mat.push(zero_row.clone());
        }
        self.dims.0 = num_rows;
    }

    fn compress_row(new_row: &Vec<E>) -> FxHashMap<usize, E> {
        let mut new_row_comp = FxHashMap::<usize, E>::default();
        for (loc, val) in new_row.iter().enumerate() {
            if *val != E::ZERO {
                new_row_comp.insert(loc, *val);
            }
        }
        new_row_comp
    } 

    pub fn add_row(&mut self, new_row: &Vec<E>) {
        if new_row.len() != self.dims.1 {
            // FIXME: add error handling
        }
        self.mat.push(Self::compress_row(&new_row));
        self.dims.0 = self.dims.0 + 1;
    }

    pub fn pad_rows(&mut self, new_row_count: usize) {
        let new_row = FxHashMap::<usize, E>::default();
        for _ in 0..new_row_count {
            self.mat.push(new_row.clone());
            self.dims.0 = self.dims.0 + 1;
        }
    }

    pub fn pad_power_two(&mut self) {
        let rows = self.dims.0;
        let cols = self.dims.1;

        self.define_cols(cols.next_power_of_two());
        self.pad_rows(rows.next_power_of_two() - rows);
    }

    pub fn make_square(&mut self) {
        let rows = self.dims.0;
        let cols = self.dims.1;
        let square_dim = rows.max(cols);

        self.define_cols(square_dim);
        self.define_rows(square_dim);
    }

    pub fn debug_print(&self) {
        println!("{}", self.name);
        for row in &self.mat {
            println!("{:?}", row);
            println!("");
        }
    }

    fn row_to_vec(&self, row: &FxHashMap<usize, E>) -> Vec<E> {
        let mut vec_form = vec![E::ZERO; self.dims.1];
        row.iter()
        .map(|(&loc, val)| vec_form[loc] = *val);
        vec_form
    }

    pub fn debug_print_bits(&self) {
        println!("{}", self.name);
        for row in &self.mat {
            print_vec_bits(&self.row_to_vec(row));
            println!("");
        }
    }
}

pub(crate) fn create_empty_matrix<E: StarkField>(name: String) -> Matrix<E> {
    Matrix {
        name,
        mat: Vec::<FxHashMap<usize, E>>::new(),
        dims: (0, 0),
    }
}

pub(crate) fn create_empty_r1cs<E: StarkField>() -> Result<R1CS<E>, R1CSError> {
    let matrix_a = create_empty_matrix("A".to_string());
    let matrix_b = create_empty_matrix("B".to_string());
    let matrix_c = create_empty_matrix("C".to_string());
    R1CS::new(matrix_a, matrix_b, matrix_c)
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct R1CS<E: StarkField> {
    #[allow(non_snake_case)]
    pub A: Matrix<E>,
    #[allow(non_snake_case)]
    pub B: Matrix<E>,
    #[allow(non_snake_case)]
    pub C: Matrix<E>,
}

// TODO Might want to change this to include checks for A, B and C.
impl<E: StarkField> R1CS<E> {
    pub fn new(
        matrix_a: Matrix<E>,
        matrix_b: Matrix<E>,
        matrix_c: Matrix<E>,
    ) -> Result<Self, R1CSError> {
        let valid = valid_r1cs(&matrix_a, &matrix_b, &matrix_c);
        match valid {
            Ok(_) => Ok(R1CS {
                A: matrix_a,
                B: matrix_b,
                C: matrix_c,
            }),
            Err(e) => Err(e),
        }
    }

    pub fn num_rows(&self) -> usize {
        self.A.dims.0
    }

    pub fn num_cols(&self) -> usize {
        self.A.dims.1
    }

    pub fn max_num_nonzero(&self) -> usize {
        self.A.l0_norm().max(self.B.l0_norm()).max(self.C.l0_norm())
    }

    pub fn get_a(&mut self) -> &mut Matrix<E> {
        &mut self.A
    }

    pub fn get_b(&mut self) -> &mut Matrix<E> {
        &mut self.B
    }

    pub fn get_c(&mut self) -> &mut Matrix<E> {
        &mut self.C
    }

    pub fn set_cols(&mut self, num_cols: usize) {
        self.A.define_cols(num_cols);
        self.B.define_cols(num_cols);
        self.C.define_cols(num_cols);
    }

    pub fn add_rows(&mut self, new_row_a: Vec<E>, new_row_b: Vec<E>, new_row_c: Vec<E>) {
        self.A.add_row(&new_row_a);
        self.B.add_row(&new_row_b);
        self.C.add_row(&new_row_c);
    }

    pub fn pad_power_two(&mut self) {
        self.A.pad_power_two();
        self.B.pad_power_two();
        self.C.pad_power_two();
    }

    pub fn make_square(&mut self) {
        self.A.make_square();
        self.B.make_square();
        self.C.make_square();
    }

    pub fn debug_print(&self) {
        println!("Dimensions: {} {}", self.A.dims.0, self.A.dims.1);
        self.A.debug_print();
        self.B.debug_print();
        self.C.debug_print();
    }

    pub fn debug_print_bits(&self) {
        println!("Dimensions: {} {}", self.A.dims.0, self.A.dims.1);
        self.A.debug_print_bits();
        self.B.debug_print_bits();
        self.C.debug_print_bits();
    }

    // pub fn debug_print_bits_horizontal(&self) {
    //     let num_rows = self.A.dims.0;
    //     if num_rows == 0 {
    //         println!("No rows in the matrix!");
    //         return;
    //     }
    //     for row_idx in 0..num_rows {
    //         print_vec_bits(&self.A.mat[row_idx]);
    //         print!("  ");
    //         print_vec_bits(&self.B.mat[row_idx]);
    //         print!("  ");
    //         print_vec_bits(&self.C.mat[row_idx]);
    //         println!("");
    //     }
    // }

    fn debug_print_row_symbolic(&self, row: &Vec<E>) {
        let mut first = true;
        for col_idx in 0..row.len() {
            let elt = row[col_idx];
            if elt != E::ZERO {
                if first {
                    first = false;
                } else {
                    print!(" + ");
                }
                if col_idx == 0 {
                    print!("{}", elt);
                } else {
                    if elt == E::ONE {
                        print!("v{}", col_idx);
                    } else if elt == E::ONE.neg() {
                        print!("-v{}", col_idx);
                    } else {
                        print!("{} v{}", elt, col_idx);
                    }
                }
            }
        }
        if first {
            // Never encountered a nonzero. So the sum itself is 0.
            print!("0");
        }
    }

    // pub fn debug_print_symbolic(&self) {
    //     let num_rows = self.A.dims.0;
    //     if num_rows == 0 {
    //         println!("No rows in the matrix!");
    //         return;
    //     }
    //     for row_idx in 0..num_rows {
    //         print!("(");
    //         self.debug_print_row_symbolic(&self.A.mat[row_idx]);
    //         print!(")  (");
    //         self.debug_print_row_symbolic(&self.B.mat[row_idx]);
    //         print!(") == ");
    //         self.debug_print_row_symbolic(&self.C.mat[row_idx]);
    //         println!("");
    //     }
    // }
}

// TODO: indexed R1CS consisting of 3 indexed matrices

// TODO: Add error here

pub fn valid_r1cs<E: StarkField>(
    a: &Matrix<E>,
    b: &Matrix<E>,
    c: &Matrix<E>,
) -> Result<bool, crate::errors::R1CSError> {
    let a_dims = a.dims;
    let b_dims = b.dims;
    let c_dims = c.dims;
    if b_dims != a_dims {
        Err(R1CSError::MatrixSizeMismatch(
            a.name.clone(),
            b.name.clone(),
        ))
    } else if c_dims != a_dims {
        Err(R1CSError::MatrixSizeMismatch(
            a.name.clone(),
            c.name.clone(),
        ))
    } else {
        Ok(true)
    }
}
#[cfg(test)]
mod localtests {
    use winter_math::{
        fields::f128::{self, BaseElement},
        FieldElement,
    };

    use super::Matrix;
    use crate::errors::R1CSError;
    #[test]
    fn test_matrix_dot() {
        let matrix = make_all_ones_matrix_f128("jim", 2, 2).unwrap();
        let output = matrix.dot(&vec![BaseElement::new(2u128), BaseElement::new(5u128)]);
        let expected = vec![BaseElement::new(7u128); 2];
        for i in 0..output.len() {
            assert_eq!(output[i], expected[i]);
        }
    }
    #[test]
    fn test_matrix_dot_2() {
        let mut mat = Vec::new();
        let first_row = vec![BaseElement::new(3u128), BaseElement::new(2u128)];
        let second_row = vec![BaseElement::new(4u128), BaseElement::new(5u128)];
        mat.push(first_row);
        mat.push(second_row);
        let matrix = Matrix::new("steve", mat).unwrap();
        let output = matrix.dot(&vec![BaseElement::new(7u128), BaseElement::new(11u128)]);
        let expected = vec![BaseElement::new(43u128), BaseElement::new(83u128)];
        for i in 0..output.len() {
            assert_eq!(output[i], expected[i]);
        }
    }

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
}
