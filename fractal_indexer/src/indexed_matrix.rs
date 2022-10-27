// TODO Implement an indexed matrix struct with row, col and val and evaluations on an eval domain

use std::convert::TryInto;

// TODO: This implementation assumes all matrices are square and all inputs are public, ie no witness. Update to accomodate this.
use crate::index::*;
use fractal_math::polynom;
use fractal_utils::polynomial_utils;
use models::r1cs::*;
use winter_math::{fft, StarkField};

#[derive(Clone, Debug)]
pub struct IndexedMatrix<E: StarkField> {
    // i_field: Vec<E>, // This is the subfield of H (below) based on which non-witness inputs are indexed
    // h_field: Vec<E>, // This is the field which indexes the matrices. Total length of input+witness = |H|.
    // k_field: Vec<E>, // This is the field used for the sparse matrix representation
    // L_field: Vec<E>, // This is the lde domain as in the fractal paper
    pub matrix: Matrix<E>,

    pub row_poly: Vec<E>,
    pub col_poly: Vec<E>,
    pub val_poly: Vec<E>,

    pub row_evals_on_l: Vec<E>,
    pub col_evals_on_l: Vec<E>,
    pub val_evals_on_l: Vec<E>,
}

// TODO: Implement commitment for the index to be used as part of the verifier key
// TODO: Add error checking
impl<E: StarkField> IndexedMatrix<E> {
    pub fn new(mat: &Matrix<E>, domains: &IndexDomains<E>) -> Self {
        index_matrix(mat, domains)
    }
}

// TODO where should we save the global domain and other values?
// Also, should the new indexed matrix be generated using something here or in the new
// function for Indexed Matrix?
// QUESTION: Should the IndexDomain struct also depend on E?
pub fn index_matrix<E: StarkField>(
    mat: &Matrix<E>,
    index_domains: &IndexDomains<E>,
) -> IndexedMatrix<E> {
    let h_size = index_domains.h_field.len().try_into().unwrap();
    let l_size = index_domains.l_field_len;
    let h_field = index_domains.h_field.clone();
    let num_rows = mat.dims.0;
    let num_cols = mat.dims.1;

    let k_field_size = index_domains.k_field.len();

    // K is chosen large enough to enumerate the nonzero elements of M.
    // H is chosen large enough to enumerate the rows (or cols) of M.
    // index : K -> H x H x L
    // We need up to k_field_size entries.
    let mut row_elts = vec![h_field[0]; k_field_size];
    let mut col_elts = vec![h_field[0]; k_field_size];
    let mut val_elts = vec![E::ZERO; k_field_size];

    let mut count = 0;

    for r_int in 0..num_rows {
        for c_int in 0..num_cols {
            if mat.mat[r_int][c_int] == E::ZERO {
                continue;
            }
            let c = h_field[c_int];
            let r = h_field[r_int];

            row_elts[count] = c;
            col_elts[count] = r;
            val_elts[count] = mat.mat[r_int][c_int]
                * polynomial_utils::compute_derivative_on_single_val(r, h_size)
                / (compute_derivative(c, h_size) * compute_derivative(r, h_size));
            count += 1;
        }
    }

    let inv_twiddles_k_elts = index_domains.inv_twiddles_k_elts.clone();

    // interpolate row_elts into a polynomial
    fft::interpolate_poly_with_offset(&mut row_elts, &inv_twiddles_k_elts, index_domains.eta_k);

    // interpolate col_elts into a polynomial
    fft::interpolate_poly_with_offset(&mut col_elts, &inv_twiddles_k_elts, index_domains.eta_k);

    // interpolate val_elts into a polynomial
    fft::interpolate_poly_with_offset(&mut val_elts, &inv_twiddles_k_elts, index_domains.eta_k);

    let twiddles_l_elts = index_domains.twiddles_l_elts.clone();

    // evaluate row_elts polynomial over l
    let mut row_evaluations = vec![E::ZERO; l_size];
    row_evaluations[..k_field_size].copy_from_slice(&row_elts);
    fft::evaluate_poly(&mut row_evaluations, &twiddles_l_elts);

    // evaluate col_elts polynomial over l
    let mut col_evaluations = vec![E::ZERO; l_size];
    col_evaluations[..k_field_size].copy_from_slice(&col_elts);
    fft::evaluate_poly(&mut col_evaluations, &twiddles_l_elts);

    // evaluate row_elts polynomial over l
    let mut val_evaluations = vec![E::ZERO; l_size];
    val_evaluations[..k_field_size].copy_from_slice(&val_elts);
    fft::evaluate_poly(&mut val_evaluations, &twiddles_l_elts);

    IndexedMatrix {
        matrix: mat.clone(),
        row_poly: row_elts,
        col_poly: col_elts,
        val_poly: val_elts,
        row_evals_on_l: row_evaluations,
        col_evals_on_l: col_evaluations,
        val_evals_on_l: val_evaluations,
    }
}

/// ***************  HELPERS *************** \\\

// This is equivalent to computing u_H(X, X) for a multiplicative group H
// of order dom_size = |H|.
// TODO: Add error checking and throwing
pub fn compute_derivative<E: StarkField>(x: E, dom_size: u128) -> E {
    let dom_size_coeff = E::from(dom_size);
    let power_u64: u64 = (dom_size - 1).try_into().unwrap();
    let power = E::PositiveInteger::from(power_u64);
    dom_size_coeff * x.exp(power)
}
