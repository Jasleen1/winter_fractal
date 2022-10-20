use crate::matrix_utils::*;
use fractal_math::FieldElement;
use std::convert::TryInto;
// TODO: Add error checking and throwing
/**
 * This is equivalent to computing v_H(X) for a multiplicative coset
 * H = eta * H_0 for a multiplicative subgroup H_0 of order dom_size
 * Note that v_H(X) = X^dom_size - eta^dom_size. If eta = 1 we're
 * in a multiplicative subgroup itself.
 **/
pub fn compute_vanishing_poly<E: FieldElement>(x: E, eta: E, dom_size: usize) -> E {
    let power_u64: u64 = dom_size.try_into().unwrap();
    let power = E::PositiveInteger::from(power_u64);
    x.exp(power) - eta.exp(power)
}

/**
 * Compute vanishing polynomial for a multiplicative subgroup. Same as above with
 * eta = ONE.
 **/
pub fn vanishing_poly_for_mult_subgroup<E: FieldElement>(x: E, dom_size: usize) -> E {
    compute_vanishing_poly(x, E::ONE, dom_size)
}

// The derivative is calculated as |H| x^{|H|-1}.
// This is equivalent to computing u_H(X, X) for a multiplicative coset H
// of order dom_size = |H|.
pub fn compute_derivative_on_single_val<E: FieldElement>(x: E, dom_size: u128) -> E {
    let dom_size_coeff = E::from(dom_size);
    let power_u64: u64 = (dom_size - 1).try_into().unwrap();
    let power = E::PositiveInteger::from(power_u64);
    // println!("   deriv: {}    = h: {}   x: {:?}  ** h1: {:?}",
    //     dom_size_coeff * x.exp(power), dom_size_coeff, x, power);
    dom_size_coeff * x.exp(power)
}

// Represents a binomial, i.e. a polynomial in two variables X and Y.
// The (i, j)th element of this binomial is the coefficient of X^i * Y^j
pub type BivariatePoly<E> = Vec<Vec<E>>;

pub fn compute_binomial_on_x<E: FieldElement>(bivariate: BivariatePoly<E>, x_val: E) -> Vec<E> {
    // Given a BivariatePoly, computes a monomial in Y, obtained by evaluating bivariate(x_val, Y)
    // represented by the output vector
    let bivariate_as_matrix = Matrix::new("binomial", bivariate).unwrap();
    let transposed_bivariate = bivariate_as_matrix.get_transpose("transposed_binomial");
    compute_binomial_on_y(transposed_bivariate.mat, x_val)
}

pub fn compute_binomial_on_y<E: FieldElement>(bivariate: BivariatePoly<E>, y_val: E) -> Vec<E> {
    // Given a BivariatePoly, computes a monomial in X, obtained by evaluating bivariate(X, y_val)
    // represented by the output vector
    // Note that since bivariate[i][j] is the coefficient of X^i Y^j, technically,
    // bivariate[i] * (Y^j)_j is the coeffient of X^i
    // Hence, evaluating the polynomial bivariate[i] on y_val gives the coeff of X^i
    let mut x_coeffs = Vec::new();
    for i in 0..bivariate.len() {
        x_coeffs.push(fractal_math::polynom::eval(&bivariate[i], y_val));
    }
    x_coeffs
}

pub fn pad_with_zeroes<E: FieldElement>(poly: &mut Vec<E>, total_len: usize) {
    if total_len <= poly.len() {
        return;
    }
    let diff = total_len - poly.len();
    for _ in 0..diff {
        poly.push(E::ZERO);
    }
}

pub fn get_to_degree_size<E: FieldElement>(poly: &mut Vec<E>) {
    let len_poly = poly.len();
    let mut count = len_poly - 1;
    while count > 0 && poly[count] == E::ZERO {
        poly.pop();
        count = count - 1;
    }
}

pub fn get_complementary_poly<E: FieldElement>(
    current_degree: usize,
    desired_degree: usize,
) -> Vec<E> {
    assert!(desired_degree >= current_degree);
    let comp_deg = desired_degree - current_degree;
    let mut out_poly = vec![E::ZERO; comp_deg];
    out_poly.push(E::ONE);
    out_poly[0] = E::ONE;
    out_poly
}

pub fn get_randomized_complementary_poly<E: FieldElement>(
    current_degree: usize,
    desired_degree: usize,
    alpha: E,
    beta: E,
) -> Vec<E> {
    assert!(desired_degree >= current_degree);
    let comp_deg = desired_degree - current_degree;
    let mut out_poly = vec![E::ZERO; comp_deg];
    out_poly.push(alpha);
    out_poly[0] = beta;
    out_poly
}
