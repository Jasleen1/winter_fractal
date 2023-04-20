use std::{cmp::max, marker::PhantomData};

// TODO: This class will include the indexes of 3 matrices
// Should domain info be in here or in a separate class?
use winter_math::{fft, FieldElement, StarkField}; // utils

type SmallFieldElement17 = fractal_math::smallprimefield::BaseElement<17, 3, 4>;

use crate::indexed_matrix::IndexedMatrix;
use models::r1cs::R1CS;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IndexParams<B: StarkField> {
    pub num_input_variables: usize,
    // num_witness_variables: usize,
    pub num_constraints: usize,
    pub num_non_zero: usize,
    pub max_degree: usize,
    pub eta: B,
    pub eta_k: B,
}
#[derive(Clone, Debug)]
pub struct Index<B: StarkField, E: FieldElement<BaseField = B>> {
    pub params: IndexParams<B>,
    pub indexed_a: IndexedMatrix<B, E>,
    pub indexed_b: IndexedMatrix<B, E>,
    pub indexed_c: IndexedMatrix<B, E>,
}

impl<B: StarkField, E: FieldElement<BaseField = B>> Index<B, E> {
    pub fn new(
        params: IndexParams<B>,
        indexed_a: IndexedMatrix<B, E>,
        indexed_b: IndexedMatrix<B, E>,
        indexed_c: IndexedMatrix<B, E>,
    ) -> Self {
        Index {
            params,
            indexed_a: indexed_a,
            indexed_b: indexed_b,
            indexed_c: indexed_c,
        }
    }
}

/// QUESTION: Currently IndexDomains is implemented over a generic FieldElement trait.
/// but do we want to keep it this way, since below the actual implementation to generate
/// indices is BaseElement
#[derive(Clone, Debug)]
pub struct IndexDomains<B: StarkField, E: FieldElement<BaseField = B>> {
    pub i_field_base: B,
    pub k_field_base: B, // Generate sufficiently large set to enumerate all nonzero matrix entries.
    pub h_field_base: B, // Generate sufficiently large set to enumerate each row or each column.
    pub l_field_base: B, // For Reed Solomon code.
    pub i_field: Vec<B>,
    pub k_field: Vec<B>,
    pub h_field: Vec<B>,
    pub l_field_len: usize,
    pub inv_twiddles_k_elts: Vec<B>,
    pub twiddles_l_elts: Vec<B>,
    pub eta: B,
    pub eta_k: B,
    pub phantom_e: PhantomData<E>,
}

/// ***************  HELPERS *************** \\\

// Currently assuming that
// 1. All the inputs to this function are powers of 2
// 2. num_input_variables is the number of inputs and num_input_variables + num_witnesses = num_constraints
// 3. 2, above implies that the matrices are all square.
/// QUESTION: This is currently built using BaseField because the trait has no generic function for
/// getting generators of a certain order. I think this would require some re-structuring.
/// Perhaps we can add a function "get_subgroup_of_size" or "get_generator_of_order"
/// Generators are needed here since we'll need those for FFT-friendly subgroups anyway.
#[cfg_attr(feature = "flame_it", flame)]
pub fn build_index_domains<B: StarkField, E: FieldElement<BaseField = B>>(
    params: IndexParams<B>,
) -> IndexDomains<B, E> {
    let num_input_variables = params.num_input_variables;
    let num_constraints = params.num_constraints;
    let num_non_zero = params.num_non_zero;
    let max_degree = params.max_degree;

    // Validate inputs.
    let ntpow2 = { |x: usize| x > 1 && (x & (x - 1) == 0) };
    assert!(
        ntpow2(num_input_variables),
        "num_input_variables {} must be nontriv power of two",
        num_input_variables
    );
    assert!(
        ntpow2(num_constraints),
        "num_constraints {} must be nontriv power of two",
        num_constraints
    );
    assert!(
        ntpow2(num_non_zero),
        "num_non_zero {} must be nontriv power of two",
        num_non_zero
    );

    assert!(
        ntpow2(max_degree),
        "max_degree {} must be nontriv power of two",
        max_degree
    );

    // Need to encode a subset of H field: indices of inputs.
    let i_field_size = num_input_variables;

    // Need to enumerate each row or each column.
    // Input variables are across the columns, constraints are across the rows.
    let h_field_size = std::cmp::max(num_input_variables, num_constraints);

    // Need to enumerate each nonzero matrix entry.
    let k_field_size = num_non_zero;

    // Find elements in F which generate each subfield.
    let i_field_base = B::get_root_of_unity(i_field_size.trailing_zeros());
    let k_field_base = B::get_root_of_unity(k_field_size.trailing_zeros());
    let h_field_base = B::get_root_of_unity(h_field_size.trailing_zeros());

    // / |L| >= 3*k_field_size - 3. For the rest of our code, we need to use powers of 2, hence we multiply by 4.
    let l_field_size = 4 * max_degree;
    let l_field_base = B::from(B::get_root_of_unity(l_field_size.trailing_zeros()));

    let i_field = winter_math::get_power_series(i_field_base, i_field_size);
    let h_field = winter_math::get_power_series_with_offset(h_field_base, params.eta, h_field_size);

    let k_field =
        winter_math::get_power_series_with_offset(k_field_base, params.eta_k, k_field_size);

    println!(
        "i: {}    k: {}    h: {}   L: {}",
        i_field_size, k_field_size, h_field_size, l_field_size
    );

    // Prepare the FFT coefficients (twiddles).
    let inv_twiddles_k_elts = fft::get_inv_twiddles::<B>(k_field_size);
    let twiddles_l_elts = fft::get_twiddles::<B>(l_field_size);

    IndexDomains {
        i_field_base: i_field_base,
        k_field_base: k_field_base,
        h_field_base: h_field_base,
        l_field_base: l_field_base,
        i_field: i_field,
        k_field: k_field,
        h_field: h_field,
        l_field_len: l_field_size,
        inv_twiddles_k_elts: inv_twiddles_k_elts,
        twiddles_l_elts: twiddles_l_elts,
        eta: params.eta,
        eta_k: params.eta_k,
        phantom_e: PhantomData::<E>,
    }
}

// Same as build_basefield_index_domains but for a prime field of size 17
pub fn build_primefield_index_domains(
    params: IndexParams<SmallFieldElement17>,
) -> IndexDomains<SmallFieldElement17, SmallFieldElement17> {
    let num_input_variables = params.num_input_variables;
    let num_constraints = params.num_constraints;
    let num_non_zero = params.num_non_zero;

    // Validate inputs.
    let ntpow2 = { |x: usize| x > 1 && (x & (x - 1) == 0) };
    assert!(
        ntpow2(num_input_variables),
        "num_input_variables {} must be nontriv power of two",
        num_input_variables
    );
    assert!(
        ntpow2(num_constraints),
        "num_constraints {} must be nontriv power of two",
        num_constraints
    );
    assert!(
        ntpow2(num_non_zero),
        "num_non_zero {} must be nontriv power of two",
        num_non_zero
    );

    // Need to encode a subset of H field: indices of inputs.
    let i_field_size = num_input_variables;

    // Need to enumerate each row or each column.
    // Input variables are across the columns, constraints are across the rows.
    let h_field_size = std::cmp::max(num_input_variables, num_constraints);

    // Need to enumerate each nonzero matrix entry.
    let k_field_size = num_non_zero;

    // Find elements in F which generate each subfield.
    let i_field_base = SmallFieldElement17::get_root_of_unity(i_field_size.trailing_zeros());
    let k_field_base = SmallFieldElement17::get_root_of_unity(k_field_size.trailing_zeros());
    let h_field_base = SmallFieldElement17::get_root_of_unity(h_field_size.trailing_zeros());

    // TODO: for p=17, cannot fit 4 here.  So use 2.  Buy a bigger prime for testing.
    //let l_field_size = 4 * num_non_zero; // this should actually be 3*k_field_size - 3 but will change later.
    let l_field_size = 2 * num_non_zero; // this should actually be 3*k_field_size - 3 but will change later.
    let l_field_base = SmallFieldElement17::get_root_of_unity(l_field_size.trailing_zeros());

    println!(
        "i: {}    k: {}    h: {}   L: {}",
        i_field_size, k_field_size, h_field_size, l_field_size
    );

    unsafe {
        // Find elements in F which generate each subfield.
        let i_field = SmallFieldElement17::get_power_series(i_field_base, i_field_size);
        let h_field = SmallFieldElement17::get_power_series(h_field_base, h_field_size);
        let k_field = SmallFieldElement17::get_power_series(k_field_base, k_field_size);

        // Prepare the FFT coefficients (twiddles).
        let inv_twiddles_k_elts = winter_math::fft::get_inv_twiddles(k_field_size);
        let twiddles_l_elts = winter_math::fft::get_twiddles(l_field_size);

        IndexDomains {
            i_field_base: i_field_base,
            k_field_base: k_field_base,
            h_field_base: h_field_base,
            l_field_base: l_field_base,
            i_field: i_field,
            k_field: k_field,
            h_field: h_field,
            l_field_len: l_field_size,
            inv_twiddles_k_elts: inv_twiddles_k_elts,
            twiddles_l_elts: twiddles_l_elts,
            eta: params.eta,
            eta_k: params.eta_k,
            phantom_e: PhantomData::<SmallFieldElement17>,
        }
    }
}

// TODO Update the new function for Index to take an R1CS instance as input.

pub fn create_index_from_r1cs<B: StarkField, E: FieldElement<BaseField = B>>(
    params: IndexParams<B>,
    r1cs_instance: R1CS<B>,
) -> Index<B, E> {
    let domains = build_index_domains(params.clone());
    let indexed_a = IndexedMatrix::new(r1cs_instance.A, &domains);
    let indexed_b = IndexedMatrix::new(r1cs_instance.B, &domains);
    let indexed_c = IndexedMatrix::new(r1cs_instance.C, &domains);
    Index::new(params, indexed_a, indexed_b, indexed_c)
}

pub fn create_primefield_index_from_r1cs(
    params: IndexParams<SmallFieldElement17>,
    r1cs_instance: R1CS<SmallFieldElement17>,
) -> Index<SmallFieldElement17, SmallFieldElement17> {
    let domains = build_primefield_index_domains(params.clone());
    let indexed_a = IndexedMatrix::new(r1cs_instance.A, &domains);
    let indexed_b = IndexedMatrix::new(r1cs_instance.B, &domains);
    let indexed_c = IndexedMatrix::new(r1cs_instance.C, &domains);
    Index::new(params, indexed_a, indexed_b, indexed_c)
}

pub fn get_max_degree(
    num_input_variables: usize,
    _num_constraints: usize,
    num_non_zero: usize,
) -> usize {
    let max_whole = max(
        num_input_variables - 1,
        max(2 * num_non_zero - 3, num_non_zero - 2),
    ) + 1;
    max_whole.next_power_of_two()
}
