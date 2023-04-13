use crate::accumulator_verifier::AccumulatorVerifier;
use crate::errors::{LincheckVerifierError, SumcheckVerifierError};

use crate::sumcheck_verifier::{verify_layered_sumcheck_proof, verify_sumcheck_proof};
use fractal_indexer::indexed_matrix::compute_derivative_xx;
use fractal_indexer::snark_keys::{ProverKey, VerifierKey};
use fractal_proofs::{
    compute_derivative_on_single_val, FieldElement, LayeredLincheckProof, LayeredSumcheckProof,
    LincheckProof, QueriedPositions,
};
use log::debug;

use winter_crypto::{ElementHasher, RandomCoin};
use winter_math::StarkField;

pub fn verify_lincheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, E, H>,
    proof: LincheckProof<B, E, H>,
    _expected_alpha: B,
    public_coin: &mut RandomCoin<B, H>,
    num_queries: usize,
) -> Result<(), LincheckVerifierError> {
    let _alpha = proof.alpha;
    println!(
        "Expected alpha vs sent alpha: {}",
        _expected_alpha == _alpha
    );
    debug!("verifier alpha: {}", &_alpha);
    let _t_alpha_commitment = proof.t_alpha_commitment;
    let _t_alpha_queried = proof.t_alpha_queried;

    let products_sumcheck_proof = proof.products_sumcheck_proof;
    debug!(
        "Lincheck verifier indexes: {:?}",
        &products_sumcheck_proof.queried_positions
    );

    let h_field_size = std::cmp::max(
        verifier_key.params.num_input_variables,
        verifier_key.params.num_constraints,
    );
    let g_degree = h_field_size - 2;
    let e_degree = h_field_size - 1;
    verify_sumcheck_proof(
        products_sumcheck_proof,
        g_degree,
        e_degree,
        public_coin,
        num_queries,
    )
    .map_err(|err| LincheckVerifierError::UnsoundProduct(err))?;

    debug!("Verified sumcheck for product");
    let _row_queried = proof.row_queried;
    let _col_queried = proof.col_queried;
    let _val_queried = proof.val_queried;

    //TODO: USE BETA
    let beta: B =
        FieldElement::as_base_elements(&[public_coin.draw::<E>().expect("failed to draw beta")])[0];

    let matrix_sumcheck_proof = proof.matrix_sumcheck_proof;
    let k_field_size = verifier_key.params.num_non_zero;
    let g_degree = k_field_size - 2;
    let e_degree = 2 * k_field_size - 3;
    verify_sumcheck_proof(
        matrix_sumcheck_proof,
        g_degree,
        e_degree,
        public_coin,
        num_queries,
    )
    .map_err(|err| LincheckVerifierError::UnsoundMatrix(err))?;
    // Need to do the checking of beta and channel passing etc.
    // Also need to make sure that the queried evals are dealt with

    Ok(())
}

pub(crate) fn verify_layered_lincheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
    verifier_key: &VerifierKey<B, E, H>,
    queried_positions: &Vec<usize>,
    proof: &LayeredLincheckProof<B, E>,
) -> Result<(), LincheckVerifierError> {
    let eta = verifier_key.params.eta;
    let h_size_u64: u64 = verifier_key.params.num_input_variables.try_into().unwrap();
    let k_size_u64: u64 = verifier_key.params.num_non_zero.try_into().unwrap();
    let l_size_u64: u64 = (verifier_key.params.max_degree * 4).try_into().unwrap();
    let l_base_elt = E::from(B::get_root_of_unity(
        l_size_u64.trailing_zeros().try_into().unwrap(),
    ));

    let v_h_alpha = compute_vanishing_poly::<B, E>(proof.alpha, h_size_u64, eta);
    let v_h_beta = compute_vanishing_poly::<B, E>(proof.beta, h_size_u64, eta);
    let eval_domain_size = verifier_key.params.max_degree * 4;
    let h_domain_size = verifier_key.params.num_input_variables;
    let k_domain_size = verifier_key.params.num_non_zero;
    accumulator_verifier.add_constraint(h_domain_size - 1);

    let mut product_sumcheck_g_decommits = Vec::<E>::new();
    let mut product_sumcheck_e_decommits = Vec::<E>::new();
    let mut product_sumcheck_numerator_decommits = Vec::<E>::new();
    let product_sumcheck_denominator_decommits = vec![E::ONE; queried_positions.len()];

    let mut matrix_sumcheck_g_decommits = Vec::<E>::new();
    let mut matrix_sumcheck_e_decommits = Vec::<E>::new();
    let mut matrix_sumcheck_numerator_decommits = Vec::<E>::new();
    let mut matrix_sumcheck_denominator_decommits = Vec::<E>::new();

    for i in 0..queried_positions.len() {
        let local_pow: u64 = queried_positions[i].try_into().unwrap();
        let current_x = l_base_elt.exp(E::PositiveInteger::from(local_pow));
        let u_alpha = compute_derivative(current_x, proof.alpha, h_size_u64);
        let f_1 = proof.f_mz_vals[i];
        let f_2 = proof.f_z_vals[i];
        let t_alpha = proof.t_alpha_vals[i];
        product_sumcheck_g_decommits.push(proof.product_sumcheck_vals[i].0);
        product_sumcheck_e_decommits.push(proof.product_sumcheck_vals[i].1);
        product_sumcheck_numerator_decommits.push((u_alpha * f_1) - (f_2 * t_alpha));

        matrix_sumcheck_g_decommits.push(proof.matrix_sumcheck_vals[i].0);
        matrix_sumcheck_e_decommits.push(proof.matrix_sumcheck_vals[i].1);
        matrix_sumcheck_numerator_decommits.push(proof.val_vals[i] * v_h_alpha * v_h_beta);
        matrix_sumcheck_denominator_decommits
            .push((proof.alpha - proof.col_vals[i]) * (proof.beta - proof.row_vals[i]));
    }

    let layered_product_sumcheck_proof = LayeredSumcheckProof {
        numerator_vals: product_sumcheck_numerator_decommits,
        denominator_vals: product_sumcheck_denominator_decommits,
        sumcheck_g_vals: product_sumcheck_g_decommits,
        sumcheck_e_vals: product_sumcheck_e_decommits,
    };

    verify_layered_sumcheck_proof::<B, E, H>(
        queried_positions,
        layered_product_sumcheck_proof,
        eval_domain_size,
        h_domain_size,
        B::ONE,
        eta,
        E::ZERO,
    )?;

    accumulator_verifier.add_constraint(h_domain_size - 2);
    accumulator_verifier.add_constraint(h_domain_size - 1);

    let layered_matrix_sumcheck_proof = LayeredSumcheckProof {
        numerator_vals: matrix_sumcheck_numerator_decommits,
        denominator_vals: matrix_sumcheck_denominator_decommits,
        sumcheck_g_vals: matrix_sumcheck_g_decommits,
        sumcheck_e_vals: matrix_sumcheck_e_decommits,
    };

    verify_layered_sumcheck_proof::<B, E, H>(
        queried_positions,
        layered_matrix_sumcheck_proof,
        eval_domain_size,
        k_domain_size,
        B::ONE,
        verifier_key.params.eta_k,
        proof.gamma,
    )?;

    accumulator_verifier.add_constraint(k_domain_size - 2);
    accumulator_verifier.add_constraint(2 * k_domain_size - 3);

    Ok(())
}

pub fn add_rational_sumcheck_verification<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    queried_positions: &Vec<usize>,
    numerator_decommits: Vec<E>,
    denominator_decommits: Vec<E>,
    g_decommits: Vec<E>,
    e_decommits: Vec<E>,
    eval_domain_size: usize,
    summing_domain_size: usize,
    eval_domain_offset: B,
    summing_domain_offset: B,
    gamma: E,
) -> Result<(), SumcheckVerifierError> {
    let summing_domain_size_u64: u64 = summing_domain_size.try_into().unwrap();
    let summing_domain_size_field = E::from(summing_domain_size_u64);
    let l_field_base = E::from(B::get_root_of_unity(
        eval_domain_size.trailing_zeros().try_into().unwrap(),
    ));
    let eta = summing_domain_offset;
    for i in 0..numerator_decommits.len() {
        let position_u64: u64 = queried_positions[i].try_into().unwrap();
        let x_val =
            l_field_base.exp(E::PositiveInteger::from(position_u64)) * E::from(eval_domain_offset);
        let denom_val = compute_vanishing_poly::<B, E>(x_val, summing_domain_size_u64, eta);
        let lhs = ((((x_val * g_decommits[i]) + (gamma / summing_domain_size_field))
            * denominator_decommits[i])
            - numerator_decommits[i])
            / denom_val;
        if lhs != e_decommits[i] {
            println!("lhs = {:?}, e = {:?}", lhs, e_decommits[i]);
            return Err(SumcheckVerifierError::ConsistentValuesErr(i));
        }
    }
    Ok(())
}

// should verify s was computed correctly and pass along the correct degree constraint
// just needs evals at queried positions?
pub fn add_lincheck_verification<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
    verifier_key: &VerifierKey<B, E, H>,
    decommit: Vec<Vec<E>>,
    row_idx: usize,
    col_idx: usize,
    val_idx: usize,
    f_z_idx: usize,
    f_mz_idx: usize,
    t_alpha_idx: usize,
    product_sumcheck_idxs: (usize, usize),
    matrix_sumcheck_idxs: (usize, usize),
    queried_positions: Vec<usize>,
    alpha: E,
    beta: E,
    gamma: E,
) -> Result<(), LincheckVerifierError> {
    let eta = verifier_key.params.eta;
    let h_size_u64: u64 = verifier_key.params.num_input_variables.try_into().unwrap();
    let l_size_u64: u64 = (verifier_key.params.max_degree * 4).try_into().unwrap();
    let l_base_elt = E::from(B::get_root_of_unity(
        l_size_u64.trailing_zeros().try_into().unwrap(),
    ));

    let v_h_alpha = compute_vanishing_poly::<B, E>(alpha, h_size_u64, eta);
    let v_h_beta = compute_vanishing_poly::<B, E>(beta, h_size_u64, eta);

    let eval_domain_size = verifier_key.params.max_degree * 4;
    let h_domain_size = verifier_key.params.num_input_variables;
    let k_domain_size = verifier_key.params.num_non_zero;

    accumulator_verifier.add_constraint(h_domain_size - 1);

    let mut product_sumcheck_g_decommits = Vec::<E>::new();
    let mut product_sumcheck_e_decommits = Vec::<E>::new();
    let mut product_sumcheck_numerator_decommits = Vec::<E>::new();
    let product_sumcheck_denominator_decommits = vec![E::ONE; queried_positions.len()];
    let mut matrix_sumcheck_g_decommits = Vec::<E>::new();
    let mut matrix_sumcheck_e_decommits = Vec::<E>::new();
    let mut matrix_sumcheck_numerator_decommits = Vec::<E>::new();
    let mut matrix_sumcheck_denominator_decommits = Vec::<E>::new();

    for i in 0..decommit.len() {
        product_sumcheck_g_decommits.push(decommit[i][product_sumcheck_idxs.0]);
        product_sumcheck_e_decommits.push(decommit[i][product_sumcheck_idxs.1]);

        let local_pow: u64 = queried_positions[i].try_into().unwrap();
        let current_x = l_base_elt.exp(E::PositiveInteger::from(local_pow));
        let u_alpha = compute_derivative(current_x, alpha, h_size_u64);

        let f_1 = decommit[i][f_mz_idx];
        let f_2 = decommit[i][f_z_idx];
        let t_alpha = decommit[i][t_alpha_idx];
        product_sumcheck_numerator_decommits.push((u_alpha * f_1) - (f_2 * t_alpha));

        matrix_sumcheck_g_decommits.push(decommit[i][matrix_sumcheck_idxs.0]);
        matrix_sumcheck_e_decommits.push(decommit[i][matrix_sumcheck_idxs.1]);
        matrix_sumcheck_numerator_decommits.push(decommit[i][val_idx] * v_h_alpha * v_h_beta);
        matrix_sumcheck_denominator_decommits
            .push((alpha - decommit[i][col_idx]) * (beta - decommit[i][row_idx]));
    }

    add_rational_sumcheck_verification::<B, E, H>(
        &queried_positions,
        product_sumcheck_numerator_decommits,
        product_sumcheck_denominator_decommits,
        product_sumcheck_g_decommits,
        product_sumcheck_e_decommits,
        eval_domain_size,
        h_domain_size,
        B::ONE,
        verifier_key.params.eta,
        E::ZERO,
    )?;

    accumulator_verifier.add_constraint(h_domain_size - 2);
    accumulator_verifier.add_constraint(h_domain_size - 1);

    println!("Checked the first sumcheck");

    add_rational_sumcheck_verification::<B, E, H>(
        &queried_positions,
        matrix_sumcheck_numerator_decommits,
        matrix_sumcheck_denominator_decommits,
        matrix_sumcheck_g_decommits,
        matrix_sumcheck_e_decommits,
        eval_domain_size,
        k_domain_size,
        B::ONE,
        verifier_key.params.eta_k,
        gamma,
    )?;

    accumulator_verifier.add_constraint(k_domain_size - 2);
    accumulator_verifier.add_constraint(2 * k_domain_size - 3);

    Ok(())
}

fn compute_vanishing_poly<B: StarkField, E: FieldElement<BaseField = B>>(
    element: E,
    size: u64,
    eta: B,
) -> E {
    let pow = E::PositiveInteger::from(size);
    element.exp(pow) - E::from(eta).exp(pow)
}

fn compute_derivative<B: StarkField, E: FieldElement<BaseField = B>>(
    x_elt: E,
    y_elt: E,
    dom_size: u64,
) -> E {
    if x_elt == y_elt {
        return compute_derivative_on_single_val(x_elt, dom_size.try_into().unwrap());
    }
    let power = E::PositiveInteger::from(dom_size);
    (x_elt.exp(power) - y_elt.exp(power)) / (x_elt - y_elt)
}

/// This function will change as we extend to also accumulate the lincheck parts
/// For now it takes in a vector of decommitted values and returns an aptly parsed decommitment.
/// It implicitly assumes that all the vectors of decommitted values are of the same length
pub(crate) fn prepare_lincheck_verifier_inputs<E: FieldElement>(
    decommits: Vec<Vec<Vec<E>>>,
) -> Vec<Vec<E>> {
    let mut return_vec = Vec::new();
    // To begin with, we use the corresponding decommitted values from the preprocessing
    let decommitted_preprocessing_row = decommits[0].clone();
    let decommitted_preprocessing_col = decommits[1].clone();
    let decommitted_preprocessing_val = decommits[2].clone();
    // Here, we first assume the first element of the input vec is the vector of f_mz evaluations
    let decommitted_fz_fmz = decommits[3].clone();
    // The second element is the evals of the t_alpha, poly_prod and sumcheck polynomial
    let decommitted_layer_2 = decommits[4].clone();
    // The third element is the evals of the matrix sumcheck protocol
    let decommitted_layer_3 = decommits[5].clone();

    for i in 0..decommitted_fz_fmz.len() {
        let mut latest_tuple = decommitted_preprocessing_row[i].clone();
        latest_tuple.push(decommitted_preprocessing_col[i][0]);
        latest_tuple.push(decommitted_preprocessing_val[i][0]);
        // Push the f_z poly
        latest_tuple.push(decommitted_fz_fmz[i][0]);
        // Push the f_mz poly
        latest_tuple.push(decommitted_fz_fmz[i][1]);
        // Push the t_alpha poly
        latest_tuple.push(decommitted_layer_2[i][0]);
        // Push the g poly from product sumcheck
        latest_tuple.push(decommitted_layer_2[i][1]);
        // Push the e poly from product sumcheck
        latest_tuple.push(decommitted_layer_2[i][2]);
        // Push the g poly from matrix sumcheck
        latest_tuple.push(decommitted_layer_3[i][0]);
        // Push the e poly from matrix sumcheck
        latest_tuple.push(decommitted_layer_3[i][1]);
        return_vec.push(latest_tuple);
    }

    return_vec
}

#[cfg(test)]
mod test {
    use crate::accumulator_verifier::AccumulatorVerifier;
    use crate::errors::TestingError;
    use crate::lincheck_verifier::{add_lincheck_verification, prepare_lincheck_verifier_inputs};
    use crate::rowcheck_verifier::add_rowcheck_verification;

    use super::verify_lincheck_proof;
    use fractal_examples2::gen_options::get_example_setup;
    use fractal_indexer::index::build_index_domains;
    use fractal_proofs::fields::QuadExtension;
    use fractal_proofs::{fft, polynom, FieldElement, SumcheckProof};
    use fractal_prover::accumulator::Accumulator;
    use fractal_prover::channel::DefaultFractalProverChannel;
    use fractal_prover::errors::ProverError;
    use fractal_prover::lincheck_prover::LincheckProver;
    use fractal_prover::prover::*;
    use fractal_prover::{FractalOptions, LayeredSubProver};
    use models::r1cs::Matrix;
    use std::ops::Add;
    use std::time::{SystemTime, UNIX_EPOCH};
    use winter_crypto::hashers::{Blake3_256, Rp64_256};
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{FriOptions, FriVerifier, ProverChannel};
    use winter_math::fields::f64::BaseElement;
    use winter_math::utils;
    use winter_math::StarkField;

    #[test]
    fn run_test_lincheck_proof() -> Result<(), TestingError> {
        test_lincheck_proof::<BaseElement, BaseElement, Rp64_256>()?;
        test_lincheck_proof::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>()?;
        Ok(())
    }

    fn test_lincheck_proof<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
    >() -> Result<(), TestingError> {
        // Here's an initial manual setup we won't be using, but could, if needed.
        /*let lde_blowup = 4;
        let num_queries = 16;
        let fri_options = FriOptions::new(lde_blowup, 4, 32);
        let max_degree: usize = 63;
        let l_field_size: usize = 4 * max_degree.next_power_of_two();
        let l_field_base = B::get_root_of_unity(l_field_size.trailing_zeros());
        let evaluation_domain = utils::get_power_series(l_field_base, l_field_size);
        let offset = B::ONE;
        let mut accumulator = Accumulator::<B,E,H>::new(evaluation_domain.len(), offset, evaluation_domain, num_queries, fri_options);

        let a = vec![0,1,2,3,4,5,6,7];
        let b = vec![2,2,2,2,2,2,2,2];
        let c = vec![0,2,4,6,8,10,12,14];

        let f_az_coeffs:Vec<> = a.iter().map(|x| B::from(*x as u128)).collect();
        let f_bz_coeffs:Vec<B> = b.iter().map(|x| B::from(*x as u128)).collect();
        let f_cz_coeffs:Vec<B> = c.iter().map(|x| B::from(*x as u128)).collect();
        */

        // SETUP TASKS

        // Let's first get the domains etc.
        let setup = get_example_setup::<B, E, H>();
        let (fractal_options, prover_key, verifier_key, wires) =
            (setup.0, setup.1, setup.2, setup.3);

        let setup_2 = get_example_setup::<B, E, H>();
        let (_, prover_key_2, verifier_key_2, wires_2) =
            (setup_2.0, setup_2.1, setup_2.2, setup_2.3);

        let evaluation_domain = fractal_options.evaluation_domain.clone();
        let eval_len = evaluation_domain.len();
        let h_domain = fractal_options.h_domain.clone();

        // PROVER TASKS
        // Actually generate the f_az, f_bz, f_cz polynomials
        // For this dummy example, we'll basically generate them randomly.
        // But remember! To be valid, f_cz must = f_az * f_bz.
        // random coefficients (todo, make it random)
        let num_coeffs = (h_domain.len()) as u128;

        let inv_twiddles_h = fft::get_inv_twiddles(wires.len());
        // 1. Generate lincheck coefficients for one matrix
        let mut z_coeffs = &mut wires.clone(); // evals
        fft::interpolate_poly_with_offset(&mut z_coeffs, &inv_twiddles_h, prover_key.params.eta); // coeffs
                                                                                                  // let matrix_a_index = prover_key.matrix_a_index;
                                                                                                  // Get the f_az coeffs
        let f_az_coeffs = &mut compute_matrix_mul_poly_coeffs::<B, E, H>(
            &prover_key.matrix_a_index.matrix,
            &wires.clone(),
            &inv_twiddles_h,
            prover_key.params.eta,
        )?;

        // Now that we have the f_Mz polynomials, we can commit to them, all in one go,
        // using the accumulator.
        let mut accumulator = Accumulator::<B, E, H>::new(
            eval_len,
            fractal_options.eta,
            evaluation_domain.clone(),
            fractal_options.num_queries,
            fractal_options.fri_options.clone(),
        );

        accumulator.add_unchecked_polynomial(z_coeffs.clone());
        accumulator.add_unchecked_polynomial(f_az_coeffs.clone());

        // Commit to the f_z and f_az polynomials before you move forward.
        let init_commit = accumulator.commit_layer()?;

        // Now the lincheck prover does its work.
        let mut lincheck_prover_a = LincheckProver::<B, E, H>::new(
            prover_key_2.matrix_a_index,
            f_az_coeffs.to_vec(),
            z_coeffs.to_vec(),
            &fractal_options,
        );
        let alpha = accumulator.draw_queries(Some(1))?[0];
        lincheck_prover_a
            .run_next_layer(alpha, &mut accumulator)
            .unwrap();
        // Now all the polynomials from the lincheck layer should be in the accumulator.
        // (spoiler: it's only three polynomials but we still need to commit it)
        let commit_layer_2 = accumulator.commit_layer()?;

        // // To show correctness, including of linking the two layers, query them at the same points
        // let decommit_fmz_polys = accumulator.decommit_layer_with_qeuries(1, layer_2_queries.clone())?;
        // let decommit_layer_2_lincheck = accumulator.decommit_layer_with_qeuries(2, layer_2_queries)?;

        let beta = accumulator.draw_queries(Some(1))?[0];

        assert!(verifier_key
            .matrix_b_commitments
            .row_poly_commitment
            .eq(&verifier_key_2.matrix_b_commitments.row_poly_commitment));

        lincheck_prover_a
            .run_next_layer(beta, &mut accumulator)
            .unwrap();

        let commit_layer_3 = accumulator.commit_layer()?;

        let layer_3_queries = accumulator.draw_query_positions()?;
        // To show correctness, including of linking the two layers, query them at the same points
        let preprocessed_values = prover_key.matrix_a_index.decommit_evals(&layer_3_queries)?;
        let decommit_layer_1_polys =
            accumulator.decommit_layer_with_qeuries(1, &layer_3_queries)?;
        let decommit_layer_2_polys =
            accumulator.decommit_layer_with_qeuries(2, &layer_3_queries)?;
        let decommit_layer_3_lincheck =
            accumulator.decommit_layer_with_qeuries(3, &layer_3_queries)?;

        let pp_0 = &preprocessed_values[0];
        let preprocessed_inputs_0 = &pp_0.0;
        // .iter()
        // .map(|val| vec![E::from(*val)])
        // .collect::<Vec<Vec<E>>>();
        let preprocessed_inputs_1 = &preprocessed_values[1].0;
        // .iter()
        // .map(|val| vec![E::from(*val)])
        // .collect::<Vec<Vec<E>>>();
        let preprocessed_inputs_2 = &preprocessed_values[2].0;
        // .iter()
        // .map(|val| vec![E::from(*val)])
        // .collect::<Vec<Vec<E>>>();

        // let preprocessed_proof_0 = pp_0.1;//prover_key_2.matrix_a_index.decommit_proof(layer_3_queries.clone(), 0)?;

        let lincheck_values = prepare_lincheck_verifier_inputs(vec![
            preprocessed_inputs_0.clone(),
            preprocessed_inputs_1.clone(),
            preprocessed_inputs_2.clone(),
            decommit_layer_1_polys.0.clone(),
            decommit_layer_2_polys.0.clone(),
            decommit_layer_3_lincheck.0.clone(),
        ]);

        // add some public input bytes VVV
        let fri_proof = accumulator.create_fri_proof()?;

        // VERIFIER TASKS
        // Instantiate the accumulator verifier to deal with all the merkle path verif.
        let mut accumulator_verifier = AccumulatorVerifier::<B, E, H>::new(
            eval_len,
            fractal_options.eta,
            evaluation_domain.clone(),
            fractal_options.num_queries,
            fractal_options.fri_options.clone(),
        );

        let query_indices = accumulator_verifier.get_query_indices(commit_layer_3);

        assert!(layer_3_queries == query_indices);

        // Check that the f_Mz decommitted values were appropriately sent by the prover
        println!("About to check accum for f_mz polynomials");
        assert!(accumulator_verifier.verify_layer_with_queries(
            init_commit,
            &query_indices,
            &decommit_layer_1_polys.0.clone(),
            &decommit_layer_1_polys.1
        ));
        // Check that the rowcheck layer decommitted values were appropriately sent.
        println!("About to check accum for everything inside the lincheck layer 1");
        assert!(accumulator_verifier.verify_layer_with_queries(
            commit_layer_2,
            &query_indices,
            &decommit_layer_2_polys.0.clone(),
            &decommit_layer_2_polys.1
        ));

        // Check that the rowcheck layer decommitted values were appropriately sent.
        println!("About to check accum for everything inside the lincheck layer 2");
        assert!(accumulator_verifier.verify_layer_with_queries(
            verifier_key.matrix_a_commitments.row_poly_commitment,
            &query_indices,
            &preprocessed_inputs_0,
            &pp_0.1,
        ));

        let gamma = lincheck_prover_a.retrieve_gamma(beta)?;
        add_lincheck_verification::<B, E, H>(
            &mut accumulator_verifier,
            &verifier_key,
            lincheck_values,
            0,
            1,
            2,
            3,
            4,
            5,
            (6, 7),
            (8, 9),
            query_indices,
            alpha,
            beta,
            gamma,
        )?;

        // Check correctness of FRI
        println!("About to check fri");
        assert!(accumulator_verifier.verify_fri_proof(commit_layer_3, fri_proof));
        /* Proof verification complete */

        // Also testing that the get_layer_commitment function is working as expected for the accum
        let first_layer_commit = accumulator.get_layer_commitment(1)?;
        let last_layer_commit = accumulator.get_layer_commitment(3)?;
        assert!(last_layer_commit == commit_layer_3);
        assert!(first_layer_commit == init_commit);

        Ok(())

        //IOP struct: vecs of commits, decommits, and a fri proof at the end
        // how does verifier know which proofs in which layers?
        // needs set of instructions: verify x constraint, move to next layer
        // as a first step, can you give it the full proof, then call functions in order?
    }

    // Multiply a matrix times a vector of evaluations, then interpolate a poly and return its coeffs.
    fn compute_matrix_mul_poly_coeffs<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
    >(
        matrix: &Matrix<B>,
        vec: &Vec<B>,
        inv_twiddles: &[B],
        eta: B,
    ) -> Result<Vec<B>, ProverError> {
        let mut product = matrix.dot(vec); // as evals
        fft::interpolate_poly_with_offset(&mut product, inv_twiddles, eta); // as coeffs
        Ok(product) // as coeffs
    }
}
