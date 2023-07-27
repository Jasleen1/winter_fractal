use crate::errors::{LincheckVerifierError, SumcheckVerifierError};
use fractal_accumulator_verifier::accumulator_verifier::AccumulatorVerifier;

use crate::sumcheck_verifier::{verify_layered_sumcheck_proof, verify_sumcheck_proof};
use fractal_indexer::indexed_matrix::compute_derivative_xx;
use fractal_indexer::snark_keys::{ProverKey, VerifierKey};
use fractal_proofs::{
    compute_derivative_on_single_val, BatchedLayeredLincheckProof, FieldElement,
    LayeredLincheckProof, LayeredSumcheckProof, LincheckProof, QueriedPositions, TopLevelProof,
};
use fractal_utils::FractalOptions;
use log::debug;

use winter_crypto::{ElementHasher, RandomCoin};
use winter_math::StarkField;

pub fn verify_lincheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
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

#[cfg_attr(feature = "flame_it", flame("lincheck_verifier"))]
pub fn verify_layered_lincheck_proof_from_top<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: VerifierKey<B, H>,
    proof: TopLevelProof<B, E, H>,
    pub_inputs_bytes: Vec<u8>,
    options: FractalOptions<B>,
) -> Result<(), LincheckVerifierError> {
    let mut accumulator_verifier: AccumulatorVerifier<B, E, H> = AccumulatorVerifier::new(
        options.evaluation_domain.len(),
        options.eta,
        options.evaluation_domain.clone(),
        options.num_queries,
        options.fri_options.clone(),
        pub_inputs_bytes.clone(),
    );

    // draw queries using only the last iop layer commit and the public input.
    // this helps keep the rngs in sync, but proper chaining of layers needs to be checked elsewhere!
    println!("layer commitment count: {}", &proof.layer_commitments.len());
    let query_seed = proof.layer_commitments[1];
    let mut coin = RandomCoin::<B, H>::new(&pub_inputs_bytes);
    coin.reseed(query_seed);

    let query_indices = coin
        .draw_integers(options.num_queries, options.evaluation_domain.len())
        .expect("failed to draw query position");

    verify_decommitments(
        &verifier_key,
        &proof,
        &query_indices,
        &mut accumulator_verifier,
    )?;

    let lincheck_proof = parse_proofs_for_subroutines(&verifier_key, &proof, &pub_inputs_bytes);
    verify_layered_lincheck_proof(
        &mut accumulator_verifier,
        &verifier_key,
        &query_indices,
        &lincheck_proof,
        1,
    )?;

    accumulator_verifier.verify_fri_proof(
        proof.layer_commitments[1],
        &proof.low_degree_proof,
        &pub_inputs_bytes,
    )?;

    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("lincheck_verifier"))]
pub fn verify_decommitments<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
    proof: &TopLevelProof<B, E, H>,
    query_indices: &Vec<usize>,
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
) -> Result<(), LincheckVerifierError> {
    // Verify that the committed preprocessing was queried correctly
    accumulator_verifier.verify_layer_with_queries(
        verifier_key.commitment,
        query_indices,
        &proof.preprocessing_decommitment.0,
        &proof.preprocessing_decommitment.1,
    )?;

    // Verify that the committed initial polynomials were queried correcly
    accumulator_verifier.verify_layer_with_queries(
        proof.initial_commitment,
        query_indices,
        &proof.initial_decommitment.0,
        &proof.initial_decommitment.1,
    )?;

    // Verify that the committed layers were queried correctly
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[0],
        query_indices,
        &proof.layer_decommitments[0].0,
        &proof.layer_decommitments[0].1,
    )?;
    accumulator_verifier.verify_layer_with_queries(
        proof.layer_commitments[1],
        query_indices,
        &proof.layer_decommitments[1].0,
        &proof.layer_decommitments[1].1,
    )?;

    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("lincheck_verifier"))]
pub fn parse_proofs_for_subroutines<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<B, H>,
    proof: &TopLevelProof<B, E, H>,
    public_inputs_bytes: &Vec<u8>,
) -> BatchedLayeredLincheckProof<B, E> {
    // Matrix A preprocessing
    let col_a = extract_vec_e(&proof.preprocessing_decommitment.0, 0);
    let row_a = extract_vec_e(&proof.preprocessing_decommitment.0, 1);
    let val_a = extract_vec_e(&proof.preprocessing_decommitment.0, 2);

    // Matrix B preprocessing
    let col_b = extract_vec_e(&proof.preprocessing_decommitment.0, 3);
    let row_b = extract_vec_e(&proof.preprocessing_decommitment.0, 4);
    let val_b = extract_vec_e(&proof.preprocessing_decommitment.0, 5);

    // Matrix C preprocessing
    let col_c = extract_vec_e(&proof.preprocessing_decommitment.0, 6);
    let row_c = extract_vec_e(&proof.preprocessing_decommitment.0, 7);
    let val_c = extract_vec_e(&proof.preprocessing_decommitment.0, 8);

    // get values from the initial polynomials
    let f_z_vals = extract_vec_e(&proof.initial_decommitment.0, 0);
    let f_az_vals = extract_vec_e(&proof.initial_decommitment.0, 1);
    let f_bz_vals = extract_vec_e(&proof.initial_decommitment.0, 2);
    let f_cz_vals = extract_vec_e(&proof.initial_decommitment.0, 3);

    // get values from the first layer
    let t_alpha_vals = extract_vec_e(&proof.layer_decommitments[0].0, 1);
    let product_sumcheck_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[0].0, 2, 3);

    // get values from the second layer
    let matrix_sumcheck_vals = extract_sumcheck_vec_e(&proof.layer_decommitments[1].0, 0, 1);

    // Sample our own alpha and beta to check the prover
    // Sample our own alpha and beta to check the prover
    let mut coin = RandomCoin::<B, H>::new(&public_inputs_bytes);
    coin.reseed(verifier_key.commitment);
    let _: E = coin.draw().expect("failed to draw FRI alpha");
    coin.reseed(proof.initial_commitment);
    let alpha: E = coin.draw().expect("failed to draw FRI alpha");
    coin.reseed(proof.layer_commitments[0]);
    let beta: E = coin.draw().expect("failed to draw FRI alpha");

    // let mut coin = RandomCoin::<B, H>::new(&public_inputs_bytes);
    // coin.reseed(proof.initial_commitment);
    // let alpha: E = coin.draw().expect("failed to draw FRI alpha");

    // coin.reseed(proof.layer_commitments[0]);
    // let beta: E = coin.draw().expect("failed to draw FRI alpha");

    let gammas = &proof.unverified_misc;

    BatchedLayeredLincheckProof {
        row_vals: [row_a, row_b, row_c],
        col_vals: [col_a, col_b, col_c],
        val_vals: [val_a, val_b, val_c],
        f_z_vals: f_z_vals,
        f_mz_vals: [f_az_vals, f_bz_vals, f_cz_vals],
        t_alpha_vals: t_alpha_vals,
        product_sumcheck_vals: product_sumcheck_vals,
        matrix_sumcheck_vals: matrix_sumcheck_vals,
        alpha,
        beta,
        gamma: gammas[0],
    }
}

fn extract_vec_e<B: StarkField, E: FieldElement<BaseField = B>>(
    vec_of_decommits: &Vec<Vec<E>>,
    position: usize,
) -> Vec<E> {
    vec_of_decommits
        .iter()
        .map(|x| x[position])
        .collect::<Vec<E>>()
}

fn extract_sumcheck_vec_e<B: StarkField, E: FieldElement<BaseField = B>>(
    vec_of_decommits: &Vec<Vec<E>>,
    position_g: usize,
    position_e: usize,
) -> Vec<(E, E)> {
    vec_of_decommits
        .iter()
        .map(|x| (x[position_g], x[position_e]))
        .collect::<Vec<(E, E)>>()
}

#[cfg_attr(feature = "flame_it", flame("lincheck_verifier"))]
pub(crate) fn verify_layered_lincheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    //todo: these params are ordered inconsistently
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
    verifier_key: &VerifierKey<B, H>,
    queried_positions: &Vec<usize>,
    proof: &BatchedLayeredLincheckProof<B, E>,
    starting_layer: usize,
) -> Result<(), LincheckVerifierError> {
    let eta = verifier_key.params.eta;
    let h_size_u64: u64 = verifier_key.params.num_input_variables.try_into().unwrap();
    let k_size_u64: u64 = verifier_key.params.num_non_zero.try_into().unwrap();
    let l_size_u64: u64 = (verifier_key.params.max_degree * 4).try_into().unwrap();
    let l_base_elt = E::from(B::get_root_of_unity(
        l_size_u64.trailing_zeros().try_into().unwrap(),
    ));

    let mut coin = RandomCoin::<B, H>::new(&[0]);
    // println!("Alpha = {:?}", proof.alpha);
    coin.reseed(H::hash(&proof.alpha.to_bytes()));
    let etas = [
        coin.draw().expect("failed to draw FRI alpha"),
        coin.draw().expect("failed to draw FRI alpha"),
        coin.draw().expect("failed to draw FRI alpha"),
    ];

    let v_h_alpha = compute_vanishing_poly::<B, E>(proof.alpha, h_size_u64, eta);
    let v_h_beta = compute_vanishing_poly::<B, E>(proof.beta, h_size_u64, eta);
    let eval_domain_size = verifier_key.params.max_degree * 4;
    let h_domain_size = verifier_key.params.num_input_variables;
    let k_domain_size = verifier_key.params.num_non_zero;
    accumulator_verifier.add_constraint(h_domain_size - 1, starting_layer);

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
        let mut f_1 = [E::ZERO; 3];
        for matrix_id in 0..3 {
            f_1[matrix_id] = etas[matrix_id] * proof.f_mz_vals[matrix_id][i];
        }

        let f_2 = proof.f_z_vals[i];
        let t_alpha = proof.t_alpha_vals[i];

        product_sumcheck_g_decommits.push(proof.product_sumcheck_vals[i].0);
        product_sumcheck_e_decommits.push(proof.product_sumcheck_vals[i].1);

        let mut product_numerator_term = E::ZERO - (f_2 * t_alpha);
        for matrix_id in 0..3 {
            product_numerator_term = product_numerator_term + (u_alpha * f_1[matrix_id]);
        }

        product_sumcheck_numerator_decommits.push(product_numerator_term);

        matrix_sumcheck_g_decommits.push(proof.matrix_sumcheck_vals[i].0);
        matrix_sumcheck_e_decommits.push(proof.matrix_sumcheck_vals[i].1);

        let mut mat_numerator_term = E::ZERO;
        let mut mat_denom_term = E::ONE;
        for matrix_id in 0..3 {
            let mat_denom_other_two = (proof.beta - proof.row_vals[(matrix_id + 1) % 3][i])
                * (proof.alpha - proof.col_vals[(matrix_id + 1) % 3][i])
                * (proof.beta - proof.row_vals[(matrix_id + 2) % 3][i])
                * (proof.alpha - proof.col_vals[(matrix_id + 2) % 3][i]);
            mat_numerator_term = mat_numerator_term
                + (proof.val_vals[matrix_id][i] * mat_denom_other_two * etas[matrix_id]);
            mat_denom_term = mat_denom_term
                * (proof.alpha - proof.col_vals[matrix_id][i])
                * (proof.beta - proof.row_vals[matrix_id][i]);
        }

        matrix_sumcheck_numerator_decommits.push(mat_numerator_term * v_h_alpha * v_h_beta);
        matrix_sumcheck_denominator_decommits.push(mat_denom_term);
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
        starting_layer,
    )?;

    // todo: g and e degree_max should be arguments to the sumcheck
    accumulator_verifier.add_constraint(h_domain_size - 2, starting_layer);
    accumulator_verifier.add_constraint(h_domain_size - 1, starting_layer);

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
        starting_layer + 1,
    )?;

    accumulator_verifier.add_constraint(k_domain_size - 2, starting_layer + 1);
    accumulator_verifier.add_constraint(6 * k_domain_size - 5, starting_layer + 1);

    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("lincheck_verifier"))]
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
#[cfg_attr(feature = "flame_it", flame("lincheck_verifier"))]
pub fn add_lincheck_verification<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
    verifier_key: &VerifierKey<B, H>,
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
    starting_layer: usize,
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

    accumulator_verifier.add_constraint(h_domain_size - 1, starting_layer);

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

    accumulator_verifier.add_constraint(h_domain_size - 2, starting_layer);
    accumulator_verifier.add_constraint(h_domain_size - 1, starting_layer);

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

    accumulator_verifier.add_constraint(k_domain_size - 2, starting_layer + 1);
    accumulator_verifier.add_constraint(2 * k_domain_size - 3, starting_layer + 1);

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
/*pub(crate) fn prepare_lincheck_verifier_inputs<E: FieldElement>(
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
}*/

#[cfg(test)]
mod test {

    use crate::errors::TestingError;
    use crate::lincheck_verifier::{
        add_lincheck_verification, verify_layered_lincheck_proof_from_top,
    };
    use crate::rowcheck_verifier::add_rowcheck_verification;
    use fractal_accumulator_verifier::accumulator_verifier::AccumulatorVerifier;

    use super::verify_lincheck_proof;
    use fractal_accumulator::accumulator::Accumulator;
    use fractal_examples2::gen_options::get_example_setup;
    use fractal_indexer::index::build_index_domains;
    use fractal_proofs::fields::QuadExtension;
    use fractal_proofs::{fft, polynom, FieldElement, SumcheckProof};
    use fractal_prover::errors::ProverError;
    use fractal_prover::lincheck_prover::LincheckProver;
    use fractal_prover::LayeredSubProver;
    use fractal_prover::{prover::*, LayeredProver};
    use fractal_utils::channel::DefaultFractalProverChannel;
    use fractal_utils::FractalOptions;
    use models::r1cs::Matrix;
    use std::ops::Add;
    use std::time::{SystemTime, UNIX_EPOCH};
    use winter_crypto::hashers::{Blake3_256, Rp64_256};
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{FriOptions, FriVerifier, ProverChannel};
    use winter_math::fields::f64::BaseElement;
    // // use winter_math::utils;
    use winter_math::StarkField;

    #[cfg_attr(feature = "flame_it", flame)]
    #[test]
    fn run_test_lincheck_proof() -> Result<(), TestingError> {
        test_lincheck_proof::<BaseElement, BaseElement, Rp64_256>()?;
        test_lincheck_proof::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>()?;
        #[cfg(feature = "flame_it")]
        flame::dump_html(&mut std::fs::File::create("stats/flame-graph.html").unwrap()).unwrap();
        Ok(())
    }

    #[cfg_attr(feature = "flame_it", flame)]
    fn test_lincheck_proof<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
    >() -> Result<(), TestingError> {
        // SETUP TASKS

        // Let's first get the domains etc.
        let setup = get_example_setup::<B, E, H>();
        let (prover_options, fractal_options, prover_key, verifier_key, wires) =
            (setup.0, setup.1, setup.2, setup.3, setup.4);

        let setup_2 = get_example_setup::<B, E, H>();
        let (_, _, prover_key_2, verifier_key_2, wires_2) =
            (setup_2.0, setup_2.1, setup_2.2, setup_2.3, setup_2.4);

        let evaluation_domain = fractal_options.evaluation_domain.clone();
        let eval_len = evaluation_domain.len();
        let h_domain = fractal_options.h_domain.clone();
        let pub_inputs_bytes = vec![];

        // PROVER TASKS
        println!("starting prover tasks");
        let num_coeffs = (h_domain.len()) as u128;

        let inv_twiddles_h = fft::get_inv_twiddles(wires.len());
        // Generate lincheck coefficients for one matrix
        let mut z_coeffs = wires.clone(); // evals
        fft::interpolate_poly_with_offset(&mut z_coeffs, &inv_twiddles_h, prover_key.params.eta);
        let f_az_coeffs = compute_matrix_mul_poly_coeffs::<B, E, H>(
            &prover_key.matrix_a_index.matrix,
            &wires.clone(),
            &inv_twiddles_h,
            prover_key.params.eta,
        )?;

        // Now the lincheck prover does its work.
        let mut lincheck_prover_a = LincheckProver::<B, E, H>::new(
            prover_key_2.matrix_a_index,
            f_az_coeffs,
            z_coeffs,
            // &fractal_options,
        );

        //flame::start("generate proof");
        let proof = lincheck_prover_a
            .generate_proof(&Some(prover_key), pub_inputs_bytes.clone(), &prover_options)
            .unwrap();
        //flame::end("generate proof");
        println!("starting verifier tasks");
        verify_layered_lincheck_proof_from_top(
            verifier_key,
            proof,
            pub_inputs_bytes,
            fractal_options,
        )
        .unwrap();

        Ok(())

        /*

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
            .commitment
            .eq(&verifier_key_2.commitment));

        lincheck_prover_a
            .run_next_layer(beta, &mut accumulator)
            .unwrap();

        let commit_layer_3 = accumulator.commit_layer()?;

        let layer_3_queries = accumulator.draw_query_positions()?;
        // To show correctness, including of linking the two layers, query them at the same points
        let preprocessed_values = prover_key.decommit_evals(&layer_3_queries)?;
        let col_a = extract_vec_e(&preprocessed_values.0, 0);
        let row_a = extract_vec_e(&preprocessed_values.0, 1);
        let val_a = extract_vec_e(&preprocessed_values.0, 2);

        let decommit_layer_1_polys =
            accumulator.decommit_layer_with_queries(1, &layer_3_queries)?;
        let decommit_layer_2_polys =
            accumulator.decommit_layer_with_queries(2, &layer_3_queries)?;
        let decommit_layer_3_lincheck =
            accumulator.decommit_layer_with_queries(3, &layer_3_queries)?;

        let pp_0 = &preprocessed_values[0];
        let preprocessed_inputs_0 = &row_a;
        // .iter()
        // .map(|val| vec![E::from(*val)])
        // .collect::<Vec<Vec<E>>>();
        let preprocessed_inputs_1 = &col_a;
        // .iter()
        // .map(|val| vec![E::from(*val)])
        // .collect::<Vec<Vec<E>>>();
        let preprocessed_inputs_2 = &val_a;
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
            pub_inputs_bytes.clone()
        );

        let query_indices = accumulator_verifier.get_query_indices(commit_layer_3, pub_inputs_bytes.clone())?;

        assert!(layer_3_queries == query_indices);

        // Check that the f_Mz decommitted values were appropriately sent by the prover
        println!("About to check accum for f_mz polynomials");
        accumulator_verifier.verify_layer_with_queries(
            init_commit,
            &query_indices,
            &decommit_layer_1_polys.0.clone(),
            &decommit_layer_1_polys.1,
        )?;
        // Check that the rowcheck layer decommitted values were appropriately sent.
        println!("About to check accum for everything inside the lincheck layer 1");
        accumulator_verifier.verify_layer_with_queries(
            commit_layer_2,
            &query_indices,
            &decommit_layer_2_polys.0.clone(),
            &decommit_layer_2_polys.1,
        )?;

        // Check that the rowcheck layer decommitted values were appropriately sent.
        println!("About to check accum for everything inside the lincheck layer 2");
        accumulator_verifier.verify_layer_with_queries(
            verifier_key.matrix_a_commitments.row_poly_commitment,
            &query_indices,
            &preprocessed_inputs_0,
            &pp_0.1,
        )?;

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
            1,
        )?;

        // Check correctness of FRI
        println!("About to check fri");
        accumulator_verifier.verify_fri_proof(commit_layer_3, fri_proof, pub_inputs_bytes.clone())?;
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
        */
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

    fn extract_vec_e<B: StarkField, E: FieldElement<BaseField = B>>(
        vec_of_decommits: &Vec<Vec<E>>,
        position: usize,
    ) -> Vec<E> {
        vec_of_decommits
            .iter()
            .map(|x| x[position])
            .collect::<Vec<E>>()
    }
}
