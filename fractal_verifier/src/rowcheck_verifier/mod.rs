use crate::{
    accumulator_verifier::AccumulatorVerifier, channel::DefaultFractalVerifierChannel,
    errors::RowcheckVerifierError, low_degree_verifier::verify_low_degree_proof,
};

use fractal_indexer::snark_keys::VerifierKey;
use fractal_proofs::{
    fft, get_complementary_poly, get_vanishing_poly, polynom, FieldElement, RowcheckProof, TryInto,
};

use log::debug;
use winter_crypto::{ElementHasher, MerkleTree, RandomCoin};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

pub fn verify_rowcheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<H, B>,
    proof: RowcheckProof<B, E, H>,
    public_coin: &mut RandomCoin<B, H>,
    initial_evals: Vec<Vec<B>>,
    num_queries: usize,
) -> Result<(), RowcheckVerifierError> {
    let indices = proof.s_proof.queried_positions.clone();
    let s_evals = proof.s_proof.unpadded_queried_evaluations.clone();

    let eval_domain_size = proof.options.blowup_factor() * verifier_key.params.max_degree;
    let h_domain_size = std::cmp::max(
        verifier_key.params.num_input_variables,
        verifier_key.params.num_constraints,
    );

    verify_low_degree_proof(
        proof.s_proof,
        verifier_key.params.max_degree - 1,
        public_coin,
        num_queries,
    )?;

    //verify_s_computation::<B, E, H>(eval_domain_size, h_domain_size, indices,
    //    E::from(verifier_key.params.eta), initial_evals, s_evals)?;

    Ok(())
}

// should verify s was computed correctly and pass along the correct degree constraint
// just needs evals at queried positions?
pub fn add_rowcheck_verification<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    accumulator_verifier: &mut AccumulatorVerifier<B, E, H>,
    verifier_key: &VerifierKey<H, B>,
    decommit: Vec<Vec<E>>,
    queried_positions: Vec<usize>, //Delete!
    f_az_idx: usize,
    f_bz_idx: usize,
    f_cz_idx: usize,
    s_idx: usize,
) -> Result<(), RowcheckVerifierError> {
    println!(
        "length of decommit: {}, {}",
        decommit.len(),
        decommit[1].len()
    );
    println!("length of queried_positions: {}", queried_positions.len());
    let initial_evals = vec![
        Vec::new(),
        decommit[f_az_idx].clone(),
        decommit[f_bz_idx].clone(),
        decommit[f_cz_idx].clone(),
    ];
    println!(
        "length of initial_evals: {}, {}",
        initial_evals.len(),
        initial_evals[0].len()
    );
    // todo: get this value from the same place consistently
    let h_domain_size = std::cmp::max(
        verifier_key.params.num_input_variables,
        verifier_key.params.num_constraints,
    );
    // The rowcheck is supposed to prove whether f_az * f_bz - f_cz = 0 on all of H.
    // Which means that the polynomial f_az * f_bz - f_cz must be divisible by the
    // vanishing polynomial for H.
    // Since the degree of f_az and f_bz is each |H| - 1, the degree of the polynomial
    // s = (f_az * f_bz - f_cz) / vanishing_H is upper bounded by |H| - 2.

    accumulator_verifier.add_constraint(h_domain_size - 2);

    let f_az_evals: Vec<E> = (0..queried_positions.len())
        .into_iter()
        .map(|i| decommit[i][f_az_idx])
        .collect();
    let f_bz_evals: Vec<E> = (0..queried_positions.len())
        .into_iter()
        .map(|i| decommit[i][f_bz_idx])
        .collect();
    let f_cz_evals: Vec<E> = (0..queried_positions.len())
        .into_iter()
        .map(|i| decommit[i][f_cz_idx])
        .collect();
    let s_evals: Vec<E> = (0..queried_positions.len())
        .into_iter()
        .map(|i| decommit[i][s_idx])
        .collect();

    verify_s_computation::<B, E, H>(
        accumulator_verifier.evaluation_domain_len,
        h_domain_size,
        queried_positions,
        E::from(verifier_key.params.eta),
        f_az_evals,
        f_bz_evals,
        f_cz_evals,
        s_evals,
    )?;

    Ok(())
}

fn verify_s_computation<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    eval_domain_size: usize,
    vanishing_domain_size: usize,
    positions: Vec<usize>,
    eta: E,
    f_az_evals: Vec<E>,
    f_bz_evals: Vec<E>,
    f_cz_evals: Vec<E>,
    s_evals: Vec<E>,
) -> Result<(), RowcheckVerifierError> {
    println!("s_evals: {:?}", s_evals);
    let eval_domain_base = E::from(B::get_root_of_unity(eval_domain_size.trailing_zeros()));
    let eval_domain_pows = positions.iter().map(|&x| x as u64).collect::<Vec<u64>>();
    let eval_domain_elts = eval_domain_pows
        .iter()
        .map(|&x| eval_domain_base.exp(E::PositiveInteger::from(x)))
        .collect::<Vec<E>>();
    let vanishing_poly = get_vanishing_poly(eta, vanishing_domain_size);

    let eval_domain_evals = polynom::eval_many(&vanishing_poly, &eval_domain_elts);
    /*let eval_domain_0 = eval_domain_evals[0];
    let vanishing_poly_0 = polynom::eval(&vanishing_poly.clone(), eval_domain_elts[0]);
    let initial_a = initial_evals[0][1];
    let initial_b = initial_evals[0][2];
    let initial_c = initial_evals[0][3];
    println!("A, b, c in verifier = {:?}, {:?}, {:?}", initial_a, initial_b, initial_c);
    println!("Vanishing poly 0 = {:?}", vanishing_poly_0);
    println!("Evals inside rowcheck verifier denom, s_coeffs = {:?}, {:?}", eval_domain_0, s_evals[0]);
    let s_0_computed = E::from(initial_a * initial_b - initial_c)/eval_domain_0;

    println!("S computed = {:?}", s_0_computed);
    println!("S array = {:?}", s_evals.clone());*/

    // todo: use a different reference for iterator
    for pos in 0..positions.len() {
        let s_val_computed =
            E::from((f_az_evals[pos] * f_bz_evals[pos]) - f_cz_evals[pos]) / eval_domain_evals[pos];
        if s_evals[pos] != s_val_computed {
            return Err(RowcheckVerifierError::ComputedValueMismatchErr(format!(
                "The computed polynomial s did not match the sent polynomial 
                at position {:?}, got {:?}, computed {:?}",
                pos, s_evals[pos], s_val_computed
            )));
        }
    }

    /*for (pos, ((init_vals, s_val), eval_domain_eval)) in initial_evals.iter().zip(s_evals).zip(eval_domain_evals).enumerate() {
        let f_az = init_vals[1];
        let f_bz = init_vals[2];
        let f_cz = init_vals[3];
        let s_val_computed = E::from((f_az * f_bz) - f_cz) / eval_domain_eval;
        if s_val != s_val_computed {
            return Err(RowcheckVerifierError::ComputedValueMismatchErr(
                format!("The computed polynomial s did not match the sent polynomial
                at position {:?}, got {:?}, computed {:?}", pos, s_val, s_val_computed)));
        }
    }*/
    Ok(())
}

#[cfg(test)]
mod test {
    use crate::accumulator_verifier::AccumulatorVerifier;
    use crate::errors::TestingError;
    use crate::rowcheck_verifier::add_rowcheck_verification;

    use super::verify_rowcheck_proof;
    use fractal_examples2::gen_options::get_example_setup;
    use fractal_indexer::index::build_index_domains;
    use fractal_proofs::fields::QuadExtension;
    use fractal_proofs::{polynom, FieldElement, SumcheckProof};
    use fractal_prover::accumulator::Accumulator;
    use fractal_prover::channel::DefaultFractalProverChannel;
    use fractal_prover::errors::ProverError;
    use fractal_prover::rowcheck_prover::RowcheckProver;
    use fractal_prover::{FractalOptions, LayeredProver};
    use std::ops::Add;
    use std::time::{SystemTime, UNIX_EPOCH};
    use winter_crypto::hashers::{Blake3_256, Rp64_256};
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{FriOptions, FriVerifier, ProverChannel};
    use winter_math::fields::f64::BaseElement;
    use winter_math::utils;
    use winter_math::StarkField;

    #[test]
    fn run_test_rowcheck_proof() -> Result<(), TestingError> {
        test_rowcheck_proof::<BaseElement, BaseElement, Rp64_256>()?;
        test_rowcheck_proof::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>()?;
        Ok(())
    }

    fn test_rowcheck_proof<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
    >() -> Result<(), TestingError> {
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

        let f_az_coeffs:Vec<B> = a.iter().map(|x| B::from(*x as u128)).collect();
        let f_bz_coeffs:Vec<B> = b.iter().map(|x| B::from(*x as u128)).collect();
        let f_cz_coeffs:Vec<B> = c.iter().map(|x| B::from(*x as u128)).collect();
        */
        let setup = get_example_setup::<B, E, H>();
        let (fractal_options, prover_key, verifier_key) = (setup.0, setup.1, setup.2);

        let evaluation_domain = fractal_options.evaluation_domain.clone();
        let eval_len = evaluation_domain.len();
        let h_domain = fractal_options.h_domain.clone();

        let mut accumulator = Accumulator::<B, E, H>::new(
            eval_len,
            fractal_options.eta,
            evaluation_domain.clone(),
            fractal_options.num_queries,
            fractal_options.fri_options.clone(),
        );

        // random coefficients (todo, make it random)
        let num_coeffs = (h_domain.len()) as u128;
        let a_coeffs: Vec<B> = (0..num_coeffs)
            .into_iter()
            .map(|i| B::from(2u128 * i + i * i + 4))
            .collect();
        let b_coeffs: Vec<B> = (0..num_coeffs)
            .into_iter()
            .map(|i| B::from(7u128 * i + i * i * i + 7))
            .collect();

        let a_evals_h = polynom::eval_many(&a_coeffs, &h_domain);
        let b_evals_h = polynom::eval_many(&b_coeffs, &h_domain);
        let c_evals_h: Vec<B> = (0..h_domain.len())
            .into_iter()
            .map(|i| *a_evals_h.get(i).unwrap() * *b_evals_h.get(i).unwrap())
            .collect();

        let c_coeffs = polynom::interpolate(&h_domain, &c_evals_h.to_vec(), true);
        println!("len(c_coeffs): {:?}", c_coeffs.len());

        /*let a_evals = polynom::eval_many(&a_coeffs, &evaluation_domain);
        let b_evals = polynom::eval_many(&b_coeffs, &evaluation_domain);
        let c_evals: Vec<B> = (0..eval_len).into_iter().map(|i| *a_evals.get(i).unwrap() * *b_evals.get(i).unwrap()).collect();

        let f_az_coeffs = polynom::interpolate(
            &evaluation_domain,
            &a_evals.to_vec(),
            true,
        );
        let f_bz_coeffs = polynom::interpolate(
            &evaluation_domain,
            &b_evals.to_vec(),
            true,
        );
        let f_cz_coeffs = polynom::interpolate(
            &evaluation_domain,
            &c_evals.to_vec(),
            true,
        );*/

        accumulator.add_unchecked_polynomial(a_coeffs.clone());
        accumulator.add_unchecked_polynomial(b_coeffs.clone());
        accumulator.add_unchecked_polynomial(c_coeffs.clone());
        let mut rowcheck_prover =
            RowcheckProver::<B, E, H>::new(a_coeffs, b_coeffs, c_coeffs, fractal_options.clone());
        let query = E::from(0u128);
        rowcheck_prover
            .run_next_layer(query, &mut accumulator)
            .unwrap();
        let commit = accumulator.commit_layer()?;
        let decommit = accumulator.decommit_layer(1)?;
        // add some public input bytes VVV
        let fri_proof = accumulator.create_fri_proof()?;

        let mut accumulator_verifier = AccumulatorVerifier::<B, E, H>::new(
            eval_len,
            fractal_options.eta,
            evaluation_domain.clone(),
            fractal_options.num_queries,
            fractal_options.fri_options.clone(),
        );

        assert!(accumulator_verifier.verify_layer(commit, decommit.0.clone(), decommit.1));

        //todo: do this in the verifier accumulator only
        let mut coin = RandomCoin::<B, H>::new(&vec![]);
        coin.reseed(commit);
        let queried_positions = coin
            .draw_integers(
                fractal_options.num_queries,
                fractal_options.evaluation_domain.len(),
            )
            .expect("failed to draw query position");

        println!("queried_positions: {:?}", &queried_positions);

        add_rowcheck_verification(
            &mut accumulator_verifier,
            &verifier_key,
            decommit.0,
            queried_positions,
            0,
            1,
            2,
            3,
        )?;

        assert!(accumulator_verifier.verify_fri_proof(commit, fri_proof));
        Ok(())

        //IOP struct: vecs of commits, decommits, and a fri proof at the end
        // how does verifier know which proofs in which layers?
        // needs set of instructions: verify x constraint, move to next layer
        // as a first step, can you give it the full proof, then call functions in order?
    }
}
