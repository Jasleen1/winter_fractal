use crate::{
    channel::DefaultFractalVerifierChannel, errors::RowcheckVerifierError,
    low_degree_verifier::verify_low_degree_proof,
};

use fractal_indexer::snark_keys::VerifierKey;
use fractal_proofs::{get_complementary_poly, polynom, FieldElement, RowcheckProof, TryInto, fft, get_vanishing_poly};

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
    let h_domain_size = std::cmp::max(verifier_key.params.num_input_variables, verifier_key.params.num_constraints);

    verify_low_degree_proof(
        proof.s_proof,
        verifier_key.params.max_degree - 1,
        public_coin,
        num_queries,
    )?;
    
    verify_s_computation::<B, E, H>(eval_domain_size, h_domain_size, indices, 
        E::from(verifier_key.params.eta), initial_evals, s_evals)?;

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
    initial_evals: Vec<Vec<B>>,
    s_evals: Vec<E>
) -> Result<(), RowcheckVerifierError> {
    let eval_domain_base = E::from(B::get_root_of_unity(eval_domain_size.trailing_zeros()));
    let eval_domain_pows = positions.iter().map(|&x| x as u64).collect::<Vec<u64>>();
    let eval_domain_elts = eval_domain_pows
        .iter()
        .map(|&x| eval_domain_base.exp(E::PositiveInteger::from(x)))
        .collect::<Vec<E>>();
    let vanishing_poly = get_vanishing_poly(eta, vanishing_domain_size);
    
    let eval_domain_evals = polynom::eval_many(&vanishing_poly, &eval_domain_elts);
    let eval_domain_0 = eval_domain_evals[0];
    let vanishing_poly_0 = polynom::eval(&vanishing_poly.clone(), eval_domain_elts[0]);
    let initial_a = initial_evals[0][1];
    let initial_b = initial_evals[0][2];
    let initial_c = initial_evals[0][3];
    println!("A, b, c in verifier = {:?}, {:?}, {:?}", initial_a, initial_b, initial_c);
    println!("Vanishing poly 0 = {:?}", vanishing_poly_0);
    println!("Evals inside rowcheck verifier denom, s_coeffs = {:?}, {:?}", eval_domain_0, s_evals[0]);
    let s_0_computed = E::from(initial_a * initial_b - initial_c)/eval_domain_0;

    println!("S computed = {:?}", s_0_computed);
    println!("S array = {:?}", s_evals.clone());

    for (pos, ((init_vals, s_val), eval_domain_eval)) in initial_evals.iter().zip(s_evals).zip(eval_domain_evals).enumerate() {
        let f_az = init_vals[1];
        let f_bz = init_vals[2];
        let f_cz = init_vals[3];
        let s_val_computed = E::from((f_az * f_bz) - f_cz) / eval_domain_eval;
        if s_val != s_val_computed {
            return Err(RowcheckVerifierError::ComputedValueMismatchErr(
                format!("The computed polynomial s did not match the sent polynomial 
                at position {:?}, got {:?}, computed {:?}", pos, s_val, s_val_computed)));
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::verify_rowcheck_proof;
    use fractal_indexer::index::build_index_domains;
    use fractal_proofs::fields::QuadExtension;
    use fractal_proofs::{FieldElement, SumcheckProof};
    use fractal_prover::{LayeredProver, FractalOptions};
    use fractal_prover::accumulator::Accumulator;
    use fractal_prover::channel::DefaultFractalProverChannel;
    use fractal_prover::rowcheck_prover::RowcheckProver;
    use fractal_examples2::gen_options::get_example_fractal_options;
    use std::time::{SystemTime, UNIX_EPOCH};
    use winter_crypto::hashers::{Rp64_256, Blake3_256};
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{FriOptions, FriVerifier, ProverChannel};
    use winter_math::fields::f64::BaseElement;
    use winter_math::utils;
    use winter_math::StarkField;

    #[test]
    fn run_test_rowcheck_proof() {
        test_rowcheck_proof::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>();
        test_rowcheck_proof::<BaseElement, BaseElement, Rp64_256>();
    }

    fn test_rowcheck_proof<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
    >() {

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
        let fractal_options = get_example_fractal_options::<B,E,H>();

        let evaluation_domain = fractal_options.evaluation_domain.clone();
        let mut accumulator = Accumulator::<B,E,H>::new(evaluation_domain.len(), fractal_options.eta, evaluation_domain, fractal_options.num_queries, fractal_options.fri_options.clone());


        let a_evals: Vec<B> = (0..31).into_iter().map(|_| B::from(2u128)).collect();
        let b_evals: Vec<B> = (0..31).into_iter().map(|_| B::from(7u128)).collect();
        let c_evals: Vec<B> = (0..31).into_iter().map(|i| *a_evals.get(i).unwrap() * *b_evals.get(i).unwrap()).collect();

        let f_az_coeffs:Vec<B> = interpolate a_evals...
        
        let mut rowcheck_prover = RowcheckProver::<B,E,H>::new(
            f_az_coeffs,
            f_bz_coeffs,
            f_cz_coeffs,
            fractal_options
        );
        let query = E::from(0u128);
        rowcheck_prover.run_next_layer(query, &mut accumulator).unwrap();
        let commit = accumulator.commit_layer();
        let fri_proof = accumulator.create_fri_proof();

    }
}