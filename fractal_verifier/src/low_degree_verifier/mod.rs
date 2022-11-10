use crate::{channel::DefaultFractalVerifierChannel, errors::LowDegreeVerifierError};

use fractal_proofs::{polynom, FieldElement, LowDegreeProof};
use fractal_utils::polynomial_utils::*;
use winter_crypto::{ElementHasher, RandomCoin};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

// public_coin is used similarly to a proving channel. Why is that?
pub fn verify_low_degree_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    proof: LowDegreeProof<B, E, H>,
    max_degree: usize,
    public_coin: &mut RandomCoin<B, H>,
    num_queries: usize,
) -> Result<(), LowDegreeVerifierError> {
    let mut channel = DefaultFractalVerifierChannel::<E, H>::new(
        proof.fri_proof,
        proof.commitments,
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;

    public_coin.reseed(proof.tree_root.clone());
    // rederive the evaluation domain size the same way as in the FRI verifier
    let eval_domain_size = proof.options.blowup_factor() * (proof.fri_max_degree + 1);
    let queried_positions = public_coin.draw_integers(num_queries, eval_domain_size).unwrap();

    let fri_verifier = FriVerifier::<B, E, DefaultFractalVerifierChannel<E, H>, H>::new(
        &mut channel,
        public_coin,
        proof.options.clone(),
        proof.fri_max_degree,
    )?;
    
    //todo, are the queried position ever checked?
    fri_verifier.verify(
        &mut channel,
        &proof.padded_queried_evaluations,
        &queried_positions,
    )?;
    if max_degree < proof.fri_max_degree {
        verify_lower_degree::<B, E, H>(
            proof.options.blowup_factor() * (proof.fri_max_degree + 1),
            max_degree,
            proof.fri_max_degree,
            proof.unpadded_queried_evaluations,
            proof.padded_queried_evaluations,
            queried_positions,
        )?;
    }
    Ok(())
}

fn verify_lower_degree<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    eval_domain_size: usize,
    original_degree: usize,
    fri_max_degree: usize,
    original_evals: Vec<E>,
    final_evals: Vec<E>,
    positions: Vec<usize>,
) -> Result<(), LowDegreeVerifierError> {
    let comp_poly = get_complementary_poly::<E>(original_degree, fri_max_degree);
    let eval_domain_base = E::from(B::get_root_of_unity(eval_domain_size.trailing_zeros()));
    let eval_domain_pows = positions.iter().map(|&x| x as u64).collect::<Vec<u64>>();
    let eval_domain_elts = eval_domain_pows
        .iter()
        .map(|&x| eval_domain_base.exp(E::PositiveInteger::from(x)))
        .collect::<Vec<E>>();
    let eval_domain_evals = polynom::eval_many(&comp_poly, &eval_domain_elts);
    for (pos, _) in eval_domain_elts.iter().enumerate() {
        if original_evals[pos].mul(eval_domain_evals[pos]) != final_evals[pos] {
            println!("Position {}", pos);
            println!("Original_evals = {:?}", original_evals);
            println!("Domain elt = {:?}", eval_domain_elts[pos]);
            println!(
                "Mul = {:?}",
                original_evals[pos].mul(eval_domain_evals[pos])
            );
            println!("Final evals = {:?}", final_evals[pos]);
            return Err(LowDegreeVerifierError::PaddingErr); //::SmallPolyAdjustmentErr());
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::verify_low_degree_proof;
    use fractal_proofs::{FieldElement, SumcheckProof};
    use fractal_prover::channel::DefaultFractalProverChannel;
    use fractal_prover::low_degree_prover::LowDegreeProver;
    use std::time::{SystemTime, UNIX_EPOCH};
    use winter_crypto::hashers::Rp64_256;
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{DefaultProverChannel, DefaultVerifierChannel, FriOptions, FriVerifier};
    use winter_math::fields::f64::BaseElement;
    use winter_math::utils;
    use winter_math::StarkField;

    #[test]
    fn run_test_low_degree_proof() {
        test_low_degree_proof::<BaseElement, BaseElement, Rp64_256>();
    }

    fn test_low_degree_proof<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher<BaseField = B>,
    >() {
        let lde_blowup = 4;
        let num_queries = 16;
        let fri_options = FriOptions::new(lde_blowup, 4, 32);
        let max_degree = 63;
        let poly = nonrand_poly(max_degree);
        let l_field_size: usize = 4 * max_degree.next_power_of_two();
        let l_field_base = B::get_root_of_unity(l_field_size.trailing_zeros());
        let evaluation_domain = utils::get_power_series(l_field_base, l_field_size);
        let pub_input_bytes = vec![0u8];
        let mut public_coin = RandomCoin::<B, H>::new(&pub_input_bytes.clone());

        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            evaluation_domain.len(),
            num_queries,
            pub_input_bytes.clone(),
        );
        //let initial_queries = channel.draw_query_positions();
        let prover = LowDegreeProver::<B, E, H>::from_polynomial(
            &poly,
            &evaluation_domain,
            max_degree,
            fri_options.clone(),
        );
        let proof = prover.generate_proof(&mut channel);
        //assert!(verify_low_degree_proof(proof, 63, &mut public_coin, num_queries).is_ok());
        verify_low_degree_proof(proof, 63, &mut public_coin, num_queries).unwrap();

        let max_degree2 = 17;
        let poly2 = nonrand_poly(max_degree2);
        let prover = LowDegreeProver::<B, E, H>::from_polynomial(
            &poly2,
            &evaluation_domain,
            max_degree2,
            fri_options.clone(),
        );
        let proof2 = prover.generate_proof(&mut channel);
        assert!(verify_low_degree_proof(proof2, 17, &mut public_coin, num_queries).is_ok());
    }

    // a random-ish polynomial that isn't actually random at all. Instead, it uses the system clock since that doesn't require a new crate import
    fn nonrand_poly<B: StarkField>(degree: usize) -> Vec<B> {
        let mut out: Vec<B> = Vec::new();
        for _ in 0..degree + 1 {
            let nanos = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .subsec_nanos();
            out.push(B::from(nanos as u128));
        }
        out
    }
}
