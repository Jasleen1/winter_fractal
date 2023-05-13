use crate::errors::LowDegreeVerifierError;

use fractal_proofs::{polynom, FieldElement, LowDegreeBatchProof};
use fractal_utils::channel::DefaultFractalVerifierChannel;
use fractal_utils::polynomial_utils::*;
use winter_crypto::{Digest, ElementHasher, MerkleTree, RandomCoin};
use winter_fri::FriVerifier;
use winter_math::StarkField;

/// Verifies that all the values that are decomitted in the LowDegreeBatchProof correspond
/// to polynomials with the specified maximum degrees
#[cfg_attr(feature = "flame_it", flame("low_degree_verifier"))]
pub fn verify_low_degree_batch_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    proof: &LowDegreeBatchProof<B, E, H>,
    max_degrees: Vec<usize>,
    public_coin: &mut RandomCoin<B, H>,
    num_queries: usize,
) -> Result<(), LowDegreeVerifierError> {
    let mut channel = DefaultFractalVerifierChannel::<E, H>::new(
        proof.fri_proof.clone(),
        proof.commitments.clone(),
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;

    //todo: need to be able to sample these throughout the protocol like for the batch verifier
    let mut alphas: Vec<E> = Vec::with_capacity(max_degrees.len());
    let mut betas: Vec<E> = Vec::with_capacity(max_degrees.len());
    for _ in 0..max_degrees.len() {
        alphas.push(public_coin.draw::<E>().unwrap());
        betas.push(public_coin.draw::<E>().unwrap());
    }

    // rederive the evaluation domain size the same way as in the FRI verifier
    let eval_domain_size = proof.options.blowup_factor() * (proof.fri_max_degree + 1);
    public_coin.reseed(proof.tree_root);
    //for root in proof.tree_roots.iter() {
    //    public_coin.reseed(*root);
    //}
    let queried_positions = public_coin
        .draw_integers(num_queries, eval_domain_size)
        .unwrap();

    flame::start("verify fri");
    let fri_verifier = FriVerifier::<B, E, DefaultFractalVerifierChannel<E, H>, H>::new(
        &mut channel,
        public_coin,
        proof.options.clone(),
        proof.fri_max_degree,
    )?;
    fri_verifier.verify(
        &mut channel,
        &proof.composed_queried_evaluations,
        &queried_positions,
    )?;
    flame::end("verify fri");

    // Verify that merkle leaves are correct
    flame::start("verify merkle leaves");
    for i in (0..queried_positions.len()).into_iter() {
        let evals_at_idx: Vec<E> = proof
            .all_unpadded_queried_evaluations
            .iter()
            .map(|poly_evals| poly_evals[i])
            .collect();
        if H::hash_elements(&evals_at_idx) != proof.tree_proof.leaves[i] {
            println!(
                "Hash_elements applied to input array elts {:?}",
                proof
                    .all_unpadded_queried_evaluations
                    .iter()
                    .map(|x| H::hash_elements(x).as_bytes())
                    .collect::<Vec<[u8; 32]>>()
            );
            println!("Leaves {:?}", proof.tree_proof.leaves);
            return Err(LowDegreeVerifierError::MerkleTreeErr);
        }
    }
    flame::end("verify merkle leaves");

    flame::start("verify merkle batch");
    MerkleTree::verify_batch(&proof.tree_root, &queried_positions, &proof.tree_proof)
        .map_err(|_e| LowDegreeVerifierError::MerkleTreeErr)?;
    flame::end("verify merkle batch");

    verify_lower_degree_batch::<B, E, H>(
        eval_domain_size,
        max_degrees,
        proof.fri_max_degree,
        &proof.all_unpadded_queried_evaluations,
        &proof.composed_queried_evaluations,
        queried_positions,
        alphas,
        betas,
    )?;
    Ok(())
}

#[cfg_attr(feature = "flame_it", flame("low_degree_verifier"))]
fn verify_lower_degree_batch<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    eval_domain_size: usize,
    original_degrees: Vec<usize>,
    fri_max_degree: usize,
    original_evals: &Vec<Vec<E>>,
    final_evals: &Vec<E>,
    positions: Vec<usize>,
    alphas: Vec<E>,
    betas: Vec<E>,
) -> Result<(), LowDegreeVerifierError> {
    let eval_domain_base = E::from(B::get_root_of_unity(eval_domain_size.trailing_zeros()));
    let eval_domain_pows = positions.iter().map(|&x| x as u64).collect::<Vec<u64>>();
    let eval_domain_elts = eval_domain_pows
        .iter()
        .map(|&x| eval_domain_base.exp(E::PositiveInteger::from(x)))
        .collect::<Vec<E>>();

    //todo: use length of queried positions here
    let mut reconstructed_evals = vec![E::ZERO; eval_domain_elts.len()];
    for pos in 0..original_degrees.len() {
        let comp_poly = get_randomized_complementary_poly::<E>(
            original_degrees[pos],
            fri_max_degree,
            alphas[pos],
            betas[pos],
        );
        let eval_domain_evals = polynom::eval_many(&comp_poly, &eval_domain_elts);
        for pos2 in 0..eval_domain_elts.len() {
            reconstructed_evals[pos2] += original_evals[pos][pos2] * eval_domain_evals[pos2];
        }
    }
    for (pos, _) in eval_domain_elts.iter().enumerate() {
        if reconstructed_evals[pos] != final_evals[pos] {
            println!("Position {}", pos);
            println!("reconstructed_evals = {:?}", reconstructed_evals);
            println!("Final evals = {:?}", final_evals[pos]);
            return Err(LowDegreeVerifierError::PaddingErr);
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::verify_low_degree_batch_proof;
    use fractal_proofs::fields::QuadExtension;
    use fractal_proofs::{FieldElement, SumcheckProof};
    use fractal_utils::channel::DefaultFractalProverChannel;
    use low_degree_prover::low_degree_batch_prover::LowDegreeBatchProver;
    use std::time::{SystemTime, UNIX_EPOCH};
    use winter_crypto::hashers::{Blake3_256, Rp64_256};
    use winter_crypto::{ElementHasher, RandomCoin};
    use winter_fri::{FriOptions, FriVerifier, ProverChannel};
    use winter_math::fields::f64::BaseElement;
    use winter_math::utils;
    use winter_math::StarkField;

    #[test]
    fn run_test_low_degree_proof() {
        test_low_degree_proof::<BaseElement, QuadExtension<BaseElement>, Blake3_256<BaseElement>>();
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
        let max_degree: usize = 63;
        let l_field_size: usize = 4 * max_degree.next_power_of_two();
        let l_field_base = B::get_root_of_unity(l_field_size.trailing_zeros());
        let evaluation_domain = utils::get_power_series(l_field_base, l_field_size);

        //TODO: RandomCoin needs to be able to sample from the extension field, but currently can't
        //let mut public_coin = RandomCoin::<B, H>::new(&[]);
        let mut public_coin = RandomCoin::<_, H>::new(&vec![]);
        let mut channel = DefaultFractalProverChannel::<B, E, H>::new(
            evaluation_domain.len(),
            num_queries,
            vec![],
        );

        // this block can be removed if it's causing problems. Currently it's here to document weird behaviour
        let queried_positions = channel.draw_query_positions();
        assert!(public_coin.draw::<E>().unwrap() != channel.draw_fri_alpha());
        channel.public_coin.reseed(H::hash(&[1, 2, 3, 4]));
        public_coin.reseed(H::hash(&[1, 2, 3, 4]));
        assert!(public_coin.draw::<E>().unwrap() == channel.draw_fri_alpha());

        let mut prover = LowDegreeBatchProver::<B, E, H>::new(
            &evaluation_domain,
            fri_options.clone(),
            max_degree,
        );

        let max_degrees: Vec<usize> = vec![14, 63, 29, 31];
        let mut polys: Vec<Vec<B>> = Vec::new();
        for degree in max_degrees.iter() {
            let poly = nonrand_poly(*degree);
            prover.add_polynomial(&poly, *degree, &mut channel);
            polys.push(poly);
        }

        //println!("queried positions: {:?}", &queried_positions);
        let queried_positions = vec![
            144, 79, 190, 228, 234, 31, 172, 50, 78, 253, 194, 44, 21, 134, 22, 140,
        ];
        let proof = prover.generate_proof(&mut channel);
        assert!(
            verify_low_degree_batch_proof(&proof, max_degrees, &mut public_coin, num_queries)
                .is_ok()
        );

        assert!(public_coin.draw::<E>().unwrap() == channel.draw_fri_alpha());

        let mut prover = LowDegreeBatchProver::<B, E, H>::new(
            &evaluation_domain,
            fri_options.clone(),
            max_degree,
        );
        let max_degrees2: Vec<usize> = vec![37, 41, 36, 9];
        let mut polys: Vec<Vec<B>> = Vec::new();
        for degree in max_degrees2.iter() {
            let poly = nonrand_poly(*degree);
            prover.add_polynomial(&poly, *degree, &mut channel);
            polys.push(poly);
        }

        let proof2 = prover.generate_proof(&mut channel);
        assert!(verify_low_degree_batch_proof(
            &proof2,
            max_degrees2,
            &mut public_coin,
            num_queries
        )
        .is_ok());

        assert!(public_coin.draw::<E>().unwrap() == channel.draw_fri_alpha());
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
