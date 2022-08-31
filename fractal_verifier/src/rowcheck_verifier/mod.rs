use crate::errors::RowcheckVerifierError;

use fractal_indexer::snark_keys::VerifierKey;
use fractal_proofs::{FieldElement, RowcheckProof, get_complementary_poly, polynom, TryInto};

use winter_crypto::{ElementHasher, RandomCoin, MerkleTree};
use winter_fri::{DefaultVerifierChannel, FriVerifier};
use winter_math::StarkField;

pub fn verify_rowcheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<H, B>,
    proof: RowcheckProof<B, E, H>,
    // Change to include public seed
) -> Result<(), RowcheckVerifierError> {
    // let mut public_coin_seed = Vec::new();
    // proof.write_into(&mut public_coin_seed);
    let mut public_coin = RandomCoin::new(&[]);

    let mut channel = DefaultVerifierChannel::new(
        proof.s_proof,
        proof.s_commitments,
        proof.num_evaluations,
        proof.options.folding_factor(),
    )?;
    let s_queried_evals = proof.s_queried_evals;
    let s_original_evals = proof.s_original_evals;
    
    let s_original_proof = proof.s_original_proof;
    MerkleTree::verify_batch(&proof.s_eval_root, &proof.queried_positions.clone(), &s_original_proof).map_err(|err| RowcheckVerifierError::MerkleTreeErr(err))?;
    verify_lower_degree::<B, E, H>(4 * verifier_key.params.max_degree, verifier_key.params.num_input_variables - 1, verifier_key.params.max_degree, s_original_evals, s_queried_evals.clone(), proof.queried_positions.clone())?;
    

    let fri_verifier = FriVerifier::<B, E, DefaultVerifierChannel<E, H>, H>::new(
        &mut channel,
        &mut public_coin,
        proof.options.clone(),
        verifier_key.params.max_degree - 1,
    )?;
    //fri_verifier.verify(&mut channel, &s_queried_evals, &proof.queried_positions)
    fri_verifier.verify(&mut channel, &s_queried_evals, &proof.queried_positions).map_err(|err| RowcheckVerifierError::FriVerifierErr(err))
}


fn verify_lower_degree<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(eval_domain_size: usize, original_degree: usize, max_degree: usize, 
    original_evals: Vec<E>, final_evals: Vec<E>, positions: Vec<usize>) -> Result<(), RowcheckVerifierError> {
    let comp_poly = get_complementary_poly::<E>(original_degree, max_degree - 1);
    let eval_domain_base = E::from(B::get_root_of_unity(eval_domain_size.trailing_zeros()));
    let eval_domain_pows = positions.iter().map(|&x| {let z: u64 = x.try_into().unwrap(); z}).collect::<Vec<u64>>();
    let eval_domain_elts = eval_domain_pows.iter().map(|&x| eval_domain_base.exp(E::PositiveInteger::from(x))).collect::<Vec<E>>();
    let eval_domain_evals = polynom::eval_many(&comp_poly, &eval_domain_elts);
    for (pos, _) in eval_domain_elts.iter().enumerate() {
        if original_evals[pos].mul(eval_domain_evals[pos]) != final_evals[pos] {
            println!("Position {}", pos);
            println!("Original_evals = {:?}", original_evals);
            println!("Domain elt = {:?}", eval_domain_elts[pos]);
            println!("Mul = {:?}", original_evals[pos].mul(eval_domain_evals[pos]));
            println!("Final evals = {:?}", final_evals[pos]);
            return Err(RowcheckVerifierError::SmallPolyAdjustmentErr());
        }
    }
    Ok(())
}

