use crate::errors::LincheckVerifierError;

use fractal_indexer::snark_keys::VerifierKey;
use fractal_proofs::{FieldElement, LincheckProof};
use fractal_sumcheck::{sumcheck_verifier::verify_sumcheck_proof};

use winter_crypto::{ElementHasher};
use winter_math::StarkField;

pub fn verify_lincheck_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: &VerifierKey<H, B>,
    proof: LincheckProof<B, E, H>,
    _expected_alpha: B,
) -> Result<(), LincheckVerifierError> {
    // LincheckProof::<B, E, H> {
    //     options: self.options.fri_options.clone(),
    //     num_evaluations: self.options.evaluation_domain.len(),
    //     alpha: self.alpha,
    //     beta,
    //     t_alpha_commitment,
    //     t_alpha_queried,
    //     products_sumcheck_proof,
    //     gamma,
    //     row_queried,
    //     col_queried,
    //     val_queried,
    //     matrix_sumcheck_proof,
    //     _e: PhantomData,
    // }
    // let mut public_coin_seed = Vec::new();
    // proof.write_into(&mut public_coin_seed);
    // let _public_coin: RandomCoin<B, H> = RandomCoin::new(&public_coin_seed);

    let _alpha = proof.alpha;
    println!("verifier alpha: {}", &_alpha);
    let _t_alpha_commitment = proof.t_alpha_commitment;
    let _t_alpha_queried = proof.t_alpha_queried;
    
    let products_sumcheck_proof = proof.products_sumcheck_proof;
    println!("Lincheck verifier indexes: {:?}", &products_sumcheck_proof.queried_positions);

    let h_field_size = std::cmp::max(verifier_key.params.num_input_variables, verifier_key.params.num_constraints);
    let g_degree = h_field_size - 2;
    let e_degree = h_field_size - 1;
    verify_sumcheck_proof(products_sumcheck_proof, g_degree, e_degree)
    .map_err(|err| LincheckVerifierError::UnsoundProduct(err))?;

    println!("Verified sumcheck for product");
    let _row_queried = proof.row_queried;
    let _col_queried = proof.col_queried;
    let _val_queried = proof.val_queried;

    let matrix_sumcheck_proof = proof.matrix_sumcheck_proof;
    let k_field_size = verifier_key.params.num_non_zero;
    let g_degree = k_field_size - 2;
    let e_degree = 2 * k_field_size - 3;
    verify_sumcheck_proof(matrix_sumcheck_proof, g_degree, e_degree)
    .map_err(|err| LincheckVerifierError::UnsoundMatrix(err))?;
    // Need to do the checking of beta and channel passing etc.
    // Also need to make sure that the queried evals are dealt with

    Ok(())
}
