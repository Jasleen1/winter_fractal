use crate::errors::FractalVerifierError;
use fractal_prover::prover_channel::FractalContext;
use log::debug;

use fractal_indexer::snark_keys::*;
use fractal_proofs::{FieldElement, FractalProof, StarkField};


use winter_crypto::{ElementHasher, RandomCoin};
use winter_utils::{collections::Vec, Serializable};

use crate::{lincheck_verifier::verify_lincheck_proof, rowcheck_verifier::verify_rowcheck_proof};

pub fn verify_fractal_proof<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher<BaseField = B>,
>(
    verifier_key: VerifierKey<H, B>,
    proof: FractalProof<B, E, H>,
    pub_inputs_bytes: Vec<u8>,
) -> Result<(), FractalVerifierError> {
    let mut coin_seed = pub_inputs_bytes;
    let fractal_context = FractalContext{
        index_commitments: verifier_key,
    };
    fractal_context.write_into(&mut coin_seed);
    let mut public_coin = RandomCoin::<_, H>::new(&coin_seed);
    let expected_alpha: B = public_coin.draw().expect("failed to draw OOD point");
    
    
    debug!("Lincheck a indexes: {:?}", &proof.lincheck_a.products_sumcheck_proof.queried_positions);
    verify_lincheck_proof(&verifier_key, proof.lincheck_a, &mut public_coin, expected_alpha)?;
    debug!("Lincheck a verified");
    verify_lincheck_proof(&verifier_key, proof.lincheck_b, &mut public_coin, expected_alpha)?;
    debug!("Lincheck b verified");
    verify_lincheck_proof(&verifier_key, proof.lincheck_c, &mut public_coin, expected_alpha)?;
    debug!("Lincheck c verified");
    verify_rowcheck_proof(&verifier_key, proof.rowcheck_proof, &mut public_coin)?;
    debug!("Rowcheck verified");
    
    Ok(())
}


/// Returns an out-of-domain point drawn uniformly at random from the public coin. 
/// This is used in lincheck to get the alpha value.
/// Alpha should be selected from Field \ h_domain, and we know h_domain is a multiplicative
/// coset, so we just check that before returning.
pub fn get_lincheck_alpha<B: StarkField, E: FieldElement<BaseField = B>, H: ElementHasher<BaseField = B>>(coin: &mut RandomCoin<B, H>, h_domain_size: usize, eta: E) -> E {
    let mut drawn_pt: E = coin.draw().expect("failed to draw OOD point");
    let h_domain_size_exp = E::PositiveInteger::from(h_domain_size as u64);
    while  (drawn_pt.div(eta)).exp(h_domain_size_exp) == E::ONE {
        coin.reseed(H::hash(&drawn_pt.to_bytes()));
        drawn_pt = coin.draw().expect("failed to draw OOD point");
    }
    drawn_pt
}
