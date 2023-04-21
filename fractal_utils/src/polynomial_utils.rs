use crate::{errors::FractalUtilError, matrix_utils::*};
use fractal_math::{fft, FieldElement, StarkField};
use std::{convert::TryInto, marker::PhantomData};
use winter_crypto::{BatchMerkleProof, Digest, ElementHasher, MerkleTree};
use winter_fri::{DefaultProverChannel, FriOptions};
use winter_utils::batch_iter_mut;
// TODO: Add error checking and throwing
/**
 * This is equivalent to computing v_H(X) for a multiplicative coset
 * H = eta * H_0 for a multiplicative subgroup H_0 of order dom_size
 * Note that v_H(X) = X^dom_size - eta^dom_size. If eta = 1 we're
 * in a multiplicative subgroup itself.
 **/
#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn compute_vanishing_poly<E: FieldElement>(x: E, eta: E, dom_size: usize) -> E {
    let power_u64: u64 = dom_size.try_into().unwrap();
    let power = E::PositiveInteger::from(power_u64);
    x.exp(power) - eta.exp(power)
}

/// This function generates the vanshing polynomial coefficients for a multiplicative
/// subgroup of size dom_size and with multiplicative factor eta.
#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn get_vanishing_poly<E: FieldElement>(eta: E, dom_size: usize) -> Vec<E> {
    let mut vanishing_poly = vec![E::ZERO; dom_size];
    vanishing_poly.push(E::ONE);
    let dom_size_32: u32 = dom_size.try_into().unwrap();
    let eta_pow = E::PositiveInteger::from(dom_size_32);
    vanishing_poly[0] = eta.exp(eta_pow).neg();
    vanishing_poly
}

/**
 * Compute vanishing polynomial for a multiplicative subgroup. Same as above with
 * eta = ONE.
 **/
pub fn vanishing_poly_for_mult_subgroup<E: FieldElement>(x: E, dom_size: usize) -> E {
    compute_vanishing_poly(x, E::ONE, dom_size)
}

// The derivative is calculated as |H| x^{|H|-1}.
// This is equivalent to computing u_H(X, X) for a multiplicative coset H
// of order dom_size = |H|.
#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn compute_derivative_on_single_val<E: FieldElement>(x: E, dom_size: u128) -> E {
    let dom_size_coeff = E::from(dom_size);
    let power_u64: u64 = (dom_size - 1).try_into().unwrap();
    let power = E::PositiveInteger::from(power_u64);
    // println!("   deriv: {}    = h: {}   x: {:?}  ** h1: {:?}",
    //     dom_size_coeff * x.exp(power), dom_size_coeff, x, power);
    dom_size_coeff * x.exp(power)
}

// Represents a binomial, i.e. a polynomial in two variables X and Y.
// The (i, j)th element of this binomial is the coefficient of X^i * Y^j
pub type BivariatePoly<E> = Vec<Vec<E>>;

#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn compute_binomial_on_x<E: FieldElement>(bivariate: BivariatePoly<E>, x_val: E) -> Vec<E> {
    // Given a BivariatePoly, computes a monomial in Y, obtained by evaluating bivariate(x_val, Y)
    // represented by the output vector
    let bivariate_as_matrix = Matrix::new("binomial", bivariate).unwrap();
    let transposed_bivariate = bivariate_as_matrix.get_transpose("transposed_binomial");
    compute_binomial_on_y(transposed_bivariate.mat, x_val)
}

#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn compute_binomial_on_y<E: FieldElement>(bivariate: BivariatePoly<E>, y_val: E) -> Vec<E> {
    // Given a BivariatePoly, computes a monomial in X, obtained by evaluating bivariate(X, y_val)
    // represented by the output vector
    // Note that since bivariate[i][j] is the coefficient of X^i Y^j, technically,
    // bivariate[i] * (Y^j)_j is the coeffient of X^i
    // Hence, evaluating the polynomial bivariate[i] on y_val gives the coeff of X^i
    let mut x_coeffs = Vec::new();
    for i in 0..bivariate.len() {
        x_coeffs.push(fractal_math::polynom::eval(&bivariate[i], y_val));
    }
    x_coeffs
}

pub fn pad_with_zeroes<E: FieldElement>(poly: &mut Vec<E>, total_len: usize) {
    if total_len <= poly.len() {
        return;
    }
    let diff = total_len - poly.len();
    let mut zeroes = vec![E::ZERO; diff];
    poly.append(&mut zeroes);
}

pub fn pad_to_next_power_of_two<E: FieldElement>(poly: &mut Vec<E>) {
    let total_len = poly.len().next_power_of_two();
    let diff = total_len - poly.len();
    let mut zeroes = vec![E::ZERO; diff];
    poly.append(&mut zeroes);
}

pub fn get_to_degree_size<E: FieldElement>(poly: &mut Vec<E>) {
    let len_poly = poly.len();
    let mut count = len_poly - 1;
    while count > 0 && poly[count] == E::ZERO {
        poly.pop();
        count = count - 1;
    }
}

#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn get_complementary_poly<E: FieldElement>(
    current_degree: usize,
    desired_degree: usize,
) -> Vec<E> {
    assert!(desired_degree >= current_degree);
    let comp_deg = desired_degree - current_degree;
    let mut out_poly = vec![E::ZERO; comp_deg];
    out_poly.push(E::ONE);
    out_poly[0] = E::ONE;
    out_poly
}

#[cfg_attr(feature = "flame_it", flame("utils"))]
pub fn get_randomized_complementary_poly<E: FieldElement>(
    current_degree: usize,
    desired_degree: usize,
    alpha: E,
    beta: E,
) -> Vec<E> {
    assert!(desired_degree >= current_degree);
    let comp_deg = desired_degree - current_degree;
    let mut out_poly = vec![E::ZERO; comp_deg];
    out_poly.push(alpha);
    out_poly[0] = beta;
    out_poly
}

pub trait MultiPoly<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
>
{
    /// This function should take as input self and commit the values of the polynomials in question.
    /// It outputs the commitment.
    fn commit_polynomial_evaluations(&mut self) -> Result<(), FractalUtilError>;
    /// This function retrieves the commitment to the polynomials.
    fn get_commitment(&self) -> Result<&H::Digest, FractalUtilError>;
    /// This function retrieves the evaluations of the polynomials in question at the given
    /// index in the evaluation domain.
    fn get_values_at(&self, index: usize) -> Result<Vec<E>, FractalUtilError>;
    /// This function retrieves the evals of the polynomials at a set of evaluation points.
    fn batch_get_values_at(&self, indices: &Vec<usize>) -> Result<Vec<Vec<E>>, FractalUtilError>;
    /// This function takes as input an index of a point in the evaluation domain and
    /// outputs the evals committed at that point and a proof.
    fn get_values_and_proof_at(
        &self,
        index: usize,
    ) -> Result<(Vec<E>, Vec<H::Digest>), FractalUtilError>;
    /// This function takes as input the indices of multiple points in the evaluation domain and
    /// returns the evaluations of all the polynomials at these points, together with a batch merkle
    /// proof showing that this eval was done correctly.
    fn batch_get_values_and_proofs_at(
        &self,
        indices: &Vec<usize>,
    ) -> Result<(Vec<Vec<E>>, BatchMerkleProof<H>), FractalUtilError>;
    /// This function takes as input the value of the polynomials at a particular index and
    /// verifies it wrt to the commitment.
    /// Note how this function is stateless, so it can be efficiently used by the verifier.
    fn verify_values_and_proof_at(
        vals: Vec<B>,
        root: &H::Digest,
        proof: &[H::Digest],
        index: usize,
    ) -> Result<(), FractalUtilError>;
    /// This function takes as input the value of the polynomials at multiple indices and
    /// verifies them wrt to the commitment.
    fn batch_verify_values_and_proofs_at(
        vals: &Vec<Vec<E>>,
        root: &H::Digest,
        proof: &BatchMerkleProof<H>,
        indices: &Vec<usize>,
    ) -> Result<(), FractalUtilError>;

    fn batch_verify_transposed_values_and_proofs_at(
        vals: Vec<[E; 1]>,
        root: &<H>::Digest,
        proof: &BatchMerkleProof<H>,
        indices: &Vec<usize>,
    ) -> Result<(), FractalUtilError>;
}

pub struct MultiEval<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub evaluations: Vec<Vec<E>>,
    pub coefficients: Vec<Vec<E>>,
    pub committed_tree: Option<MerkleTree<H>>,
    _e: PhantomData<E>,
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > MultiEval<B, E, H>
{
    /// This function takes as input a set of polynomials in coefficient form,
    /// an evaluation domain and evaluates the polynomials at that domain.
    /// Note that coefficients is semantically of the form <poly_1, ..., poly_n>
    /// that is, each element of the vector coefficients is the vector of coefficients
    /// for one of the polynomials in question.
    pub fn new(
        coefficients_b: Vec<Vec<B>>,
        coefficients_e: Vec<Vec<E>>,
        evaluation_domain_len: usize,
        // TODO: offset is not used. Currently fine as the offset for eval_domain is ONE
        offset: B,
    ) -> Self {
        let eval_twiddles = fft::get_twiddles(evaluation_domain_len);

        let mut accumulated_evals = Vec::<Vec<E>>::new();
        for (_, poly) in coefficients_b.iter().enumerate() {
            accumulated_evals.push(
                eval_on_domain(poly, evaluation_domain_len, &eval_twiddles)
                    .into_iter()
                    .map(|i| E::from(i))
                    .collect(),
            );
        }

        for (_, poly) in coefficients_e.iter().enumerate() {
            accumulated_evals.push(eval_on_domain(poly, evaluation_domain_len, &eval_twiddles));
        }

        let mut coefficients = coefficients_e;
        for (_, poly) in coefficients_b.into_iter().enumerate() {
            coefficients.push(poly.into_iter().map(|i| E::from(i)).collect());
        }

        let evaluations = Self::zip_evals(accumulated_evals, evaluation_domain_len);
        let committed_tree: Option<MerkleTree<H>> = Option::None;
        Self {
            evaluations,
            coefficients,
            committed_tree,
            _e: PhantomData,
        }
    }

    // Todo: Bug. This function does not use zip_evals, and so probably pushes values incorrectly
    // luckily, doesn't seem to be used got anything right now
    /*pub fn add_polynomial(&mut self, coefficients: Vec<B>, evaluation_domain_len: usize) -> () {
        let eval_twiddles = fft::get_twiddles(evaluation_domain_len);
        let evaluations = eval_on_domain(&coefficients, evaluation_domain_len, &eval_twiddles);
        self.coefficients.push(
            coefficients
                .into_iter()
                .map(|i| E::from(i))
                .collect::<Vec<E>>(),
        );
        self.evaluations.push(
            evaluations
                .into_iter()
                .map(|i| E::from(i))
                .collect::<Vec<E>>(),
        );
        self.committed_tree = Option::None;
    }*/

    /// This is mostly a helper function to evaluate the polynomials on a domain of given length
    /// for which twiddles are already computed.
    /*pub fn eval_on_domain(
        coefficients: &[E],
        evaluation_domain_len: usize,
        eval_twiddles: &[B],
    ) -> Vec<E> {
        let mut eval = coefficients.clone();
        pad_with_zeroes(&mut eval, evaluation_domain_len);

        fft::evaluate_poly(&mut eval, eval_twiddles);

        eval.to_vec()
    }*/

    /// Helper function to zip the evaluations so that each element of the output is of the
    /// form [poly_1(e), ..., poly_n(e)] i.e. evaluations of all the polynomials are included
    /// in the same array.
    #[cfg_attr(feature = "flame_it", flame("utils"))]
    fn zip_evals(separate_evals: Vec<Vec<E>>, evaluation_domain_len: usize) -> Vec<Vec<E>> {
        let mut zipped_evals = vec![Vec::<E>::new(); evaluation_domain_len];
        for (_, eval) in separate_evals.iter().enumerate() {
            for (loc, &val) in eval.iter().enumerate() {
                zipped_evals[loc].push(val);
            }
        }
        zipped_evals
    }
}

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > MultiPoly<B, E, H> for MultiEval<B, E, H>
{
    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn commit_polynomial_evaluations(&mut self) -> Result<(), FractalUtilError> {
        // todo!()
        let eval_hashes = self
            .evaluations
            .iter()
            .map(|evals| H::hash_elements(evals))
            .collect::<Vec<_>>();
        let com_tree = MerkleTree::new(eval_hashes).map_err(|e| {
            FractalUtilError::MultiPolyErr(format!(
                "Got an error when committing to the evals: {e}"
            ))
        })?;
        self.committed_tree = Some(com_tree);
        Ok(())
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn get_commitment(&self) -> Result<&<H as winter_crypto::Hasher>::Digest, FractalUtilError> {
        match &self.committed_tree {
            Some(merkle_tree) => Ok(merkle_tree.root()),
            None => Err(FractalUtilError::MultiPolyErr(
                "The Merkle Tree in the poly commit is None.".to_string(),
            )),
        }
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn get_values_at(&self, index: usize) -> Result<Vec<E>, FractalUtilError> {
        Ok(self.evaluations[index].clone())
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn batch_get_values_at(&self, indices: &Vec<usize>) -> Result<Vec<Vec<E>>, FractalUtilError> {
        let mut output_vals = Vec::<Vec<E>>::new();
        for (_, &index) in indices.iter().enumerate() {
            output_vals.push(self.evaluations[index].clone());
        }
        Ok(output_vals)
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn get_values_and_proof_at(
        &self,
        index: usize,
    ) -> Result<(Vec<E>, Vec<<H>::Digest>), FractalUtilError> {
        let value = self.evaluations[index].clone();
        let proof = match &self.committed_tree {
            None => Err(FractalUtilError::MultiPolyErr(
                "Nothing committed yet!".to_string(),
            )),
            Some(tree) => tree.prove(index).map_err(|e| {
                FractalUtilError::MultiPolyErr(format!(
                    "Got an error when committing to the evals: {e}"
                ))
            }),
        }?;
        Ok((value, proof))
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn batch_get_values_and_proofs_at(
        &self,
        indices: &Vec<usize>,
    ) -> Result<(Vec<Vec<E>>, BatchMerkleProof<H>), FractalUtilError> {
        let values = self.batch_get_values_at(indices)?;
        let proof = match &self.committed_tree {
            None => Err(FractalUtilError::MultiPolyErr(
                "Nothing committed yet!".to_string(),
            )),
            Some(tree) => tree.prove_batch(&indices).map_err(|e| {
                FractalUtilError::MultiPolyErr(format!(
                    "Got an error when committing to the evals: {e}"
                ))
            }),
        }?;
        Ok((values, proof))
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn verify_values_and_proof_at(
        vals: Vec<B>,
        root: &<H>::Digest,
        proof: &[<H>::Digest],
        index: usize,
    ) -> Result<(), FractalUtilError> {
        let r = index & 1;
        if proof[r] != H::hash_elements(&vals) {
            return Err(FractalUtilError::MultiPolyErr(
                "The proof's value does not match the sent value".to_string(),
            ));
        }
        MerkleTree::<H>::verify(*root, index, proof).map_err(|e| {
            FractalUtilError::MultiPolyErr(format!(
                "Got an error when committing to the evals: {e}"
            ))
        })
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn batch_verify_values_and_proofs_at(
        vals: &Vec<Vec<E>>,
        root: &<H>::Digest,
        proof: &BatchMerkleProof<H>,
        indices: &Vec<usize>,
    ) -> Result<(), FractalUtilError> {
        //for (index, i) in indices.iter().enumerate() {
        for i in (0..indices.len()).into_iter() {
            if H::hash_elements(&vals[i]) != proof.leaves[i] {
                println!(
                    "Hash_elements applied to input array elts {:?}",
                    vals.iter()
                        .map(|x| H::hash_elements(x).as_bytes())
                        .collect::<Vec<[u8; 32]>>()
                );
                println!("Leaves {:?}", proof.leaves);
                return Err(FractalUtilError::MultiPolyErr(
                    "The proof's value does not match the sent value".to_string(),
                ));
            }
        } // TODO: still need to check this but currently leaves is private
        MerkleTree::verify_batch(root, &indices, proof).map_err(|e| {
            FractalUtilError::MultiPolyErr(format!(
                "Got an error when committing to the evals: {e}"
            ))
        })
    }

    #[cfg_attr(feature = "flame_it", flame("MultiEval"))]
    fn batch_verify_transposed_values_and_proofs_at(
        vals: Vec<[E; 1]>,
        root: &<H>::Digest,
        proof: &BatchMerkleProof<H>,
        indices: &Vec<usize>,
    ) -> Result<(), FractalUtilError> {
        //for (index, i) in indices.iter().enumerate() {
        for i in (0..indices.len()).into_iter() {
            if H::hash_elements(&vals[i]) != proof.leaves[i] {
                println!(
                    "Hash_elements applied to input array elts {:?}",
                    vals.iter()
                        .map(|x| H::hash_elements(x).as_bytes())
                        .collect::<Vec<[u8; 32]>>()
                );
                println!("Leaves {:?}", proof.leaves);
                return Err(FractalUtilError::MultiPolyErr(
                    "The proof's value does not match the sent value".to_string(),
                ));
            }
        } // TODO: still need to check this but currently leaves is private
        MerkleTree::verify_batch(root, &indices, proof).map_err(|e| {
            FractalUtilError::MultiPolyErr(format!(
                "Got an error when committing to the evals: {e}"
            ))
        })
    }
}

#[cfg_attr(feature = "flame_it", flame("polynomial_utils"))]
pub fn eval_on_domain<B, E>(
    coefficients: &[E],
    evaluation_domain_len: usize,
    eval_twiddles: &[B],
) -> Vec<E>
where
    B: StarkField,
    E: FieldElement<BaseField = B>,
{
    let mut eval = Vec::from(coefficients);
    pad_with_zeroes(&mut eval, evaluation_domain_len);
    fft::evaluate_poly(&mut eval, eval_twiddles);

    eval
}
