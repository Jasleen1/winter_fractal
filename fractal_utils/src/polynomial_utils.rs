use crate::{errors::FractalUtilError, matrix_utils::*};
use fractal_math::{fft, FieldElement, StarkField};
use winter_fri::{DefaultProverChannel, FriOptions};
use std::{convert::TryInto, marker::PhantomData};
use winter_crypto::{BatchMerkleProof, ElementHasher, MerkleTree};
use winter_utils::batch_iter_mut;
// TODO: Add error checking and throwing
/**
 * This is equivalent to computing v_H(X) for a multiplicative coset
 * H = eta * H_0 for a multiplicative subgroup H_0 of order dom_size
 * Note that v_H(X) = X^dom_size - eta^dom_size. If eta = 1 we're
 * in a multiplicative subgroup itself.
 **/
pub fn compute_vanishing_poly<E: FieldElement>(x: E, eta: E, dom_size: usize) -> E {
    let power_u64: u64 = dom_size.try_into().unwrap();
    let power = E::PositiveInteger::from(power_u64);
    x.exp(power) - eta.exp(power)
}

/// This function generates the vanshing polynomial coefficients for a multiplicative
/// subgroup of size dom_size and with multiplicative factor eta.
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

pub fn compute_binomial_on_x<E: FieldElement>(bivariate: BivariatePoly<E>, x_val: E) -> Vec<E> {
    // Given a BivariatePoly, computes a monomial in Y, obtained by evaluating bivariate(x_val, Y)
    // represented by the output vector
    let bivariate_as_matrix = Matrix::new("binomial", bivariate).unwrap();
    let transposed_bivariate = bivariate_as_matrix.get_transpose("transposed_binomial");
    compute_binomial_on_y(transposed_bivariate.mat, x_val)
}

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
    for _ in 0..diff {
        poly.push(E::ZERO);
    }
}

pub fn get_to_degree_size<E: FieldElement>(poly: &mut Vec<E>) {
    let len_poly = poly.len();
    let mut count = len_poly - 1;
    while count > 0 && poly[count] == E::ZERO {
        poly.pop();
        count = count - 1;
    }
}

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
    fn get_values_at(&self, index: usize) -> Result<Vec<B>, FractalUtilError>;
    /// This function retrieves the evals of the polynomials at a set of evaluation points.
    fn batch_get_values_at(&self, indices: Vec<usize>) -> Result<Vec<Vec<B>>, FractalUtilError>;
    /// This function takes as input an index of a point in the evaluation domain and
    /// outputs the evals committed at that point and a proof.
    fn get_values_and_proof_at(
        &self,
        index: usize,
    ) -> Result<(Vec<B>, Vec<H::Digest>), FractalUtilError>;
    /// This function takes as input the indices of multiple points in the evaluation domain and
    /// returns the evaluations of all the polynomials at these points, together with a batch merkle
    /// proof showing that this eval was done correctly.
    fn batch_get_values_and_proofs_at(
        &self,
        indices: Vec<usize>,
    ) -> Result<(Vec<Vec<B>>, BatchMerkleProof<H>), FractalUtilError>;
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
        vals: Vec<Vec<B>>,
        root: &H::Digest,
        proof: &BatchMerkleProof<H>,
        indices: Vec<usize>,
    ) -> Result<(), FractalUtilError>;
}

pub struct MultiEval<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub evaluations: Vec<Vec<B>>,
    pub coefficients: Vec<Vec<B>>,
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
    pub fn new(coefficients: Vec<Vec<B>>, evaluation_domain_len: usize, offset: B) -> Self {
        let eval_twiddles = fft::get_twiddles(evaluation_domain_len);
        let mut accumulated_evals = Vec::<Vec<B>>::new();
        for (_, poly) in coefficients.iter().enumerate() {
            accumulated_evals.push(Self::eval_on_domain(
                poly.to_vec(),
                evaluation_domain_len,
                eval_twiddles.clone(),
            ));
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

    pub fn add_polynomial(&mut self, coefficients: Vec<B>, evaluation_domain_len: usize) -> (){
        let eval_twiddles = fft::get_twiddles(evaluation_domain_len);
        let evaluations = Self::eval_on_domain(
            coefficients.clone(),
            evaluation_domain_len,
            eval_twiddles.clone(),
        );
        self.coefficients.push(coefficients);
        self.evaluations.push(evaluations);
        self.committed_tree =  Option::None;
    }

    /// This is mostly a helper function to evaluate the polynomials on a domain of given length
    /// for which twiddles are already computed.
    pub fn eval_on_domain(
        coefficients: Vec<B>,
        evaluation_domain_len: usize,
        eval_twiddles: Vec<B>,
    ) -> Vec<B> {
        let mut eval = coefficients.clone();
        pad_with_zeroes(&mut eval, evaluation_domain_len);

        fft::evaluate_poly(&mut eval, &mut eval_twiddles.clone());

        eval
    }

    /// Helper function to zip the evaluations so that each element of the output is of the
    /// form [poly_1(e), ..., poly_n(e)] i.e. evaluations of all the polynomials are included
    /// in the same array.
    fn zip_evals(separate_evals: Vec<Vec<B>>, evaluation_domain_len: usize) -> Vec<Vec<B>> {
        let mut zipped_evals = vec![Vec::<B>::new(); evaluation_domain_len];
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

    fn get_commitment(&self) -> Result<&<H as winter_crypto::Hasher>::Digest, FractalUtilError> {
        match &self.committed_tree {
            Some(merkle_tree) => Ok(merkle_tree.root()),
            None => Err(FractalUtilError::MultiPolyErr(
                "The Merkle Tree in the poly commit is None.".to_string(),
            )),
        }
    }

    fn get_values_at(&self, index: usize) -> Result<Vec<B>, FractalUtilError> {
        Ok(self.evaluations[index].clone())
    }

    fn batch_get_values_at(&self, indices: Vec<usize>) -> Result<Vec<Vec<B>>, FractalUtilError> {
        let mut output_vals = Vec::<Vec<B>>::new();
        for (_, &index) in indices.iter().enumerate() {
            output_vals.push(self.evaluations[index].clone());
        }
        Ok(output_vals)
    }

    fn get_values_and_proof_at(
        &self,
        index: usize,
    ) -> Result<(Vec<B>, Vec<<H>::Digest>), FractalUtilError> {
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

    fn batch_get_values_and_proofs_at(
        &self,
        indices: Vec<usize>,
    ) -> Result<(Vec<Vec<B>>, BatchMerkleProof<H>), FractalUtilError> {
        let values = self.batch_get_values_at(indices.clone())?;
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

    fn batch_verify_values_and_proofs_at(
        _vals: Vec<Vec<B>>,
        root: &<H>::Digest,
        proof: &BatchMerkleProof<H>,
        indices: Vec<usize>,
    ) -> Result<(), FractalUtilError> {
        // for index in indices {
        //     if H::hash_elements(&vals[index]) != proof.leaves[index] {
        //         return Err(FractalUtilError::MultiPolyErr("The proof's value does not match the sent value".to_string()));
        //     }
        // } // TODO: still need to check this but currently leaves is private
        MerkleTree::verify_batch(root, &indices, proof).map_err(|e| {
            FractalUtilError::MultiPolyErr(format!(
                "Got an error when committing to the evals: {e}"
            ))
        })
    }
}