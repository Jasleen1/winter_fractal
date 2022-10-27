use crate::{
    errors::*,
    index::{create_index_from_r1cs, Index, IndexParams},
    indexed_matrix::IndexedMatrix,
};
//use fri::utils::hash_values;
use models::r1cs::{Matrix, R1CS};
use winter_crypto::{ElementHasher, MerkleTree};
use winter_fri::utils::hash_values;
use winter_math::{polynom, FieldElement, StarkField};
use winter_utils::transpose_slice;

#[derive(Debug)] // Clone
pub struct ProverIndexPolynomial<H: ElementHasher + ElementHasher<BaseField = E>, E: FieldElement> {
    pub polynomial: Vec<E>,
    pub evaluations: Vec<E>,
    pub tree: MerkleTree<H>,
}

impl<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> ProverIndexPolynomial<H, B> {
    // TODO Add error checking, currently assumes index is
    // within range.
    pub fn get_eval_at_index(&self, index: usize) -> B {
        self.evaluations[index]
    }

    pub fn get_eval_at_point(&self, point: B) -> B {
        polynom::eval(&self.polynomial, point)
        //unimplemented!()
    }
}

#[derive(Debug)] // Clone
pub struct ProverMatrixIndex<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> {
    pub matrix: Matrix<B>,
    pub row_poly: ProverIndexPolynomial<H, B>,
    pub col_poly: ProverIndexPolynomial<H, B>,
    pub val_poly: ProverIndexPolynomial<H, B>,
}

impl<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> ProverMatrixIndex<H, B> {
    pub fn get_val_eval(&self, point: B) -> B {
        self.val_poly.get_eval_at_point(point)
    }
    pub fn get_val_eval_at_index(&self, index: usize) -> B {
        self.val_poly.get_eval_at_index(index)
    }

    pub fn get_col_eval(&self, point: B) -> B {
        self.col_poly.get_eval_at_point(point)
    }
    pub fn get_col_eval_at_index(&self, index: usize) -> B {
        self.col_poly.get_eval_at_index(index)
    }

    pub fn get_row_eval(&self, point: B) -> B {
        self.row_poly.get_eval_at_point(point)
    }
    pub fn get_row_eval_at_index(&self, index: usize) -> B {
        self.row_poly.get_eval_at_index(index)
    }
}

#[derive(Debug)] // Clone
pub struct ProverKey<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> {
    pub params: IndexParams<B>,
    pub matrix_a_index: ProverMatrixIndex<H, B>,
    pub matrix_b_index: ProverMatrixIndex<H, B>,
    pub matrix_c_index: ProverMatrixIndex<H, B>,
}

#[derive(Debug, Clone)]
pub struct VerifierMatrixIndex<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> {
    pub row_poly_commitment: H::Digest,
    pub col_poly_commitment: H::Digest,
    pub val_poly_commitment: H::Digest,
}

#[derive(Debug, Clone)]
pub struct VerifierKey<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> {
    pub params: IndexParams<B>,
    pub matrix_a_commitments: VerifierMatrixIndex<H, B>,
    pub matrix_b_commitments: VerifierMatrixIndex<H, B>,
    pub matrix_c_commitments: VerifierMatrixIndex<H, B>,
}

// QUESTION: Currently using the utils hash_values function which uses quartic folding.
// Is there any drawback to doing this here, where there's no layering?
pub fn commit_polynomial_evaluations<
    H: ElementHasher + ElementHasher<BaseField = B>,
    B: StarkField,
    const N: usize,
>(
    evaluations: &Vec<B>,
) -> Result<MerkleTree<H>, IndexerError> {
    let transposed_evaluations = transpose_slice(evaluations);
    let hashed_evaluations = hash_values::<H, B, N>(&transposed_evaluations);
    Ok(MerkleTree::<H>::new(hashed_evaluations)?)
}

pub fn generate_prover_and_verifier_matrix_index<
    H: ElementHasher + ElementHasher<BaseField = B>,
    B: StarkField,
    const N: usize,
>(
    indexed: IndexedMatrix<B>,
) -> Result<(ProverMatrixIndex<H, B>, VerifierMatrixIndex<H, B>), IndexerError> {
    let matrix = indexed.matrix;
    let row_polynomial = indexed.row_poly;
    let col_polynomial = indexed.col_poly;
    let val_polynomial = indexed.val_poly;
    let row_evals = indexed.row_evals_on_l;
    let col_evals = indexed.col_evals_on_l;
    let val_evals = indexed.val_evals_on_l;
    let row_tree = commit_polynomial_evaluations::<H, B, N>(&row_evals)?;
    let col_tree = commit_polynomial_evaluations::<H, B, N>(&col_evals)?;
    let val_tree = commit_polynomial_evaluations::<H, B, N>(&val_evals)?;
    let row_poly_commitment = *row_tree.root();
    let col_poly_commitment = *col_tree.root();
    let val_poly_commitment = *val_tree.root();

    let row_poly = ProverIndexPolynomial {
        polynomial: row_polynomial,
        evaluations: row_evals,
        tree: row_tree,
    };
    let col_poly = ProverIndexPolynomial {
        polynomial: col_polynomial,
        evaluations: col_evals,
        tree: col_tree,
    };
    let val_poly = ProverIndexPolynomial {
        polynomial: val_polynomial,
        evaluations: val_evals,
        tree: val_tree,
    };
    let prover_matrix_index = ProverMatrixIndex {
        matrix,
        row_poly,
        col_poly,
        val_poly,
    };
    let verifier_matrix_index = VerifierMatrixIndex {
        row_poly_commitment,
        col_poly_commitment,
        val_poly_commitment,
    };
    Ok((prover_matrix_index, verifier_matrix_index))
}

pub fn generate_prover_and_verifier_keys<
    H: ElementHasher + ElementHasher<BaseField = B>,
    B: StarkField,
    const N: usize,
>(
    Index {
        params,
        indexed_a,
        indexed_b,
        indexed_c,
    }: Index<B>,
) -> Result<(ProverKey<H, B>, VerifierKey<H, B>), IndexerError> {
    let (matrix_a_index, matrix_a_commitments) =
        generate_prover_and_verifier_matrix_index::<H, B, N>(indexed_a)?;
    let (matrix_b_index, matrix_b_commitments) =
        generate_prover_and_verifier_matrix_index::<H, B, N>(indexed_b)?;
    let (matrix_c_index, matrix_c_commitments) =
        generate_prover_and_verifier_matrix_index::<H, B, N>(indexed_c)?;
    Ok((
        ProverKey {
            params: params.clone(),
            matrix_a_index,
            matrix_b_index,
            matrix_c_index,
        },
        VerifierKey {
            params,
            matrix_a_commitments,
            matrix_b_commitments,
            matrix_c_commitments,
        },
    ))
}

pub fn generate_basefield_keys<
    H: ElementHasher + ElementHasher<BaseField = B>,
    B: StarkField,
    const N: usize,
>(
    params: IndexParams<B>,
    r1cs_instance: R1CS<B>,
) -> Result<(ProverKey<H, B>, VerifierKey<H, B>), IndexerError> {
    let index = create_index_from_r1cs(params, r1cs_instance);
    generate_prover_and_verifier_keys::<H, B, N>(index)
}
