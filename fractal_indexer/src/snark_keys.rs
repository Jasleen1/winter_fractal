use std::marker::PhantomData;

use crate::{
    errors::*,
    index::{create_index_from_r1cs, Index, IndexParams},
    indexed_matrix::IndexedMatrix,
};
//use fri::utils::hash_values;
use models::r1cs::{Matrix, R1CS};
use winter_crypto::{BatchMerkleProof, ElementHasher, Hasher, MerkleTree, MerkleTreeError};
use winter_fri::utils::hash_values;
use winter_math::{polynom, FieldElement, StarkField};
use winter_utils::transpose_slice;

#[derive(Debug, Clone)] // Clone
pub struct ProverIndexPolynomial<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub polynomial: Vec<B>,
    pub evaluations: Vec<E>,
    pub tree: MerkleTree<H>,
}

// impl<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> Clone for ProverIndexPolynomial<H, B> {
//     fn clone(&self) -> Self {
//         ProverIndexPolynomial {
//             polynomial: self.polynomial.clone(),
//             evaluations: self.evaluations.clone(),
//             tree: self.tree.clone()
//         }
//     }
// }

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > ProverIndexPolynomial<B, E, H>
{
    // TODO Add error checking, currently assumes index is
    // within range.
    pub fn get_eval_at_index(&self, index: usize) -> E {
        self.evaluations[index]
    }

    pub fn get_eval_at_point(&self, point: E) -> E {
        polynom::eval(&self.polynomial, point)
        //unimplemented!()
    }
}

#[derive(Debug, Clone)] // Clone
pub struct ProverMatrixIndex<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub matrix: Matrix<B>,
    pub row_poly: ProverIndexPolynomial<B, E, H>,
    pub col_poly: ProverIndexPolynomial<B, E, H>,
    pub val_poly: ProverIndexPolynomial<B, E, H>,
}

// impl<H: ElementHasher + ElementHasher<BaseField = B>, B: StarkField> Clone for ProverMatrixIndex<H, B> {
//     fn clone(&self) -> Self {
//         ProverMatrixIndex {
//             matrix: self.matrix.clone(),
//             row_poly: self.row_poly.clone(),
//             col_poly: self.col_poly.clone(),
//             val_poly: self.val_poly.clone()
//         }
//     }
// }

impl<
        B: StarkField,
        E: FieldElement<BaseField = B>,
        H: ElementHasher + ElementHasher<BaseField = B>,
    > ProverMatrixIndex<B, E, H>
{
    pub fn get_val_eval(&self, point: E) -> E {
        self.val_poly.get_eval_at_point(point)
    }
    pub fn get_val_eval_at_index(&self, index: usize) -> E {
        self.val_poly.get_eval_at_index(index)
    }

    pub fn get_col_eval(&self, point: E) -> E {
        self.col_poly.get_eval_at_point(point)
    }
    pub fn get_col_eval_at_index(&self, index: usize) -> E {
        self.col_poly.get_eval_at_index(index)
    }

    pub fn get_row_eval(&self, point: E) -> E {
        self.row_poly.get_eval_at_point(point)
    }
    pub fn get_row_eval_at_index(&self, index: usize) -> E {
        self.row_poly.get_eval_at_index(index)
    }

    pub fn decommit_row_eval(
        &self,
        queries: &Vec<usize>,
    ) -> Result<(Vec<E>, BatchMerkleProof<H>), MerkleTreeError> {
        let values = queries
            .iter()
            .map(|q| self.row_poly.evaluations[*q])
            .collect::<Vec<E>>()
            .to_vec();
        let proof = self.row_poly.tree.prove_batch(queries)?;
        Ok((values, proof))
    }

    pub fn decommit_transposed_row_eval(
        &self,
        queries: &Vec<usize>,
    ) -> Result<(Vec<[E; 1]>, BatchMerkleProof<H>), MerkleTreeError> {
        let evals: Vec<[E; 1]> = transpose_slice(&self.row_poly.evaluations.clone());
        let values = queries
            .iter()
            .map(|q| evals[*q])
            .collect::<Vec<[E; 1]>>()
            .to_vec();
        let proof = self.row_poly.tree.prove_batch(queries)?;
        Ok((values, proof))
    }

    pub fn decommit_col_eval(
        &self,
        queries: &Vec<usize>,
    ) -> Result<(Vec<E>, BatchMerkleProof<H>), MerkleTreeError> {
        let values = queries
            .iter()
            .map(|q| self.col_poly.evaluations[*q])
            .collect::<Vec<E>>()
            .to_vec();
        let proof = self.col_poly.tree.prove_batch(queries)?;
        Ok((values, proof))
    }

    pub fn decommit_val_eval(
        &self,
        queries: &Vec<usize>,
    ) -> Result<(Vec<E>, BatchMerkleProof<H>), MerkleTreeError> {
        let values = queries
            .iter()
            .map(|q| self.val_poly.evaluations[*q])
            .collect::<Vec<E>>()
            .to_vec();
        let proof = self.val_poly.tree.prove_batch(queries)?;
        Ok((values, proof))
    }

    pub fn decommit_evals(
        &self,
        queries: Vec<usize>,
    ) -> Result<Vec<(Vec<E>, BatchMerkleProof<H>)>, MerkleTreeError> {
        let row_evals = self.decommit_row_eval(&queries)?;
        let col_evals = self.decommit_col_eval(&queries)?;
        let val_evals = self.decommit_val_eval(&queries)?;
        Ok(vec![row_evals, col_evals, val_evals])
    }
    pub fn decommit_proofs(
        &self,
        queries: Vec<usize>,
    ) -> Result<Vec<BatchMerkleProof<H>>, MerkleTreeError> {
        let row_evals = self.row_poly.tree.prove_batch(&queries)?;
        let col_evals = self.col_poly.tree.prove_batch(&queries)?;
        let val_evals = self.val_poly.tree.prove_batch(&queries)?;
        Ok(vec![row_evals, col_evals, val_evals])
    }

    pub fn decommit_proof(
        &self,
        queries: Vec<usize>,
        matrix_idx: usize,
    ) -> Result<BatchMerkleProof<H>, MerkleTreeError> {
        match matrix_idx {
            0 => Ok(self.row_poly.tree.prove_batch(&queries)?),
            1 => Ok(self.col_poly.tree.prove_batch(&queries)?),
            2 => Ok(self.val_poly.tree.prove_batch(&queries)?),
            _ => Err(MerkleTreeError::InvalidProof),
        }
    }
}

#[derive(Debug, Clone)] // Clone
pub struct ProverKey<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub params: IndexParams<B>,
    pub matrix_a_index: ProverMatrixIndex<B, E, H>,
    pub matrix_b_index: ProverMatrixIndex<B, E, H>,
    pub matrix_c_index: ProverMatrixIndex<B, E, H>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifierMatrixIndex<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub row_poly_commitment: H::Digest,
    pub col_poly_commitment: H::Digest,
    pub val_poly_commitment: H::Digest,
    pub _phantom_e: PhantomData<E>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifierKey<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
> {
    pub params: IndexParams<B>,
    pub matrix_a_commitments: VerifierMatrixIndex<B, E, H>,
    pub matrix_b_commitments: VerifierMatrixIndex<B, E, H>,
    pub matrix_c_commitments: VerifierMatrixIndex<B, E, H>,
}

// QUESTION: Currently using the utils hash_values function which uses quartic folding.
// Is there any drawback to doing this here, where there's no layering?
pub fn commit_polynomial_evaluations<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
>(
    evaluations: &Vec<E>,
) -> Result<MerkleTree<H>, IndexerError> {
    let transposed_evaluations = transpose_slice(evaluations);
    let hashed_evaluations = hash_values::<H, E, N>(&transposed_evaluations);
    Ok(MerkleTree::<H>::new(hashed_evaluations)?)
}

pub fn generate_prover_and_verifier_matrix_index<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
>(
    indexed: IndexedMatrix<B, E>,
) -> Result<(ProverMatrixIndex<B, E, H>, VerifierMatrixIndex<B, E, H>), IndexerError> {
    let matrix = indexed.matrix;
    let row_polynomial = indexed.row_poly;
    let col_polynomial = indexed.col_poly;
    let val_polynomial = indexed.val_poly;
    let row_evals = indexed.row_evals_on_l;
    let col_evals = indexed.col_evals_on_l;
    let val_evals = indexed.val_evals_on_l;
    let row_tree = commit_polynomial_evaluations::<B, E, H, N>(&row_evals)?;
    let col_tree = commit_polynomial_evaluations::<B, E, H, N>(&col_evals)?;
    let val_tree = commit_polynomial_evaluations::<B, E, H, N>(&val_evals)?;
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
        _phantom_e: PhantomData::<E>,
    };
    Ok((prover_matrix_index, verifier_matrix_index))
}

pub fn generate_prover_and_verifier_keys<
    B: StarkField,
    E: FieldElement<BaseField = B>,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
>(
    Index {
        params,
        indexed_a,
        indexed_b,
        indexed_c,
    }: Index<B, E>,
) -> Result<(ProverKey<B, E, H>, VerifierKey<B, E, H>), IndexerError> {
    let (matrix_a_index, matrix_a_commitments) =
        generate_prover_and_verifier_matrix_index::<B, E, H, N>(indexed_a)?;
    let (matrix_b_index, matrix_b_commitments) =
        generate_prover_and_verifier_matrix_index::<B, E, H, N>(indexed_b)?;
    let (matrix_c_index, matrix_c_commitments) =
        generate_prover_and_verifier_matrix_index::<B, E, H, N>(indexed_c)?;
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
    B: StarkField,
    H: ElementHasher + ElementHasher<BaseField = B>,
    const N: usize,
>(
    params: IndexParams<B>,
    r1cs_instance: R1CS<B>,
) -> Result<(ProverKey<B, B, H>, VerifierKey<B, B, H>), IndexerError> {
    let index = create_index_from_r1cs(params, r1cs_instance);
    generate_prover_and_verifier_keys::<B, B, H, N>(index)
}
