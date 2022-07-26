pub mod errors;

#[cfg(test)]
mod tests;

pub fn create_index(r1cs: R1CS, params: IndexParams) -> Index {
    let domains = build_index_domains(params);
    Index(
        params,
        IndexedMatrix::new(r1cs.A, domains),
        IndexedMatrix::new(r1cs.B, domains),
        IndexedMatrix::new(r1cs.C, domains)
    )
}
