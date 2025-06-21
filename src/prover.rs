use ark_bn254::Fr;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain, Polynomial};
use ark_poly::univariate::DensePolynomial;

use crate::circuit::structure::Circuit;
use crate::constraints::{permutation_check_polynomial, zero_test_check};
use crate::transcript::PlonkTranscript;
use crate::kzg::KZG;

pub struct ProofData {
    pub commitment_a: Vec<u8>,
    pub commitment_b: Vec<u8>,
    pub commitment_c: Vec<u8>,
    pub commitment_z: Vec<u8>,
    pub commitment_t: Vec<u8>, // quotient polynomial commitment
    pub proof_openings: Vec<u8>, // combined opening proof (could be multiple in practice)
}

pub fn prove(
    circuit: &mut Circuit,
    kzg: &KZG,
    max_degree: usize,
) -> ProofData {
    // 1. Setup domain and synthesize
    let n = circuit.gates.len();
    let domain = Radix2EvaluationDomain::<Fr>::new(n)
        .expect("domain size must be power of two and >= circuit size");

    circuit.synthesize_r1cs();

    // 2. Convert R1CS vectors to polynomials
    let (a_poly, b_poly, c_poly) = circuit.r1cs_to_polynomials(domain);

    // 3. Commit to a, b, c
    let powers = kzg.trim(max_degree);

    let commitment_a = kzg.commit(&a_poly, &powers);
    let commitment_b = kzg.commit(&b_poly, &powers);
    let commitment_c = kzg.commit(&c_poly, &powers);

    // 4. Fiat-Shamir transcript
    let mut transcript = PlonkTranscript::new(b"plonk");

    transcript.append_message(b"commitment_a", &commitment_a.to_bytes());
    transcript.append_message(b"commitment_b", &commitment_b.to_bytes());
    transcript.append_message(b"commitment_c", &commitment_c.to_bytes());

    let beta = transcript.get_challenge(b"beta");
    let gamma = transcript.get_challenge(b"gamma");
    let alpha = transcript.get_challenge(b"alpha");
    let zeta = transcript.get_challenge(b"zeta");

    // 5. Prepare permutation polynomials (for demo, use a, b, c clones)
    let s1_poly = a_poly.clone();
    let s2_poly = b_poly.clone();
    let s3_poly = c_poly.clone();

    // 6. Compute grand product polynomial Z (replace with your actual method)
    let z_poly = DensePolynomial::from_coefficients_vec(vec![Fr::one(); n]);

    // 7. Commit to Z
    let commitment_z = kzg.commit(&z_poly, &powers);
    transcript.append_message(b"commitment_z", &commitment_z.to_bytes());

    // 8. Compute permutation check polynomial
    let perm_check_poly = permutation_check_polynomial(
        &a_poly,
        &b_poly,
        &c_poly,
        &s1_poly,
        &s2_poly,
        &s3_poly,
        &z_poly,
        domain,
        beta,
        gamma,
        alpha,
    );

    // 9. Check zero test (optional in prover, but good sanity check)
    assert!(zero_test_check(perm_check_poly.clone(), domain), "Permutation check failed zero test");

    // 10. Compute quotient polynomial t(x) = (arithmetic + permutation + etc) / Z_H(x)
    // Here you would sum all your constraint polynomials and divide by vanishing polynomial
    // For demo, let's pretend perm_check_poly is your quotient
    let quotient_poly = perm_check_poly; // TODO: add arithmetic gate checks too!

    // 11. Commit to quotient polynomial t(x)
    let commitment_t = kzg.commit(&quotient_poly, &powers);
    transcript.append_message(b"commitment_t", &commitment_t.to_bytes());

    // 12. Open polynomials at challenge point zeta
    // These calls generate opening proofs of committed polynomials at point zeta
    let proof_a = kzg.open(&a_poly, zeta, &powers);
    let proof_b = kzg.open(&b_poly, zeta, &powers);
    let proof_c = kzg.open(&c_poly, zeta, &powers);
    let proof_z = kzg.open(&z_poly, zeta, &powers);
    let proof_t = kzg.open(&quotient_poly, zeta, &powers);

    // 13. Combine or return proofs (combine in practice, return vector here)
    // You could aggregate proofs for optimization (not shown here)
    let mut proof_openings = Vec::new();
    proof_openings.extend(proof_a.to_bytes());
    proof_openings.extend(proof_b.to_bytes());
    proof_openings.extend(proof_c.to_bytes());
    proof_openings.extend(proof_z.to_bytes());
    proof_openings.extend(proof_t.to_bytes());

    ProofData {
        commitment_a: commitment_a.to_bytes(),
        commitment_b: commitment_b.to_bytes(),
        commitment_c: commitment_c.to_bytes(),
        commitment_z: commitment_z.to_bytes(),
        commitment_t: commitment_t.to_bytes(),
        proof_openings,
    }
}
