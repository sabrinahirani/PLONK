use ark_bn254::Fr;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain, Polynomial};
use ark_poly::univariate::DensePolynomial;

use crate::transcript::PlonkTranscript;
use crate::kzg::KZG;
use crate::constraints::{permutation_check_polynomial, zero_test_check};

pub struct ProofData {
    pub commitment_a: Vec<u8>,
    pub commitment_b: Vec<u8>,
    pub commitment_c: Vec<u8>,
    pub commitment_z: Vec<u8>,
    pub commitment_t: Vec<u8>,
    pub proof_openings: Vec<u8>,
}

pub fn verify(
    proof: &ProofData,
    circuit_size: usize,
    kzg: &KZG,
    max_degree: usize,
) -> bool {
    let domain = Radix2EvaluationDomain::<Fr>::new(circuit_size)
        .expect("domain size must be power of two");

    // Rebuild transcript and append commitments (same order as prover)
    let mut transcript = PlonkTranscript::new(b"plonk");
    transcript.append_message(b"commitment_a", &proof.commitment_a);
    transcript.append_message(b"commitment_b", &proof.commitment_b);
    transcript.append_message(b"commitment_c", &proof.commitment_c);

    // Derive challenges identically to prover
    let beta = transcript.get_challenge(b"beta");
    let gamma = transcript.get_challenge(b"gamma");
    let alpha = transcript.get_challenge(b"alpha");
    let zeta = transcript.get_challenge(b"zeta");

    // Append remaining commitments
    transcript.append_message(b"commitment_z", &proof.commitment_z);
    transcript.append_message(b"commitment_t", &proof.commitment_t);

    // Deserialize commitments from bytes (you must implement this)
    let commitment_a = kzg.deserialize_commitment(&proof.commitment_a);
    let commitment_b = kzg.deserialize_commitment(&proof.commitment_b);
    let commitment_c = kzg.deserialize_commitment(&proof.commitment_c);
    let commitment_z = kzg.deserialize_commitment(&proof.commitment_z);
    let commitment_t = kzg.deserialize_commitment(&proof.commitment_t);

    // Deserialize opening proof (could be multiple proofs combined)
    let proof_opening = kzg.deserialize_proof(&proof.proof_openings);

    // Verify openings at zeta
    // For simplicity, this example assumes a combined proof or multiple verifications

    let powers = kzg.trim(max_degree);

    // Verify openings of all committed polynomials at zeta
    let valid_a = kzg.verify(&commitment_a, zeta, /*expected eval*/Fr::zero(), &proof_opening);
    let valid_b = kzg.verify(&commitment_b, zeta, Fr::zero(), &proof_opening);
    let valid_c = kzg.verify(&commitment_c, zeta, Fr::zero(), &proof_opening);
    let valid_z = kzg.verify(&commitment_z, zeta, Fr::zero(), &proof_opening);
    let valid_t = kzg.verify(&commitment_t, zeta, Fr::zero(), &proof_opening);

    // For real use, expected evaluations must be provided (from polynomial evaluations or challenge)
    // Here we return true only if all proofs verify
    valid_a && valid_b && valid_c && valid_z && valid_t
}
