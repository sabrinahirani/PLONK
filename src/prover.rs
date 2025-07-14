use std::collections::HashMap;
use ark_ff::{Field, FftField};
use ark_poly::{
    EvaluationDomain, GeneralEvaluationDomain,
    univariate::DensePolynomial, Polynomial,
};
use ark_poly_commit::{
    marlin_pc::MarlinKZG10,
    LabeledPolynomial,
    PolynomialCommitment,
};
use ark_crypto_primitives::sponge::poseidon::{PoseidonSponge, PoseidonConfig};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use ark_std::rand::Rng;
use ark_std::test_rng;
use ark_serialize::CanonicalSerialize;

use crate::circuit::{WireColumn, WirePosition, PermutationLayout};
use crate::transcript::{PlonkTranscript, absorb_public_inputs};

use ark_bn254::{Bn254, Fr};

use ark_poly::DenseUVPolynomial;

// type aliases for convenience
type UniPoly = DensePolynomial<Fr>;
type PCS = MarlinKZG10<Bn254, UniPoly>;
pub type KZGCommitment = <PCS as PolynomialCommitment<Fr, UniPoly>>::Commitment;
pub type KZGOpeningProof = <PCS as PolynomialCommitment<Fr, UniPoly>>::Proof;
pub type KZGProverKey = <PCS as PolynomialCommitment<Fr, UniPoly>>::CommitterKey;
pub type KZGVerifierKey = <PCS as PolynomialCommitment<Fr, UniPoly>>::VerifierKey;

/// Create a test sponge for the MarlinKZG10 scheme
fn test_sponge<F: ark_ff::PrimeField>() -> PoseidonSponge<F> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut v = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();
        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}

pub struct PlonkProof<F: Field> {

    // Evaluations at zeta:
    pub eval_a: F,
    pub eval_b: F,
    pub eval_c: F,
    pub eval_q_add: F,
    pub eval_q_mul: F,
    pub eval_s_id: F,
    pub eval_s_sigma: F,
    pub eval_z: F,
    pub eval_t: F,

    // Commitments:
    pub comm_a: KZGCommitment,
    pub comm_b: KZGCommitment,
    pub comm_c: KZGCommitment,
    pub comm_q_add: KZGCommitment,
    pub comm_q_mul: KZGCommitment,
    pub comm_s_id: KZGCommitment,
    pub comm_s_sigma: KZGCommitment,
    pub comm_z: KZGCommitment,
    pub comm_t: KZGCommitment,

    // Opening proofs:
    pub proof_a: KZGOpeningProof,
    pub proof_b: KZGOpeningProof,
    pub proof_c: KZGOpeningProof,
    pub proof_q_add: KZGOpeningProof,
    pub proof_q_mul: KZGOpeningProof,
    pub proof_s_id: KZGOpeningProof,
    pub proof_s_sigma: KZGOpeningProof,
    pub proof_z: KZGOpeningProof,
    pub proof_t: KZGOpeningProof,
}

pub fn evaluate_polynomials_at_zeta<F: Field>(
    zeta: F,
    a: &DensePolynomial<F>,
    b: &DensePolynomial<F>,
    c: &DensePolynomial<F>,
    s_id: &DensePolynomial<F>,
    s_sigma: &DensePolynomial<F>,
    z: &DensePolynomial<F>,
    t: &DensePolynomial<F>,
) -> (F, F, F, F, F, F, F) {
    (
        a.evaluate(&zeta),
        b.evaluate(&zeta),
        c.evaluate(&zeta),
        s_id.evaluate(&zeta),
        s_sigma.evaluate(&zeta),
        z.evaluate(&zeta),
        t.evaluate(&zeta),
    )
}

/// Create a PLONK proof using non-interactive challenge generation via Fiat-Shamir
pub fn create_plonk_proof_with_transcript<R: Rng>(
    pk: &KZGProverKey,
    a: &DensePolynomial<Fr>,
    b: &DensePolynomial<Fr>,
    c: &DensePolynomial<Fr>,
    q_add: &DensePolynomial<Fr>,
    q_mul: &DensePolynomial<Fr>,
    s_id: &DensePolynomial<Fr>,
    s_sigma: &DensePolynomial<Fr>,
    public_inputs: &[Fr],
    witness_flat: &[Fr],
    sigma: &[usize],
    domain: &ark_poly::GeneralEvaluationDomain<Fr>,
    rng: &mut R,
) -> PlonkProof<Fr> {
    // Initialize transcript
    let mut transcript = PlonkTranscript::new(b"plonk_proof");
    
    // Absorb public inputs first
    absorb_public_inputs(transcript.inner(), public_inputs);
    
    // Create labeled polynomials for selector and witness polynomials
    let labeled_a = LabeledPolynomial::new("a".to_string(), a.clone(), None, Some(1));
    let labeled_b = LabeledPolynomial::new("b".to_string(), b.clone(), None, Some(1));
    let labeled_c = LabeledPolynomial::new("c".to_string(), c.clone(), None, Some(1));
    let labeled_q_add = LabeledPolynomial::new("q_add".to_string(), q_add.clone(), None, Some(1));
    let labeled_q_mul = LabeledPolynomial::new("q_mul".to_string(), q_mul.clone(), None, Some(1));
    let labeled_s_id = LabeledPolynomial::new("s_id".to_string(), s_id.clone(), None, Some(1));
    let labeled_s_sigma = LabeledPolynomial::new("s_sigma".to_string(), s_sigma.clone(), None, Some(1));

    // Commit to selector and witness polynomials first
    let (comms_a, states_a) = PCS::commit(pk, [&labeled_a], Some(rng)).unwrap();
    let (comms_b, states_b) = PCS::commit(pk, [&labeled_b], Some(rng)).unwrap();
    let (comms_c, states_c) = PCS::commit(pk, [&labeled_c], Some(rng)).unwrap();
    let (comms_q_add, states_q_add) = PCS::commit(pk, [&labeled_q_add], Some(rng)).unwrap();
    let (comms_q_mul, states_q_mul) = PCS::commit(pk, [&labeled_q_mul], Some(rng)).unwrap();
    let (comms_s_id, states_s_id) = PCS::commit(pk, [&labeled_s_id], Some(rng)).unwrap();
    let (comms_s_sigma, states_s_sigma) = PCS::commit(pk, [&labeled_s_sigma], Some(rng)).unwrap();

    let comm_a = comms_a[0].commitment().clone();
    let comm_b = comms_b[0].commitment().clone();
    let comm_c = comms_c[0].commitment().clone();
    let comm_q_add = comms_q_add[0].commitment().clone();
    let comm_q_mul = comms_q_mul[0].commitment().clone();
    let comm_s_id = comms_s_id[0].commitment().clone();
    let comm_s_sigma = comms_s_sigma[0].commitment().clone();

    // Absorb commitments into transcript in the correct order
    // Use raw transcript methods to avoid trait bound issues
    let mut buf = Vec::new();
    comm_q_add.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"q_add", &buf);
    
    let mut buf = Vec::new();
    comm_q_mul.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"q_mul", &buf);
    
    let mut buf = Vec::new();
    comm_a.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"a", &buf);
    
    let mut buf = Vec::new();
    comm_b.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"b", &buf);
    
    let mut buf = Vec::new();
    comm_c.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"c", &buf);
    
    let mut buf = Vec::new();
    comm_s_id.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"s_id", &buf);
    
    let mut buf = Vec::new();
    comm_s_sigma.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"s_sigma", &buf);

    // Derive beta and gamma challenges for permutation argument
    let beta = transcript.challenge_beta::<Fr>();
    let gamma = transcript.challenge_gamma::<Fr>();

    println!("=== Fiat-Shamir Challenges ===");
    println!("beta: {}", beta);
    println!("gamma: {}", gamma);

    // Now compute the grand product polynomial using the derived challenges
    let s_id_vals: Vec<Fr> = (0..3 * domain.size()).map(|i| Fr::from(i as u64)).collect();
    let z_poly = crate::circuit::Circuit::<Fr>::build_grand_product(witness_flat, sigma, domain.clone(), beta, gamma, &s_id_vals);
    
    // Commit to the grand product polynomial
    let labeled_z = LabeledPolynomial::new("z".to_string(), z_poly.clone(), None, Some(1));
    let (comms_z, states_z) = PCS::commit(pk, [&labeled_z], Some(rng)).unwrap();
    let comm_z = comms_z[0].commitment().clone();
    
    // Absorb grand product commitment
    let mut buf = Vec::new();
    comm_z.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"z", &buf);

    // Derive alpha challenge for quotient polynomial
    let alpha = transcript.challenge_alpha::<Fr>();
    println!("alpha: {}", alpha);

    // Now compute the quotient polynomial using all the derived challenges
    // We need to create a temporary circuit to compute the quotient polynomial
    let mut temp_circuit = crate::circuit::Circuit::<Fr>::from_builder(
        crate::circuit::CircuitBuilder::<Fr>::new(), 
        domain.clone()
    );
    temp_circuit.witness = crate::circuit::WitnessTable {
        a_col: a.coeffs.clone(),
        b_col: b.coeffs.clone(),
        c_col: c.coeffs.clone(),
        q_add: domain.fft(&q_add.coeffs),
        q_mul: domain.fft(&q_mul.coeffs),
    };
    
    // Set the permutation argument with the derived challenges
    let z_vals = domain.fft(&z_poly.coeffs);
    temp_circuit.permutation_argument = Some(crate::circuit::PermutationArgument {
        s_id_vals,
        s_sigma_vals: (0..3 * domain.size()).map(|i| s_sigma.evaluate(&domain.element(i))).collect(),
        z_vals,
        beta,
        gamma,
        alpha,
    });
    
    // Compute the quotient polynomial
    let t_poly = temp_circuit.build_quotient_polynomial(sigma);
    println!("Quotient polynomial degree: {}", t_poly.degree());

    // Commit to quotient polynomial
    let labeled_t = LabeledPolynomial::new("t".to_string(), t_poly.clone(), None, Some(1));
    let (comms_t, states_t) = PCS::commit(pk, [&labeled_t], Some(rng)).unwrap();
    let comm_t = comms_t[0].commitment().clone();
    
    // Absorb quotient commitment
    let mut buf = Vec::new();
    comm_t.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"t", &buf);

    // Derive zeta challenge for evaluation point
    let zeta = transcript.challenge_zeta::<Fr>();
    println!("zeta: {}", zeta);

    // Open polynomials at zeta, unpacking evaluation + proof
    let eval_a = a.evaluate(&zeta);
    let eval_b = b.evaluate(&zeta);
    let eval_c = c.evaluate(&zeta);
    let eval_q_add = q_add.evaluate(&zeta);
    let eval_q_mul = q_mul.evaluate(&zeta);
    let eval_s_id = s_id.evaluate(&zeta);
    let eval_s_sigma = s_sigma.evaluate(&zeta);
    let eval_z = z_poly.evaluate(&zeta);
    let eval_t = t_poly.evaluate(&zeta);

    // Absorb evaluations into transcript
    transcript.absorb_evaluations(
        &eval_a, &eval_b, &eval_c,
        &eval_q_add, &eval_q_mul,
        &eval_s_id, &eval_s_sigma,
        &eval_z, &eval_t,
    );

    // Create opening proofs using actual MarlinKZG10
    let mut sponge_a = test_sponge::<Fr>();
    let mut sponge_b = test_sponge::<Fr>();
    let mut sponge_c = test_sponge::<Fr>();
    let mut sponge_q_add = test_sponge::<Fr>();
    let mut sponge_q_mul = test_sponge::<Fr>();
    let mut sponge_s_id = test_sponge::<Fr>();
    let mut sponge_s_sigma = test_sponge::<Fr>();
    let mut sponge_z = test_sponge::<Fr>();
    let mut sponge_t = test_sponge::<Fr>();

    let proof_a = PCS::open(&pk, [&labeled_a], &comms_a, &zeta, &mut sponge_a, &states_a, None).unwrap();
    let proof_b = PCS::open(&pk, [&labeled_b], &comms_b, &zeta, &mut sponge_b, &states_b, None).unwrap();
    let proof_c = PCS::open(&pk, [&labeled_c], &comms_c, &zeta, &mut sponge_c, &states_c, None).unwrap();
    let proof_q_add = PCS::open(&pk, [&labeled_q_add], &comms_q_add, &zeta, &mut sponge_q_add, &states_q_add, None).unwrap();
    let proof_q_mul = PCS::open(&pk, [&labeled_q_mul], &comms_q_mul, &zeta, &mut sponge_q_mul, &states_q_mul, None).unwrap();
    let proof_s_id = PCS::open(&pk, [&labeled_s_id], &comms_s_id, &zeta, &mut sponge_s_id, &states_s_id, None).unwrap();
    let proof_s_sigma = PCS::open(&pk, [&labeled_s_sigma], &comms_s_sigma, &zeta, &mut sponge_s_sigma, &states_s_sigma, None).unwrap();
    let proof_z = PCS::open(&pk, [&labeled_z], &comms_z, &zeta, &mut sponge_z, &states_z, None).unwrap();
    let proof_t = PCS::open(&pk, [&labeled_t], &comms_t, &zeta, &mut sponge_t, &states_t, None).unwrap();

    PlonkProof {
        eval_a,
        eval_b,
        eval_c,
        eval_q_add,
        eval_q_mul,
        eval_s_id,
        eval_s_sigma,
        eval_z,
        eval_t,

        comm_a,
        comm_b,
        comm_c,
        comm_q_add,
        comm_q_mul,
        comm_s_id,
        comm_s_sigma,
        comm_z,
        comm_t,

        proof_a,
        proof_b,
        proof_c,
        proof_q_add,
        proof_q_mul,
        proof_s_id,
        proof_s_sigma,
        proof_z,
        proof_t,
    }
}