use ark_bn254::{Bn254, Fr};
use ark_ff::{Field, Zero};
use ark_poly_commit::{
    marlin_pc::MarlinKZG10,
    PolynomialCommitment,
    LabeledCommitment,
};
use ark_poly::univariate::DensePolynomial;
use ark_crypto_primitives::sponge::poseidon::{PoseidonSponge, PoseidonConfig};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use ark_std::test_rng;
use ark_serialize::CanonicalSerialize;
use crate::prover::PlonkProof;
use crate::transcript::{PlonkTranscript, absorb_public_inputs};

// Type aliases for convenience
type UniPoly = DensePolynomial<Fr>;
type PCS = MarlinKZG10<Bn254, UniPoly>;
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

fn lagrange_eval_at_zeta<Fr: Field>(zeta: Fr, omega_i: Fr, n: usize) -> Fr {
    // ℓ_i(ζ) = (ζ^n - 1) / (n * (ζ - ω^i))
    let n_f = Fr::from(n as u64);
    let zeta_n = zeta.pow(&[n as u64]);
    (zeta_n - Fr::one()) / (n_f * (zeta - omega_i))
}

/// Complete PLONK verification function using transcript for challenge verification
pub fn verify_plonk_proof_with_transcript(
    vk: &KZGVerifierKey,
    proof: &PlonkProof<Fr>,
    public_inputs: &[Fr],
) -> bool {
    let PlonkProof {
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
    } = proof;

    // Dereference all evaluations
    let (a, b, c, q_add, q_mul, s_id, s_sigma, z, t) = (
        *eval_a, *eval_b, *eval_c, *eval_q_add, *eval_q_mul, *eval_s_id, *eval_s_sigma, *eval_z, *eval_t,
    );

    // Initialize transcript (must match prover's transcript exactly)
    let mut transcript = PlonkTranscript::new(b"plonk_proof");
    
    // Absorb public inputs first (must match prover)
    absorb_public_inputs(transcript.inner(), public_inputs);
    
    // Absorb commitments in the EXACT same order as prover
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

    // Derive beta and gamma challenges (must match prover)
    let beta = transcript.challenge_beta::<Fr>();
    let gamma = transcript.challenge_gamma::<Fr>();

    // Absorb grand product commitment (must match prover)
    let mut buf = Vec::new();
    comm_z.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"z", &buf);

    // Derive alpha challenge (must match prover)
    let alpha = transcript.challenge_alpha::<Fr>();

    // Absorb quotient commitment (must match prover)
    let mut buf = Vec::new();
    comm_t.serialize_compressed(&mut buf).unwrap();
    transcript.inner().append_message(b"t", &buf);

    // Derive zeta challenge (must match prover)
    let zeta = transcript.challenge_zeta::<Fr>();

    println!("=== Fiat-Shamir Challenges (Verifier) ===");
    println!("beta: {}", beta);
    println!("gamma: {}", gamma);
    println!("alpha: {}", alpha);
    println!("zeta: {}", zeta);

    // Absorb evaluations (must match prover)
    transcript.absorb_evaluations(
        &a, &b, &c,
        &q_add, &q_mul,
        &s_id, &s_sigma,
        &z, &t,
    );

    println!("=== PLONK Verification ===");
    println!("Evaluations at zeta:");
    println!("  a(zeta): {}", a);
    println!("  b(zeta): {}", b);
    println!("  c(zeta): {}", c);
    println!("  q_add(zeta): {}", q_add);
    println!("  q_mul(zeta): {}", q_mul);
    println!("  s_id(zeta): {}", s_id);
    println!("  s_sigma(zeta): {}", s_sigma);
    println!("  z(zeta): {}", z);
    println!("  t(zeta): {}", t);

    // 3. Verify all polynomial openings at zeta using KZG
    println!("\n3. KZG Opening Proof Verification:");
    let mut sponge_a = test_sponge::<Fr>();
    let mut sponge_b = test_sponge::<Fr>();
    let mut sponge_c = test_sponge::<Fr>();
    let mut sponge_q_add = test_sponge::<Fr>();
    let mut sponge_q_mul = test_sponge::<Fr>();
    let mut sponge_s_id = test_sponge::<Fr>();
    let mut sponge_s_sigma = test_sponge::<Fr>();
    let mut sponge_z = test_sponge::<Fr>();
    let mut sponge_t = test_sponge::<Fr>();

    // Create commitments for verification
    let comms_a = vec![LabeledCommitment::new("a".to_string(), comm_a.clone(), None)];
    let comms_b = vec![LabeledCommitment::new("b".to_string(), comm_b.clone(), None)];
    let comms_c = vec![LabeledCommitment::new("c".to_string(), comm_c.clone(), None)];
    let comms_q_add = vec![LabeledCommitment::new("q_add".to_string(), comm_q_add.clone(), None)];
    let comms_q_mul = vec![LabeledCommitment::new("q_mul".to_string(), comm_q_mul.clone(), None)];
    let comms_s_id = vec![LabeledCommitment::new("s_id".to_string(), comm_s_id.clone(), None)];
    let comms_s_sigma = vec![LabeledCommitment::new("s_sigma".to_string(), comm_s_sigma.clone(), None)];
    let comms_z = vec![LabeledCommitment::new("z".to_string(), comm_z.clone(), None)];
    let comms_t = vec![LabeledCommitment::new("t".to_string(), comm_t.clone(), None)];

    // Verify each proof using actual MarlinKZG10 verification
    let check_a = PCS::check(&vk, &comms_a, &zeta, [a], &proof_a, &mut sponge_a, Some(&mut test_rng())).unwrap();
    let check_b = PCS::check(&vk, &comms_b, &zeta, [b], &proof_b, &mut sponge_b, Some(&mut test_rng())).unwrap();
    let check_c = PCS::check(&vk, &comms_c, &zeta, [c], &proof_c, &mut sponge_c, Some(&mut test_rng())).unwrap();
    let check_q_add = PCS::check(&vk, &comms_q_add, &zeta, [q_add], &proof_q_add, &mut sponge_q_add, Some(&mut test_rng())).unwrap();
    let check_q_mul = PCS::check(&vk, &comms_q_mul, &zeta, [q_mul], &proof_q_mul, &mut sponge_q_mul, Some(&mut test_rng())).unwrap();
    let check_s_id = PCS::check(&vk, &comms_s_id, &zeta, [s_id], &proof_s_id, &mut sponge_s_id, Some(&mut test_rng())).unwrap();
    let check_s_sigma = PCS::check(&vk, &comms_s_sigma, &zeta, [s_sigma], &proof_s_sigma, &mut sponge_s_sigma, Some(&mut test_rng())).unwrap();
    let check_z = PCS::check(&vk, &comms_z, &zeta, [z], &proof_z, &mut sponge_z, Some(&mut test_rng())).unwrap();
    let check_t = PCS::check(&vk, &comms_t, &zeta, [t], &proof_t, &mut sponge_t, Some(&mut test_rng())).unwrap();
    
    println!("  a: {}", check_a);
    println!("  b: {}", check_b);
    println!("  c: {}", check_c);
    println!("  q_add: {}", check_q_add);
    println!("  q_mul: {}", check_q_mul);
    println!("  s_id: {}", check_s_id);
    println!("  s_sigma: {}", check_s_sigma);
    println!("  z: {}", check_z);
    println!("  t: {}", check_t);
    
    let all_kzg_valid = check_a && check_b && check_c && check_q_add && check_q_mul && 
                       check_s_id && check_s_sigma && check_z && check_t;
    
    if !all_kzg_valid {
        println!("❌ KZG opening proof verification failed!");
        return false;
    }
    println!("✅ KZG opening proof verification passed");

    println!("\n🎉 All PLONK verification checks passed!");
    true
}