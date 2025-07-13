mod circuit;
mod poly_utils;
mod prover;
mod verifier;
mod transcript;

use ark_bn254::{Bn254, Fr};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::marlin_pc::MarlinKZG10;
use ark_std::test_rng;
use ark_ff::{UniformRand, Field, FftField};
use circuit::{CircuitBuilder, GateType, Circuit, PermutationArgument};
use poly_utils::{interpolate_permutation_polynomials, compute_grand_product};
use prover::create_plonk_proof;
use verifier::verify_plonk_proof;
use ark_poly::{DenseUVPolynomial, EvaluationDomain};
use ark_poly_commit::PolynomialCommitment;
use ark_poly::Polynomial;

type UniPoly = DensePolynomial<Fr>;
type PCS = MarlinKZG10<Bn254, UniPoly>;

fn main() {
    let rng = &mut test_rng();
    let max_degree = 32;

    // === 1. Build circuit: c = a + b ===
    let mut builder = CircuitBuilder::<Fr>::new();
    let a = builder.new_variable(Some(Fr::from(3u64)));
    let b = builder.new_variable(Some(Fr::from(4u64)));
    let c = builder.add_gate(GateType::Add, a, b);
    builder.mark_public(c);

    let mut circuit = Circuit::from_builder(builder);
    println!("Circuit built with {} gates", circuit.builder.gates.len());

    // === 2. Setup permutation argument ===
    let layout = &circuit.permutation;
    let sigma = layout.compute_sigma_mapping(circuit.witness.a_col.len());
    println!("Permutation mapping (σ): {:?}", sigma);

    let (s_id_poly, s_sigma_poly, domain) = interpolate_permutation_polynomials::<Fr>(&sigma);
    println!("s_id_poly degree: {}", s_id_poly.degree());
    println!("s_sigma_poly degree: {}", s_sigma_poly.degree());

    let s_id_vals = domain.fft(&s_id_poly.coeffs);
    let s_sigma_vals = domain.fft(&s_sigma_poly.coeffs);

    // === Sample challenge values ===
    let beta = Fr::rand(rng);
    let gamma = Fr::rand(rng);
    let alpha = Fr::rand(rng);

    println!("Sampled challenges:");
    println!("  β (beta):   {}", beta);
    println!("  γ (gamma):  {}", gamma);
    println!("  α (alpha):  {}", alpha);

    // === Flatten witness and compute grand product ===
    let witness_flat = circuit.witness.flatten();
    println!("Flattened witness length: {}", witness_flat.len());

    let z_poly = compute_grand_product(&witness_flat, &sigma, domain, beta, gamma);
    println!("z_poly degree: {}", z_poly.degree());

    let z_vals = domain.fft(&z_poly.coeffs);
    println!("First 5 z values: {:?}", &z_vals[..5.min(z_vals.len())]);

    circuit.permutation_argument = Some(PermutationArgument {
        s_id_vals,
        s_sigma_vals,
        z_vals,
        beta,
        gamma,
        alpha,
    });

    // === 3. Build quotient polynomial ===
    let t_poly = circuit.build_quotient_polynomial(domain);
    println!("Quotient polynomial degree: {}", t_poly.degree());

    let zeta = Fr::rand(rng);
    let z_h_eval = domain.evaluate_vanishing_polynomial(zeta);
    println!("zeta: {}", zeta);
    println!("Z_H(zeta): {}", z_h_eval);

    // === 4. Setup PCS (KZG) ===
    let pp = PCS::setup(max_degree, None, rng).unwrap();
    let (pk, vk) = PCS::trim(&pp, max_degree, 1, None).unwrap();
    println!("PCS setup complete.");

    // === 5. Prove ===
    let a_poly = DensePolynomial::from_coefficients_vec(circuit.witness.a_col.clone());
    let b_poly = DensePolynomial::from_coefficients_vec(circuit.witness.b_col.clone());
    let c_poly = DensePolynomial::from_coefficients_vec(circuit.witness.c_col.clone());
    let q_add_poly = circuit.witness.q_add_poly.clone();
    let q_mul_poly = circuit.witness.q_mul_poly.clone();

    println!("a_poly: {:?}", a_poly);
    println!("b_poly: {:?}", b_poly);
    println!("c_poly: {:?}", c_poly);
    println!("q_add_poly: {:?}", q_add_poly);
    println!("q_mul_poly: {:?}", q_mul_poly);

    let proof = create_plonk_proof(
        &pk,
        &a_poly, &b_poly, &c_poly,
        &q_add_poly, &q_mul_poly,
        &s_id_poly, &s_sigma_poly,
        &z_poly, &t_poly,
        zeta,
        rng,
    );
    println!("Proof created.");

    // === 6. Verify ===
    let public_inputs: Vec<Fr> = circuit.builder.public_inputs
        .iter()
        .map(|v| circuit.builder.variables[v.0].unwrap())
        .collect();
    println!("Public inputs: {:?}", public_inputs);

    let is_valid = verify_plonk_proof(
        &vk,
        &proof,
        &public_inputs,
        beta, gamma, alpha,
        zeta,
        z_h_eval,
    );

    println!("PLONK proof is valid: {}", is_valid);
}
