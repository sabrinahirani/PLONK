mod circuit;
mod poly_utils;
mod prover;
mod verifier;
mod transcript;
use ark_ff::One;

use ark_bn254::{Bn254, Fr};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::marlin_pc::MarlinKZG10;
use ark_std::test_rng;
use ark_ff::{UniformRand, Field, FftField, Zero};
use circuit::{CircuitBuilder, GateType, Circuit, PermutationArgument};
use poly_utils::interpolate_permutation_polynomials;
use prover::create_plonk_proof_with_transcript;
use verifier::verify_plonk_proof_with_transcript;
use ark_poly::{DenseUVPolynomial, EvaluationDomain};
use ark_poly_commit::PolynomialCommitment;
use ark_poly::Polynomial;

type UniPoly = DensePolynomial<Fr>;
type PCS = MarlinKZG10<Bn254, UniPoly>;

fn main() {
    let rng = &mut test_rng();
    let max_degree = 32;

    // === 1. Build circuit: chain 4 gates ===
    let mut builder = CircuitBuilder::<Fr>::new();
    // Create variables only for unique logical values
    let a = builder.new_variable(Some(Fr::from(3u64)));
    let b = builder.new_variable(Some(Fr::from(4u64)));
    let c = builder.add_gate(GateType::Add, a, b); // c = a + b = 7
    let d = builder.new_variable(Some(Fr::from(2u64)));
    let e = builder.add_gate(GateType::Mul, c, d); // e = c * d = 14
    let f = builder.new_variable(Some(Fr::from(5u64)));
    let g = builder.add_gate(GateType::Add, e, f); // g = e + f = 19
    // Instead of creating a new variable for 3, reuse 'a' for the final multiplication
    let out = builder.add_gate(GateType::Mul, g, a); // out = g * a = 19 * 3 = 57
    // builder.mark_public(out);

    // 2. Determine domain size (next power of two >= number of gates)
    let num_gates = builder.gates.len();
    let domain_size = num_gates.next_power_of_two();
    let domain = ark_poly::GeneralEvaluationDomain::<Fr>::new(domain_size).unwrap();
    let permutation_domain_size = 3 * domain.size();
    let permutation_domain = ark_poly::GeneralEvaluationDomain::<Fr>::new(permutation_domain_size).unwrap();

    // 3. Build the circuit with the correct domain size
    let mut circuit = Circuit::from_builder(builder, domain.clone());
    println!("Circuit built with {} gates", circuit.builder.gates.len());
    
    // Debug: Print all gates
    for (i, gate) in circuit.builder.gates.iter().enumerate() {
        println!("Gate {}: {:?} (inputs: [{}, {}], output: {})", 
                i, gate.gate_type, gate.inputs[0].0, gate.inputs[1].0, gate.output.0);
    }

    // 4. Compute permutation layout and sigma
    let layout = &circuit.permutation;
    let mut sigma_for_grand_product = layout.compute_sigma_mapping(domain.size());
    while sigma_for_grand_product.len() < 3 * domain.size() {
        sigma_for_grand_product.push(sigma_for_grand_product.len());
    }
    println!("Permutation mapping (σ): {:?}", sigma_for_grand_product);
    println!("sigma_for_grand_product.len(): {} expected_len: {}", sigma_for_grand_product.len(), 3 * domain.size());
    for (idx, &val) in sigma_for_grand_product.iter().enumerate() {
        if val >= 3 * domain.size() {
            panic!("sigma_for_grand_product[{}] = {} is out of range for expected_len {}", idx, val, 3 * domain.size());
        }
    }
    // Print the full permutation layout for debugging
    println!("Permutation layout (variable index -> wire positions):");
    for (var, positions) in &layout.positions {
        println!("  var {}: {:?}", var, positions);
    }
    
    // Create sigma for polynomial interpolation (padded to permutation domain size)
    let mut sigma_for_interpolation = sigma_for_grand_product.clone();
    while sigma_for_interpolation.len() < permutation_domain.size() {
        sigma_for_interpolation.push(sigma_for_interpolation.len());
    }
    println!("sigma_for_interpolation.len(): {}", sigma_for_interpolation.len());
    let (s_id_poly, s_sigma_poly) = interpolate_permutation_polynomials::<Fr>(&sigma_for_interpolation, permutation_domain);
    println!("s_id_poly degree: {}", s_id_poly.degree());
    println!("s_sigma_poly degree: {}", s_sigma_poly.degree());

    // Generate s_id_vals directly: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] for 3*4=12 values
    let s_id_vals: Vec<Fr> = (0..3 * domain.size()).map(|i| Fr::from(i as u64)).collect();
    // Generate s_sigma_vals by evaluating the sigma polynomial at the domain points
    let s_sigma_vals: Vec<Fr> = (0..3 * domain.size()).map(|i| s_sigma_poly.evaluate(&domain.element(i))).collect();
    println!("s_id_vals len: {} first 5: {:?}", s_id_vals.len(), &s_id_vals[..5.min(s_id_vals.len())]);
    println!("s_sigma_vals len: {} first 5: {:?}", s_sigma_vals.len(), &s_sigma_vals[..5.min(s_sigma_vals.len())]);

    // === Flatten witness and compute grand product ===
    let mut witness_flat = circuit.witness.flatten();
    // Pad witness_flat to 3 * domain.size() (not permutation_domain.size())
    while witness_flat.len() < 3 * domain.size() {
        witness_flat.push(Fr::zero());
    }
    println!("Flattened witness length: {} first 5: {:?}", witness_flat.len(), &witness_flat[..5.min(witness_flat.len())]);

    // sigma_for_grand_product is already padded to 3 * domain.size()
    println!("sigma_for_grand_product (after padding) len: {} first 5: {:?}", sigma_for_grand_product.len(), &sigma_for_grand_product[..5.min(sigma_for_grand_product.len())]);

    // Debug domain and vector sizes
    println!("=== Domain and Vector Size Debug ===");
    println!("domain.size(): {}", domain.size());
    println!("permutation_domain.size(): {}", permutation_domain.size());
    println!("witness_flat.len(): {}", witness_flat.len());
    println!("sigma_for_grand_product.len(): {}", sigma_for_grand_product.len());
    println!("s_id_vals.len(): {}", s_id_vals.len());
    println!("s_sigma_vals.len(): {}", s_sigma_vals.len());
    println!("Expected: 3 * domain.size() = {}", 3 * domain.size());
    
    // Verify consistency
    assert_eq!(witness_flat.len(), 3 * domain.size(), "witness_flat length mismatch");
    assert_eq!(sigma_for_grand_product.len(), 3 * domain.size(), "sigma_for_grand_product length mismatch");
    assert_eq!(s_id_vals.len(), 3 * domain.size(), "s_id_vals length mismatch");
    assert_eq!(s_sigma_vals.len(), 3 * domain.size(), "s_sigma_vals length mismatch");

    // Note: The grand product polynomial (z_poly) will be computed in the prover using transcript-derived challenges

    // === Detailed wire-by-wire debug analysis ===
    println!("\n=== Wire-by-Wire Analysis ===");
    for row in 0..witness_flat.len()/3 {
        println!("Row {}:", row);
        let mut row_numerator = Fr::one();
        let mut row_denominator = Fr::one();
        
        for wire_idx in 0..3 {
            let w_val = witness_flat[3*row + wire_idx];
            let s_val = s_id_vals[3*row + wire_idx];
            let sigma_val = sigma_for_grand_product[3*row + wire_idx];
            // Note: We'll use actual challenges from transcript later
            let beta = Fr::rand(rng); // Placeholder
            let gamma = Fr::rand(rng); // Placeholder
            let num_factor = w_val + beta * s_val + gamma;
            let den_factor = w_val + beta * Fr::from(sigma_val as u64) + gamma;
            
            row_numerator = row_numerator * num_factor;
            row_denominator = row_denominator * den_factor;
            
            println!("  wire {}: w={}, s_id={}, sigma={}, num_factor={}, den_factor={}", 
                    wire_idx, w_val, s_val, sigma_val, num_factor, den_factor);
        }
        println!("  Row {}: numerator_product={}, denominator_product={}", row, row_numerator, row_denominator);
        println!("  Row {}: ratio = {}", row, row_numerator * row_denominator.inverse().unwrap());
    }

    // === Permutation consistency check ===
    println!("\n=== Permutation Consistency Check ===");
    for i in 0..witness_flat.len() {
        let permuted_idx = sigma_for_grand_product[i];
        let original_val = witness_flat[i];
        let permuted_val = witness_flat[permuted_idx];
        println!("Position {}: original_val={}, permuted_to={}, permuted_val={}, match={}", 
                i, original_val, permuted_idx, permuted_val, original_val == permuted_val);
    }

    // === Witness and sigma ordering verification ===
    println!("\n=== Witness and Sigma Ordering Verification ===");
    println!("Witness flat: {:?}", witness_flat);
    println!("Sigma mapping: {:?}", sigma_for_grand_product);
    println!("Expected ordering: [A_0, B_0, C_0, A_1, B_1, C_1, A_2, B_2, C_2, A_3, B_3, C_3]");

    // Note: We'll set the permutation argument later when we have challenges from transcript
    // For now, use placeholder values
    let z_vals_placeholder: Vec<Fr> = (0..domain.size()).map(|_| Fr::one()).collect();
    circuit.permutation_argument = Some(PermutationArgument {
        s_id_vals,
        s_sigma_vals,
        z_vals: z_vals_placeholder,
        beta: Fr::rand(rng), // Placeholder
        gamma: Fr::rand(rng), // Placeholder
        alpha: Fr::rand(rng), // Placeholder
    });

    // Note: The quotient polynomial will be computed in the prover using transcript-derived challenges

    // === 4. Setup PCS (KZG) ===
    let pp = PCS::setup(max_degree, None, rng).unwrap();
    let (pk, vk) = PCS::trim(&pp, max_degree, 1, None).unwrap();
    println!("PCS setup complete.");

    // === 5. Create polynomials for proof generation ===
    let a_poly = DensePolynomial::from_coefficients_vec(circuit.witness.a_col.clone());
    let b_poly = DensePolynomial::from_coefficients_vec(circuit.witness.b_col.clone());
    let c_poly = DensePolynomial::from_coefficients_vec(circuit.witness.c_col.clone());
    let q_add_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&circuit.witness.q_add));
    let q_mul_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&circuit.witness.q_mul));

    println!("Polynomial degrees:");
    println!("  a_poly: {}", a_poly.degree());
    println!("  b_poly: {}", b_poly.degree());
    println!("  c_poly: {}", c_poly.degree());
    println!("  q_add_poly: {}", q_add_poly.degree());
    println!("  q_mul_poly: {}", q_mul_poly.degree());
    println!("  s_id_poly: {}", s_id_poly.degree());
    println!("  s_sigma_poly: {}", s_sigma_poly.degree());

    // Get public inputs
    let public_inputs: Vec<Fr> = circuit.builder.public_inputs
        .iter()
        .map(|v| circuit.builder.variables[v.0].unwrap())
        .collect();
    println!("Public inputs: {:?}", public_inputs);

    // === 6. Prove using transcript-based proof generation ===
    println!("\n=== Generating PLONK Proof ===");
    let proof = create_plonk_proof_with_transcript(
        &pk,
        &a_poly, &b_poly, &c_poly,
        &q_add_poly, &q_mul_poly,
        &s_id_poly, &s_sigma_poly,
        &public_inputs,
        &witness_flat,
        &sigma_for_grand_product,
        &domain,
        rng,
    );
    println!("✅ Proof created using transcript-based generation.");

    // === 7. Verify using transcript-based verification ===
    println!("\n=== Verifying PLONK Proof ===");
    let is_valid = verify_plonk_proof_with_transcript(
        &vk,
        &proof,
        &public_inputs,
    );

    println!("🎉 PLONK proof verification result: {}", is_valid);

    println!("\n=== Summary ===");
    println!("Circuit gates: {}", circuit.builder.gates.len());
    println!("Domain size: {}", domain.size());
    println!("Public inputs: {}", public_inputs.len());
    println!("Proof generated: ✅");
    println!("Proof verified: {}", if is_valid { "✅" } else { "❌" });
}
