use ark_bn254::{Bn254, Fr};
use ark_poly::{
    univariate::DensePolynomial,
    DenseUVPolynomial,
    Polynomial,
};
use ark_poly_commit::{
    marlin_pc::MarlinKZG10,
    LabeledPolynomial,
    PolynomialCommitment,
};
use ark_crypto_primitives::sponge::poseidon::{PoseidonSponge, PoseidonConfig};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::UniformRand;
use ark_std::test_rng;

// Type aliases for convenience
type UniPoly = DensePolynomial<Fr>;
type PCS = MarlinKZG10<Bn254, UniPoly>;

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

/// Example usage
fn main() {
    
    let rng = &mut test_rng();
    let max_degree = 16;

    // 1. Setup - generate public parameters
    let pp = PCS::setup(max_degree, None, rng).unwrap();

    // Create a test polynomial
    let degree = 10;
    let poly = UniPoly::rand(degree, rng);
    let point = Fr::rand(rng);

    // Create labeled polynomial
    let label = String::from("test_poly");
    let labeled_poly = LabeledPolynomial::new(
        label.clone(),
        poly.clone(),
        Some(degree),
        Some(1), // we will open at one point
    );

    // Create sponge for challenge generation
    let mut poseidon = test_sponge::<Fr>();

    // 2. Trim - create prover and verifier keys
    let (ck, vk) = PCS::trim(&pp, degree, 1, Some(&[degree])).unwrap();

    // 3. Commit
    let (comms, states) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();

    // 4. Open - create opening proof
    let proof = PCS::open(&ck, [&labeled_poly], &comms, &point, &mut poseidon, &states, None).unwrap();

    // 5. Check - verify the proof
    let value = poly.evaluate(&point);
    let mut poseidon2 = test_sponge::<Fr>();
    let is_valid = PCS::check(&vk, &comms, &point, [value], &proof, &mut poseidon2, Some(rng)).unwrap();

    println!("KZG proof is valid: {}", is_valid);
}
