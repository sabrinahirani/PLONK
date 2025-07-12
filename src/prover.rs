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
    Commitment,
    Proof,
};
use ark_crypto_primitives::sponge::poseidon::{PoseidonSponge, PoseidonConfig};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use ark_std::rand::Rng;
use ark_std::test_rng;

use ark_bn254::{Bn254, Fr};

use ark_poly::DenseUVPolynomial;
use ark_poly::UVPolynomial;

// Type aliases for convenience
type UniPoly = DensePolynomial<Fr>;
type PCS = MarlinKZG10<Bn254, UniPoly, PoseidonSponge<Fr>>;
pub type KZGCommitment = Commitment<Bn254>;
pub type KZGOpeningProof = Proof<Bn254>;
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

// === Types for wiring and permutation ===

#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum WireColumn {
    A,
    B,
    C,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct WirePosition {
    pub col: WireColumn,
    pub row: usize,
}

#[derive(Debug)]
pub struct PermutationLayout {
    pub positions: HashMap<usize, Vec<WirePosition>>, // Variable index → used positions
}

// === Wire flattening and sigma mapping ===

pub fn flatten_wire_positions(num_rows: usize) -> Vec<WirePosition> {
    let mut positions = Vec::with_capacity(3 * num_rows);
    for row in 0..num_rows {
        positions.push(WirePosition { col: WireColumn::A, row });
        positions.push(WirePosition { col: WireColumn::B, row });
        positions.push(WirePosition { col: WireColumn::C, row });
    }
    positions
}

pub fn compute_sigma_mapping(layout: &PermutationLayout, num_rows: usize) -> Vec<usize> {
    let flat_positions = flatten_wire_positions(num_rows);
    let mut sigma = vec![0usize; 3 * num_rows];

    let mut position_to_index: HashMap<WirePosition, usize> = HashMap::new();
    for (i, pos) in flat_positions.iter().enumerate() {
        position_to_index.insert(*pos, i);
    }

    for (_var, uses) in layout.positions.iter() {
        let n = uses.len();
        for i in 0..n {
            let from = uses[i];
            let to = uses[(i + 1) % n];
            let from_idx = position_to_index[&from];
            let to_idx = position_to_index[&to];
            sigma[from_idx] = to_idx;
        }
    }

    sigma
}

// === Interpolate permutation polynomials ===

pub fn interpolate_permutation_polynomials<F: Field + FftField>(
    sigma: &[usize],
) -> (DensePolynomial<F>, DensePolynomial<F>, GeneralEvaluationDomain<F>) {
    let n = sigma.len();
    let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

    let s_id: Vec<F> = domain.elements().collect();
    let s_sigma: Vec<F> = sigma.iter().map(|&i| domain.element(i)).collect();

    let s_id_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&s_id));
    let s_sigma_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&s_sigma));

    (s_id_poly, s_sigma_poly, domain)
}

// === Compute Z(x) — Grand Product Polynomial ===

pub fn compute_grand_product<F: Field + FftField>(
    witness_flat: &[F],
    sigma: &[usize],
    domain: GeneralEvaluationDomain<F>,
    beta: F,
    gamma: F,
) -> DensePolynomial<F> {
    let n = witness_flat.len();
    let omega_powers: Vec<F> = domain.elements().collect();
    let s_id = &omega_powers;
    let s_sigma: Vec<F> = sigma.iter().map(|&i| omega_powers[i]).collect();

    let mut z_vals = vec![F::one(); n + 1];

    for i in 0..n {
        let num = witness_flat[i] + beta * s_id[i] + gamma;
        let den = witness_flat[i] + beta * s_sigma[i] + gamma;
        z_vals[i + 1] = z_vals[i] * num / den;
    }

    let z_vals = z_vals[..n].to_vec();
    DensePolynomial::from_coefficients_vec(domain.ifft(&z_vals))
}

// === Constraint Building ===

pub fn build_gate_constraint<F: Field>(
    a: &DensePolynomial<F>,
    b: &DensePolynomial<F>,
    c: &DensePolynomial<F>,
    q_add: &DensePolynomial<F>,
    q_mul: &DensePolynomial<F>,
) -> DensePolynomial<F> {
    let add_expr = q_add * &(a.clone() + b.clone() - c.clone());
    let mul_expr = q_mul * &((a.clone() * b.clone()) - c.clone());
    add_expr + mul_expr
}

pub fn build_permutation_constraint<F: Field + FftField>(
    a_vals: &[F],
    b_vals: &[F],
    c_vals: &[F],
    s_id_vals: &[F],
    s_sigma_vals: &[F],
    z_vals: &[F],
    beta: F,
    gamma: F,
    alpha: F,
    domain: GeneralEvaluationDomain<F>,
) -> DensePolynomial<F> {
    let n = a_vals.len();
    let mut lhs = vec![F::zero(); n];
    let mut rhs = vec![F::zero(); n];

    for i in 0..n {
        let a_term = a_vals[i] + beta * s_id_vals[i] + gamma;
        let b_term = b_vals[i] + beta * s_id_vals[i] + gamma;
        let c_term = c_vals[i] + beta * s_id_vals[i] + gamma;
        lhs[i] = z_vals[i] * a_term * b_term * c_term;

        let a_perm = a_vals[i] + beta * s_sigma_vals[i] + gamma;
        let b_perm = b_vals[i] + beta * s_sigma_vals[i] + gamma;
        let c_perm = c_vals[i] + beta * s_sigma_vals[i] + gamma;
        rhs[i] = z_vals[(i + 1) % n] * a_perm * b_perm * c_perm;
    }

    let evals: Vec<F> = lhs.iter().zip(rhs.iter()).map(|(l, r)| alpha * (*l - *r)).collect();
    DensePolynomial::from_coefficients_vec(domain.ifft(&evals))
}

pub fn build_public_input_poly<F: Field + FftField>(
    a_vals: &[F],
    public_inputs: &[F],
    alpha: F,
    domain: GeneralEvaluationDomain<F>,
) -> DensePolynomial<F> {
    let mut constraint = vec![F::zero(); a_vals.len()];
    for (i, pi) in public_inputs.iter().enumerate() {
        constraint[i] = alpha * (a_vals[i] - *pi);
    }
    DensePolynomial::from_coefficients_vec(domain.ifft(&constraint))
}

// === Quotient Polynomial ===

pub fn build_quotient_polynomial<F: FftField>(
    a_poly: &DensePolynomial<F>,
    b_poly: &DensePolynomial<F>,
    c_poly: &DensePolynomial<F>,
    q_add: &DensePolynomial<F>,
    q_mul: &DensePolynomial<F>,
    s_id_vals: &[F],
    s_sigma_vals: &[F],
    z_vals: &[F],
    a_vals: &[F],
    b_vals: &[F],
    c_vals: &[F],
    public_inputs: &[F],
    beta: F,
    gamma: F,
    alpha: F,
    domain: GeneralEvaluationDomain<F>,
) -> DensePolynomial<F> {
    let gate_poly = build_gate_constraint(a_poly, b_poly, c_poly, q_add, q_mul);
    let perm_poly = build_permutation_constraint(
        a_vals, b_vals, c_vals,
        s_id_vals, s_sigma_vals,
        z_vals, beta, gamma, alpha, domain,
    );
    let pub_poly = build_public_input_poly(a_vals, public_inputs, alpha, domain);

    let t_num = &gate_poly + &perm_poly + pub_poly.clone();

    let z_h = domain.vanishing_polynomial().into();

    let (t_quot, rem) = t_num.divide_with_q_and_r(&z_h).unwrap();
    assert!(rem.is_zero(), "Quotient division remainder not zero");

    t_quot
}

// === Evaluation and Proof Struct ===

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

pub fn create_plonk_proof<R: Rng>(
    pk: &KZGProverKey,
    a: &DensePolynomial<Fr>,
    b: &DensePolynomial<Fr>,
    c: &DensePolynomial<Fr>,
    q_add: &DensePolynomial<Fr>,
    q_mul: &DensePolynomial<Fr>,
    s_id: &DensePolynomial<Fr>,
    s_sigma: &DensePolynomial<Fr>,
    z: &DensePolynomial<Fr>,
    t: &DensePolynomial<Fr>,
    zeta: Fr,
    rng: &mut R,
) -> PlonkProof<Fr> {
    // Create labeled polynomials
    let labeled_a = LabeledPolynomial::new("a".to_string(), a.clone(), None, Some(1));
    let labeled_b = LabeledPolynomial::new("b".to_string(), b.clone(), None, Some(1));
    let labeled_c = LabeledPolynomial::new("c".to_string(), c.clone(), None, Some(1));
    let labeled_q_add = LabeledPolynomial::new("q_add".to_string(), q_add.clone(), None, Some(1));
    let labeled_q_mul = LabeledPolynomial::new("q_mul".to_string(), q_mul.clone(), None, Some(1));
    let labeled_s_id = LabeledPolynomial::new("s_id".to_string(), s_id.clone(), None, Some(1));
    let labeled_s_sigma = LabeledPolynomial::new("s_sigma".to_string(), s_sigma.clone(), None, Some(1));
    let labeled_z = LabeledPolynomial::new("z".to_string(), z.clone(), None, Some(1));
    let labeled_t = LabeledPolynomial::new("t".to_string(), t.clone(), None, Some(1));

    // Commit to polynomials
    let (comms_a, states_a) = PCS::commit(pk, [&labeled_a], Some(rng)).unwrap();
    let (comms_b, states_b) = PCS::commit(pk, [&labeled_b], Some(rng)).unwrap();
    let (comms_c, states_c) = PCS::commit(pk, [&labeled_c], Some(rng)).unwrap();
    let (comms_q_add, states_q_add) = PCS::commit(pk, [&labeled_q_add], Some(rng)).unwrap();
    let (comms_q_mul, states_q_mul) = PCS::commit(pk, [&labeled_q_mul], Some(rng)).unwrap();
    let (comms_s_id, states_s_id) = PCS::commit(pk, [&labeled_s_id], Some(rng)).unwrap();
    let (comms_s_sigma, states_s_sigma) = PCS::commit(pk, [&labeled_s_sigma], Some(rng)).unwrap();
    let (comms_z, states_z) = PCS::commit(pk, [&labeled_z], Some(rng)).unwrap();
    let (comms_t, states_t) = PCS::commit(pk, [&labeled_t], Some(rng)).unwrap();

    let comm_a = comms_a[0].clone();
    let comm_b = comms_b[0].clone();
    let comm_c = comms_c[0].clone();
    let comm_q_add = comms_q_add[0].clone();
    let comm_q_mul = comms_q_mul[0].clone();
    let comm_s_id = comms_s_id[0].clone();
    let comm_s_sigma = comms_s_sigma[0].clone();
    let comm_z = comms_z[0].clone();
    let comm_t = comms_t[0].clone();

    // Open polynomials at zeta, unpacking evaluation + proof
    let eval_a = a.evaluate(&zeta);
    let eval_b = b.evaluate(&zeta);
    let eval_c = c.evaluate(&zeta);
    let eval_q_add = q_add.evaluate(&zeta);
    let eval_q_mul = q_mul.evaluate(&zeta);
    let eval_s_id = s_id.evaluate(&zeta);
    let eval_s_sigma = s_sigma.evaluate(&zeta);
    let eval_z = z.evaluate(&zeta);
    let eval_t = t.evaluate(&zeta);

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
