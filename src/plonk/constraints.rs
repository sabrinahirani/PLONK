use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Radix2EvaluationDomain};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::Field;

fn zero_test_check(
    poly: DensePolynomial<Fr>,
    domain: Radix2EvaluationDomain<Fr>,
) -> bool {
    let z_h = domain.vanishing_polynomial(); // Z_H(X) = X^n - 1
    let (quotient, remainder) = DensePolynomial::divide_with_q_and_r(&poly, &z_h).unwrap();
    remainder.is_zero()
}

fn compute_grand_product(
    a_eval: &[Fr],
    b_eval: &[Fr],
    c_eval: &[Fr],
    s1_eval: &[Fr],
    s2_eval: &[Fr],
    s3_eval: &[Fr],
    beta: Fr,
    gamma: Fr,
) -> Vec<Fr> {
    let n = a_eval.len();
    let mut z = vec![Fr::one(); n + 1];

    for i in 0..n {
        let num = (a_eval[i] + beta * s1_eval[i] + gamma)
            * (b_eval[i] + beta * s2_eval[i] + gamma)
            * (c_eval[i] + beta * s3_eval[i] + gamma);

        let denom = (a_eval[i] + beta * Fr::from(i as u64) + gamma)
            * (b_eval[i] + beta * Fr::from(i as u64) + gamma)
            * (c_eval[i] + beta * Fr::from(i as u64) + gamma);

        z[i + 1] = z[i] * (num * denom.inverse().unwrap());
    }

    z
}

/// Constructs the permutation constraint polynomial to include in the quotient.
fn permutation_check_polynomial(
    a: &DensePolynomial<Fr>,
    b: &DensePolynomial<Fr>,
    c: &DensePolynomial<Fr>,
    s1: &DensePolynomial<Fr>,
    s2: &DensePolynomial<Fr>,
    s3: &DensePolynomial<Fr>,
    z: &DensePolynomial<Fr>,
    domain: Radix2EvaluationDomain<Fr>,
    beta: Fr,
    gamma: Fr,
    alpha: Fr,
) -> DensePolynomial<Fr> {
    // Evaluate all polynomials over the domain
    let a_evals = domain.fft(a);
    let b_evals = domain.fft(b);
    let c_evals = domain.fft(c);
    let s1_evals = domain.fft(s1);
    let s2_evals = domain.fft(s2);
    let s3_evals = domain.fft(s3);
    let z_evals = domain.fft(z);

    let mut permutation_quotient_evals = vec![Fr::zero(); domain.size()];
    let n = domain.size();

    for i in 0..n {
        let id = Fr::from(i as u64);

        // left side (Z(x))
        let lhs = z_evals[i]
            * (a_evals[i] + beta * id + gamma)
            * (b_evals[i] + beta * id + gamma)
            * (c_evals[i] + beta * id + gamma);

        // right side (Z(ωx))
        let omega_idx = (i + 1) % n;
        let rhs = z_evals[omega_idx]
            * (a_evals[i] + beta * s1_evals[i] + gamma)
            * (b_evals[i] + beta * s2_evals[i] + gamma)
            * (c_evals[i] + beta * s3_evals[i] + gamma);

        // Store the difference: α · (LHS - RHS)
        permutation_quotient_evals[i] = alpha * (lhs - rhs);
    }

    // Enforce Z(1) = 1 using L_0(x)
    let l_0 = domain.evaluate_lagrange_basis(0); // L_0(x) over domain
    for i in 0..n {
        permutation_quotient_evals[i] += alpha * l_0[i] * (z_evals[i] - Fr::one());
    }

    // Interpolate to get the actual polynomial in coefficient form
    DensePolynomial::from_coefficients_vec(domain.ifft(&permutation_quotient_evals))
}
