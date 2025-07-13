use ark_ff::{Field};
use ark_std::rand::Rng;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, DenseUVPolynomial};

/// Performs interpolation to generate selector polynomial.
pub fn interpolate_selector<F: Field + ark_ff::FftField>(values: &[F]) -> DensePolynomial<F> {
    let domain = GeneralEvaluationDomain::<F>::new(values.len()).unwrap();
    let coeffs = domain.ifft(values);
    DensePolynomial::from_coefficients_vec(coeffs)
}

/// TODO
pub fn interpolate_permutation_polynomials<F: Field + ark_ff::FftField>(
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


/// TODO
pub fn compute_grand_product<F: Field + ark_ff::FftField>(
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





