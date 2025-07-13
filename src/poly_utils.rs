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

/// Interpolates the identity and sigma permutation polynomials over the given domain.
pub fn interpolate_permutation_polynomials<F: Field + ark_ff::FftField>(
    sigma: &[usize],
    domain: GeneralEvaluationDomain<F>,
) -> (DensePolynomial<F>, DensePolynomial<F>) {
    let mut sigma_padded = sigma.to_vec();
    while sigma_padded.len() < domain.size() {
        sigma_padded.push(sigma_padded.len());
    }
    assert_eq!(sigma_padded.len(), domain.size(), "Sigma length must match domain size");

    let s_id: Vec<F> = domain.elements().collect();
    let s_sigma: Vec<F> = sigma_padded.iter().map(|&i| domain.element(i)).collect();

    let s_id_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&s_id));
    let s_sigma_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&s_sigma));

    (s_id_poly, s_sigma_poly)
}