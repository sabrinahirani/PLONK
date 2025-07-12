use ark_poly::{
    univariate::DensePolynomial,
    EvaluationDomain,
    GeneralEvaluationDomain,
    DenseUVPolynomial, // ⬅️ This is the missing trait!
};

use ark_ff::{Field};
use ark_std::rand::Rng;

pub fn interpolate_selector<F: Field + ark_ff::FftField>(values: &[F]) -> DensePolynomial<F> {
    let domain = GeneralEvaluationDomain::<F>::new(values.len()).unwrap();
    let coeffs = domain.ifft(values);
    DensePolynomial::from_coefficients_vec(coeffs)
}

pub fn blind_poly<F: Field, R: Rng>(
    poly: &DensePolynomial<F>,
    degree: usize,
    rng: &mut R,
) -> DensePolynomial<F> {
    let blinded = poly.clone();

    let r1 = F::rand(rng);
    let r2 = F::rand(rng);

    // Extend coeffs to degree n+2 if needed
    let mut coeffs = blinded.coeffs.clone();
    if coeffs.len() <= degree {
        coeffs.resize(degree + 2, F::zero());
    }

    coeffs[degree] += r1;
    coeffs[degree + 1] += r2;

    DensePolynomial::from_coefficients_vec(coeffs)
}

