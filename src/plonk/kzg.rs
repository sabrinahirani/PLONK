use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ff::UniformRand;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use rand::thread_rng;

// reference: https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf
// reference: https://youtu.be/WyT5KkKBJUw?si=2-2WM20qGsS3O3qP

/// Generate the SRS (trusted setup)
fn setup(degree: usize) -> (Vec<G1Projective>, G2Projective) {
    let mut rng = thread_rng();
    let s = Fr::rand(&mut rng);
    let g1_gen = <Bls12_381 as Pairing>::G1::generator();
    let g2_gen = <Bls12_381 as Pairing>::G2::generator();

    // Powers of s in G1: [g1, s*g1, s^2*g1, ..., s^d*g1]
    let powers_of_s_g1: Vec<_> = (0..=degree)
        .map(|i| g1_gen.mul(s.pow([i as u64])))
        .collect();

    // s in G2: s * g2
    let s_g2 = g2_gen.mul(s);

    (powers_of_s_g1, s_g2)
}

/// Commit to the polynomial using G1 powers
fn commit(poly: &DensePolynomial<Fr>, srs_g1: &[G1Projective]) -> G1Projective {
    poly.coeffs
        .iter()
        .zip(srs_g1)
        .map(|(coeff, g1)| g1.mul(*coeff))
        .sum()
}

/// Open the commitment at a point z
fn open(
    poly: &DensePolynomial<Fr>,
    z: Fr,
    srs_g1: &[G1Projective],
) -> (Fr, G1Projective) {
    let y = poly.evaluate(&z);
    let divisor = DensePolynomial::from_coefficients_vec(vec![-z, Fr::one()]);
    let quotient = poly / &divisor;
    let proof = commit(&quotient, &srs_g1[..=quotient.degree()]);
    (y, proof)
}

/// Verify the opening using pairings
fn verify(
    commitment: G1Projective,
    z: Fr,
    y: Fr,
    proof: G1Projective,
    srs_g2: G2Projective,
) -> bool {
    let g2_gen = <Bls12_381 as Pairing>::G2::generator();
    let lhs = Bls12_381::pairing(commitment - G1Projective::generator().mul(y), g2_gen);
    let rhs = Bls12_381::pairing(proof, srs_g2 - g2_gen.mul(z));
    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::One;

    #[test]
    fn test_kzg() {
        let degree = 4;
        let (srs_g1, srs_g2) = setup(degree);

        // Define f(x) = 3 + 2x + x^2
        let poly = DensePolynomial::from_coefficients_vec(vec![
            Fr::from(3u64),
            Fr::from(2u64),
            Fr::from(1u64),
        ]);

        let commitment = commit(&poly, &srs_g1);

        let z = Fr::from(10u64); // Evaluation point
        let (y, proof) = open(&poly, z, &srs_g1);

        assert!(verify(commitment, z, y, proof, srs_g2));
    }
}