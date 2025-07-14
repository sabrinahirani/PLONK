use ark_ff::{PrimeField};
use merlin::Transcript;
use ark_ec::AffineRepr;
use ark_serialize::CanonicalSerialize;

/// Helper to absorb field elements into the transcript
pub fn absorb_field<F: PrimeField>(transcript: &mut Transcript, label: &'static [u8], x: &F) {
    let mut buf = Vec::new();
    x.serialize_compressed(&mut buf).unwrap();
    transcript.append_message(label, &buf);
}

/// Extract a challenge from the transcript
pub fn challenge_scalar<F: PrimeField>(transcript: &mut Transcript, label: &'static [u8]) -> F {
    let mut buf = [0u8; 64];
    transcript.challenge_bytes(label, &mut buf);
    F::from_le_bytes_mod_order(&buf)
}

/// Helper to absorb group elements (commitments) into the transcript
pub fn absorb_commitment<G: AffineRepr>(transcript: &mut Transcript, label: &'static [u8], comm: &G) {
    let mut buf = Vec::new();
    comm.serialize_compressed(&mut buf).unwrap();
    transcript.append_message(label, &buf);
}

/// Helper to absorb any serializable commitment into the transcript
pub fn absorb_serializable_commitment<C: CanonicalSerialize>(
    transcript: &mut Transcript, 
    label: &'static [u8], 
    comm: &C
) {
    let mut buf = Vec::new();
    comm.serialize_compressed(&mut buf).unwrap();
    transcript.append_message(label, &buf);
}

/// Helper to absorb multiple field elements into the transcript
pub fn absorb_field_elements<F: PrimeField>(transcript: &mut Transcript, label: &'static [u8], elements: &[F]) {
    let mut buf = Vec::new();
    for element in elements {
        element.serialize_compressed(&mut buf).unwrap();
    }
    transcript.append_message(label, &buf);
}

/// Helper to absorb public inputs into the transcript
pub fn absorb_public_inputs<F: PrimeField>(transcript: &mut Transcript, public_inputs: &[F]) {
    if !public_inputs.is_empty() {
        absorb_field_elements(transcript, b"public_inputs", public_inputs);
    }
}

/// PLONK-specific transcript operations
pub struct PlonkTranscript {
    transcript: Transcript,
}

impl PlonkTranscript {
    /// Create a new PLONK transcript with a domain separator
    pub fn new(domain_separator: &'static [u8]) -> Self {
        Self {
            transcript: Transcript::new(domain_separator),
        }
    }

    /// Absorb selector polynomial commitments
    pub fn absorb_selector_commitments<C: CanonicalSerialize>(
        &mut self,
        comm_q_add: &C,
        comm_q_mul: &C,
    ) {
        absorb_serializable_commitment(&mut self.transcript, b"q_add", comm_q_add);
        absorb_serializable_commitment(&mut self.transcript, b"q_mul", comm_q_mul);
    }

    /// Absorb witness polynomial commitments
    pub fn absorb_witness_commitments<C: CanonicalSerialize>(
        &mut self,
        comm_a: &C,
        comm_b: &C,
        comm_c: &C,
    ) {
        absorb_serializable_commitment(&mut self.transcript, b"a", comm_a);
        absorb_serializable_commitment(&mut self.transcript, b"b", comm_b);
        absorb_serializable_commitment(&mut self.transcript, b"c", comm_c);
    }

    /// Absorb permutation polynomial commitments
    pub fn absorb_permutation_commitments<C: CanonicalSerialize>(
        &mut self,
        comm_s_id: &C,
        comm_s_sigma: &C,
    ) {
        absorb_serializable_commitment(&mut self.transcript, b"s_id", comm_s_id);
        absorb_serializable_commitment(&mut self.transcript, b"s_sigma", comm_s_sigma);
    }

    /// Absorb grand product polynomial commitment
    pub fn absorb_grand_product_commitment<C: CanonicalSerialize>(&mut self, comm_z: &C) {
        absorb_serializable_commitment(&mut self.transcript, b"z", comm_z);
    }

    /// Absorb quotient polynomial commitment
    pub fn absorb_quotient_commitment<C: CanonicalSerialize>(&mut self, comm_t: &C) {
        absorb_serializable_commitment(&mut self.transcript, b"t", comm_t);
    }

    /// Derive beta challenge for permutation argument
    pub fn challenge_beta<F: PrimeField>(&mut self) -> F {
        challenge_scalar::<F>(&mut self.transcript, b"beta")
    }

    /// Derive gamma challenge for permutation argument
    pub fn challenge_gamma<F: PrimeField>(&mut self) -> F {
        challenge_scalar::<F>(&mut self.transcript, b"gamma")
    }

    /// Derive alpha challenge for quotient polynomial
    pub fn challenge_alpha<F: PrimeField>(&mut self) -> F {
        challenge_scalar::<F>(&mut self.transcript, b"alpha")
    }

    /// Derive zeta challenge for evaluation point
    pub fn challenge_zeta<F: PrimeField>(&mut self) -> F {
        challenge_scalar::<F>(&mut self.transcript, b"zeta")
    }

    /// Absorb evaluations at zeta
    pub fn absorb_evaluations<F: PrimeField>(
        &mut self,
        eval_a: &F,
        eval_b: &F,
        eval_c: &F,
        eval_q_add: &F,
        eval_q_mul: &F,
        eval_s_id: &F,
        eval_s_sigma: &F,
        eval_z: &F,
        eval_t: &F,
    ) {
        absorb_field(&mut self.transcript, b"eval_a", eval_a);
        absorb_field(&mut self.transcript, b"eval_b", eval_b);
        absorb_field(&mut self.transcript, b"eval_c", eval_c);
        absorb_field(&mut self.transcript, b"eval_q_add", eval_q_add);
        absorb_field(&mut self.transcript, b"eval_q_mul", eval_q_mul);
        absorb_field(&mut self.transcript, b"eval_s_id", eval_s_id);
        absorb_field(&mut self.transcript, b"eval_s_sigma", eval_s_sigma);
        absorb_field(&mut self.transcript, b"eval_z", eval_z);
        absorb_field(&mut self.transcript, b"eval_t", eval_t);
    }

    /// Get the underlying transcript for custom operations
    pub fn inner(&mut self) -> &mut Transcript {
        &mut self.transcript
    }
}
