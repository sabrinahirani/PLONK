use ark_ff::{PrimeField};
use merlin::Transcript;

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
