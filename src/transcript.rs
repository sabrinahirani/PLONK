use ark_bn254::Fr;
use ark_ff::PrimeField;
use ark_ff::to_bytes;
use merlin::Transcript;

/// A wrapper around a Merlin transcript for PLONK Fiat-Shamir challenges
pub struct PlonkTranscript {
    transcript: Transcript,
}

impl PlonkTranscript {
    /// Creates a new transcript with a given label (e.g., b"plonk")
    pub fn new(label: &[u8]) -> Self {
        Self {
            transcript: Transcript::new(label),
        }
    }

    /// Append a message (e.g., a commitment or challenge)
    pub fn append_message(&mut self, label: &'static [u8], msg: &[u8]) {
        self.transcript.append_message(label, msg);
    }

    /// Append a field element to the transcript
    pub fn append_fr(&mut self, label: &'static [u8], fr: &Fr) {
        let bytes = to_bytes![fr].expect("Failed to serialize Fr");
        self.append_message(label, &bytes);
    }

    /// Generate a challenge scalar from the transcript, domain-separated by label
    pub fn get_challenge(&mut self, label: &'static [u8]) -> Fr {
        let mut buf = [0u8; 64];
        self.transcript.challenge_bytes(label, &mut buf);
        Fr::from_le_bytes_mod_order(&buf)
    }
}
