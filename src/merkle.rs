use ark_bn254::Fr;
use ark_crypto_primitives::crh::poseidon::{PoseidonCRH, PoseidonParameters};
use ark_crypto_primitives::CRH;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct MerkleTree<T: Clone + CanonicalSerialize + CanonicalDeserialize> {
    pub depth: usize,
    pub leaves: Vec<Option<T>>,
    pub hash_leaves: Vec<Fr>,
    pub root: Fr,
    pub parameters: PoseidonParameters<Fr>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MerkleProof {
    pub root: Fr,
    pub path: Vec<Fr>,
    pub index: usize,
}

impl<T: Clone + CanonicalSerialize + CanonicalDeserialize> MerkleTree<T> {
    pub fn new(depth: usize) -> Self {
        let param = ark_crypto_primitives::crh::poseidon::Parameters::<Fr>::default();
        let size = 1 << depth;
        Self {
            depth,
            leaves: vec![None; size],
            hash_leaves: vec![Fr::zero(); size],
            root: Fr::zero(),
            parameters: param,
        }
    }

    pub fn get_leaf(&self, index: usize) -> Option<T> {
        self.leaves.get(index).cloned().flatten()
    }

    pub fn update(&mut self, index: usize, value: T) -> Result<(), &'static str> {
        if index >= self.leaves.len() {
            return Err("Index out of bounds");
        }

        self.leaves[index] = Some(value.clone());

        // Hash the serialized value
        let mut bytes = vec![];
        value
            .serialize(&mut bytes)
            .map_err(|_| "Serialization failed")?;
        let leaf_hash = Fr::from_le_bytes_mod_order(&bytes);
        self.hash_leaves[index] = leaf_hash;

        self.root = self.compute_root();
        Ok(())
    }

    pub fn root(&self) -> Fr {
        self.root
    }

    pub fn compute_root(&self) -> Fr {
        let mut layer = self.hash_leaves.clone();
        let mut next = vec![];

        for _ in 0..self.depth {
            for i in (0..layer.len()).step_by(2) {
                let left = layer[i];
                let right = layer.get(i + 1).copied().unwrap_or(Fr::zero());

                let hash = PoseidonCRH::evaluate(&self.parameters, &[left, right]).unwrap();
                next.push(hash);
            }
            layer = next;
            next = vec![];
        }

        layer[0]
    }

    pub fn prove(&self, index: usize) -> Result<MerkleProof, &'static str> {
        if index >= self.hash_leaves.len() {
            return Err("Index out of bounds");
        }

        let mut path = Vec::new();
        let mut idx = index;
        let mut layer = self.hash_leaves.clone();

        for _ in 0..self.depth {
            let sibling = if idx % 2 == 0 {
                layer.get(idx + 1).copied().unwrap_or(Fr::zero())
            } else {
                layer[idx - 1]
            };
            path.push(sibling);

            // Move to parent index
            idx /= 2;

            // Compute parent layer
            let mut next = vec![];
            for i in (0..layer.len()).step_by(2) {
                let l = layer[i];
                let r = layer.get(i + 1).copied().unwrap_or(Fr::zero());
                let hash = PoseidonCRH::evaluate(&self.parameters, &[l, r]).unwrap();
                next.push(hash);
            }
            layer = next;
        }

        Ok(MerkleProof {
            root: self.root,
            path,
            index,
        })
    }
}
