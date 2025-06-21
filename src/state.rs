use ark_bn254::Fr;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use serde::{Deserialize, Serialize};

use crate::merkle::{MerkleTree, MerkleProof};

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize)]
pub struct Account {
    pub address: Fr,
    pub balance: Fr,
    pub nonce: u64,
}

impl Account {
    pub fn empty() -> Self {
        Self {
            address: Fr::zero(),
            balance: Fr::zero(),
            nonce: 0,
        }
    }
}

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize)]
pub struct Transaction {
    pub from_index: usize,
    pub to_index: usize,
    pub amount: Fr,
    pub nonce: u64,
}

#[derive(Debug)]
pub struct State {
    pub tree: MerkleTree<Account>,
}

impl State {
    pub fn new(depth: usize) -> Self {
        Self {
            tree: MerkleTree::new(depth),
        }
    }

    pub fn get_account(&self, index: usize) -> Account {
        self.tree.get_leaf(index).unwrap_or(Account::empty())
    }

    pub fn apply_tx(&mut self, tx: Transaction) -> Result<(Fr, Fr, MerkleProof, MerkleProof), &'static str> {
        let mut from = self.get_account(tx.from_index);
        let mut to = self.get_account(tx.to_index);

        if from.balance < tx.amount {
            return Err("Insufficient balance");
        }

        if from.nonce != tx.nonce {
            return Err("Invalid nonce");
        }

        // Create proofs before update
        let from_proof = self.tree.prove(tx.from_index)?;
        let to_proof = self.tree.prove(tx.to_index)?;

        // Update accounts
        from.balance -= tx.amount;
        from.nonce += 1;
        to.balance += tx.amount;

        self.tree.update(tx.from_index, from.clone())?;
        self.tree.update(tx.to_index, to.clone())?;

        Ok((
            self.tree.root(),
            from_proof.root,
            from_proof,
            to_proof
        ))
    }
}
