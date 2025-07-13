use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly::DenseUVPolynomial;

use std::collections::HashMap;

use crate::poly_utils::interpolate_selector;


// reference: https://vitalik.eth.limo/general/2016/12/10/qap.html


// Circuit construction (arithmetization) for PLONK.
//
// This module defines the logic for building arithmetic circuits used in PLONK.
// It currently only supports + and *.
//
// Key components:
// - `Variable`: Represents a variable in the circuit.
// - `Gate`: Represents a constraint in the circuit.
// - `WitnessTable`: Contains witness values (and selector polynomials) for all constraints in the circuit.
// - `PermutationLayout`: Represents the wiring of the circuit.
// - `CircuitBuilder`: Manages the construction of the circuit.



/// A variable in the circuit (identified by an index).
#[derive(Copy, Clone, Debug)]
pub struct Variable(pub usize);

// types of gates supported in the circuit: + or *.
#[derive(Debug)]
pub enum GateType {
    Add,
    Mul
}

/// A gate (constraint) in the circuit.
#[derive(Debug)]
pub struct Gate<F: Field> {
    pub gate_type: GateType,    // either + or *
    pub inputs: [Variable; 2],  // input variables (2)
    pub output: Variable,       // output variable (1)
    pub selector_row: usize,    // input into selector polynomial
    pub constant: F,            // to support gates with constants (not currently supported)
}



// ---



// corresponds to A, B, and C in constraints:
// q_add * (A + B - C) + q_mul * (A * B - C) = 0
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum WireColumn {
    A,
    B,
    C,
}

/// Represents a position in the circuit.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct WirePosition {
    pub col: WireColumn,
    pub row: usize,
}

/// Represents the wiring of the circuit.
#[derive(Debug)]
pub struct PermutationLayout {
    pub positions: HashMap<usize, Vec<WirePosition>>, // variable index → used positions
}

impl PermutationLayout {
    // TODO
    /// Computes the sigma permutation mapping for the permutation argument.
    pub fn compute_sigma_mapping(&self, num_rows: usize) -> Vec<usize> {

        // constructs all possible wire positions in the circuit in linearized order
        let mut flat_positions = Vec::with_capacity(3 * num_rows);
        for row in 0..num_rows {
            flat_positions.push(WirePosition { col: WireColumn::A, row });
            flat_positions.push(WirePosition { col: WireColumn::B, row });
            flat_positions.push(WirePosition { col: WireColumn::C, row });
        }

        // creates a lookup from wire positions to indices in the flattened vector
        let mut position_to_index: HashMap<WirePosition, usize> = HashMap::new();
        for (i, pos) in flat_positions.iter().enumerate() {
            position_to_index.insert(*pos, i);
        }

        // update sigma by cycling through wire positions
        let mut sigma = vec![0usize; 3 * num_rows];
        for (_var, uses) in self.positions.iter() {
            let n = uses.len();
            for i in 0..n {
                let from = uses[i];
                let to = uses[(i + 1) % n];
                let from_idx = position_to_index[&from];
                let to_idx = position_to_index[&to];
                sigma[from_idx] = to_idx;
            }
        }

        sigma
    }
}



// --- 



/// Contains witness values (and selector polynomials) for all constraints in the circuit.
pub struct WitnessTable<F: Field> {
    pub a_col: Vec<F>,
    pub b_col: Vec<F>,
    pub c_col: Vec<F>,
    pub q_add: Vec<F>,
    pub q_mul: Vec<F>,
    pub q_add_poly: DensePolynomial<F>, // selector polynomial for +
    pub q_mul_poly: DensePolynomial<F>, // selector polynomial for *
}

impl<F: Field> WitnessTable<F> {

    /// Flattens into linearized form: [A_0, B_0, C_0, A_1, B_1, C_1, ..., A_n, B_n, C_n]
    pub fn flatten(&self) -> Vec<F> {
        let mut flat = Vec::new();
        for i in 0..self.a_col.len() {
            flat.push(self.a_col[i]);
            flat.push(self.b_col[i]);
            flat.push(self.c_col[i]);
        }
        flat
    }
}



// ---



/// Manages the construction of the circuit.
pub struct CircuitBuilder<F: Field> {
    pub next_var_index: usize,          // tracks next available index for a new variable
    pub variables: Vec<Option<F>>,      // stores value assigned to each variable (none until witness generation)
    pub public_inputs: Vec<Variable>,   // marks which variables are public inputs to the circuit
    pub gates: Vec<Gate<F>>,            // constraints list
}

impl<F: Field + ark_ff::FftField> CircuitBuilder<F> {
    pub fn new() -> Self {
        Self {
            next_var_index: 0,
            variables: vec![],
            public_inputs: vec![],
            gates: vec![],
        }
    }

    /// Adds a new variable.
    pub fn new_variable(&mut self, value: Option<F>) -> Variable {
        let var = Variable(self.next_var_index);
        self.next_var_index += 1;
        self.variables.push(value);
        var
    }

    /// Adds a new gate.
    pub fn add_gate(&mut self, gate_type: GateType, a: Variable, b: Variable) -> Variable {
        let eval_out = match gate_type {
            GateType::Add => {
                self.variables[a.0].unwrap() + self.variables[b.0].unwrap()
            }
            GateType::Mul => {
                self.variables[a.0].unwrap() * self.variables[b.0].unwrap()
            }
        };
        let out = self.new_variable(Some(eval_out));
        let row = self.gates.len();

        self.gates.push(Gate {
            gate_type,
            inputs: [a, b],
            output: out,
            selector_row: row,
            constant: F::zero(), // placeholder
        });

        out
    }

    /// Marks a variable as public.
    pub fn mark_public(&mut self, var: Variable) {
        self.public_inputs.push(var);
    }

    /// Evaluates a gate.
    fn evaluate_gate(&self, gate: &Gate<F>) -> F {
        let a = self.variables[gate.inputs[0].0].unwrap();
        let b = self.variables[gate.inputs[1].0].unwrap();
        match gate.gate_type {
            GateType::Add => a + b,
            GateType::Mul => a * b,
        }
    }

    /// Generates the witness table.
    pub fn generate_witness_table(&self) -> WitnessTable<F> {
        let mut a_col = Vec::new();
        let mut b_col = Vec::new();
        let mut c_col = Vec::new();
        let mut q_add = Vec::new();
        let mut q_mul = Vec::new();

        for gate in &self.gates {
            let a = self.variables[gate.inputs[0].0].unwrap();
            let b = self.variables[gate.inputs[1].0].unwrap();
            let c = self.variables[gate.output.0].unwrap();

            a_col.push(a);
            b_col.push(b);
            c_col.push(c);

            match gate.gate_type {
                GateType::Add => {
                    q_add.push(F::one());
                    q_mul.push(F::zero());
                }
                GateType::Mul => {
                    q_add.push(F::zero());
                    q_mul.push(F::one());
                }
            }
        }

        // interpolate to generate selector polynomials
        let q_add_poly = interpolate_selector::<F>(&q_add);
        let q_mul_poly = interpolate_selector::<F>(&q_mul);

        WitnessTable { a_col, b_col, c_col, q_add, q_mul, q_add_poly, q_mul_poly }
    }

    /// Computes the permuation layout.
    pub fn compute_permutation_layout(&self) -> PermutationLayout {
        let mut layout: HashMap<usize, Vec<WirePosition>> = HashMap::new();

        for (row, gate) in self.gates.iter().enumerate() {
            let a = gate.inputs[0].0;
            let b = gate.inputs[1].0;
            let c = gate.output.0;

            layout.entry(a).or_default().push(WirePosition { col: WireColumn::A, row });
            layout.entry(b).or_default().push(WirePosition { col: WireColumn::B, row });
            layout.entry(c).or_default().push(WirePosition { col: WireColumn::C, row });
        }

        PermutationLayout { positions: layout }
    }

}


// ---

pub struct PermutationArgument<F: Field> {
    pub s_id_vals: Vec<F>,
    pub s_sigma_vals: Vec<F>,
    pub z_vals: Vec<F>,
    pub beta: F,
    pub gamma: F,
    pub alpha: F,
}

/// Represents the complete circuit.
pub struct Circuit<F: Field> {
    pub builder: CircuitBuilder<F>,
    pub witness: WitnessTable<F>,
    pub permutation: PermutationLayout,
    pub permutation_argument: Option<PermutationArgument<F>>,
}

impl<F: Field + ark_ff::FftField> Circuit<F> {

    /// Builds circuit using `CircuitBuilder`.
    pub fn from_builder(builder: CircuitBuilder<F>) -> Self {
        let witness = builder.generate_witness_table();
        let permutation = builder.compute_permutation_layout();
        Self { builder, witness, permutation, permutation_argument: None }
    }

     /// Builds gate constraint polynomial.
    pub fn build_gate_constraint(&self) -> DensePolynomial<F> {

        // interpolate to convert to polynomials
        let a_poly = DensePolynomial::from_coefficients_vec(self.witness.a_col.clone());
        let b_poly = DensePolynomial::from_coefficients_vec(self.witness.b_col.clone());
        let c_poly = DensePolynomial::from_coefficients_vec(self.witness.c_col.clone());
        let q_add_poly = self.witness.q_add_poly.clone();
        let q_mul_poly = self.witness.q_mul_poly.clone();

        // q_add * (A + B - C) + q_mul * (A * B - C) = 0
        let add_expr = &q_add_poly * &(a_poly.clone() + b_poly.clone() - c_poly.clone());
        let mul_expr = &q_mul_poly * &((a_poly.clone() * b_poly.clone()) - c_poly.clone());
        add_expr + mul_expr
    }

    /// Builds permutation constraint polynomial.
    pub fn build_permutation_constraint(&self, domain: GeneralEvaluationDomain<F>) -> DensePolynomial<F> {
        let pa = self.permutation_argument.as_ref().expect("Permutation argument not set");
        let a_col = &self.witness.a_col;
        let b_col = &self.witness.b_col;
        let c_col = &self.witness.c_col;
        let n = a_col.len();
        let mut lhs = vec![F::zero(); n];
        let mut rhs = vec![F::zero(); n];

        for i in 0..n {
            let a_term = a_col[i] + pa.beta * pa.s_id_vals[i] + pa.gamma;
            let b_term = b_col[i] + pa.beta * pa.s_id_vals[i] + pa.gamma;
            let c_term = c_col[i] + pa.beta * pa.s_id_vals[i] + pa.gamma;
            lhs[i] = pa.z_vals[i] * a_term * b_term * c_term;

            let a_perm = a_col[i] + pa.beta * pa.s_sigma_vals[i] + pa.gamma;
            let b_perm = b_col[i] + pa.beta * pa.s_sigma_vals[i] + pa.gamma;
            let c_perm = c_col[i] + pa.beta * pa.s_sigma_vals[i] + pa.gamma;
            rhs[i] = pa.z_vals[(i + 1) % n] * a_perm * b_perm * c_perm;
        }

        let evals: Vec<F> = lhs.iter().zip(rhs.iter()).map(|(l, r)| pa.alpha * (*l - *r)).collect();
        DensePolynomial::from_coefficients_vec(domain.ifft(&evals))
    }

    /// Builds public input polynomial.
    pub fn build_public_input_poly(&self, domain: GeneralEvaluationDomain<F>) -> DensePolynomial<F> {
        let pa = self.permutation_argument.as_ref().expect("Permutation argument not set");
        let a_vals = &self.witness.a_col;
        let public_inputs: Vec<F> = self.builder.public_inputs
            .iter()
            .map(|v| self.builder.variables[v.0].unwrap())
            .collect();
        let mut constraint = vec![F::zero(); a_vals.len()];
        for (i, pi) in public_inputs.iter().enumerate() {
            constraint[i] = pa.alpha * (a_vals[i] - *pi);
        }
        DensePolynomial::from_coefficients_vec(domain.ifft(&constraint))
    }

    /// Builds the quotient polynomial.
    pub fn build_quotient_polynomial(&self, domain: GeneralEvaluationDomain<F>) -> DensePolynomial<F> {
        let gate_poly = self.build_gate_constraint();
        let perm_poly = self.build_permutation_constraint(domain);
        let pub_poly = self.build_public_input_poly(domain);

        let t_num = &gate_poly + &perm_poly + &pub_poly;
        // For now, return the numerator polynomial as a placeholder
        // TODO: Implement proper polynomial division with vanishing polynomial
        t_num
    }
}