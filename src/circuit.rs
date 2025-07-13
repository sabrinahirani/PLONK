use ark_ff::{Field, FftField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly::DenseUVPolynomial;
use ark_poly::polynomial::Polynomial;

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

        // Initialize sigma as identity mapping
        let mut sigma: Vec<usize> = (0..3 * num_rows).collect();

        // Overwrite only the positions involved in cycles
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

        // sigma is now guaranteed to have length 3 * num_rows and include all positions
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
    pub fn generate_witness_table(&self, domain_size: usize) -> WitnessTable<F> {
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

        // Pad all columns to domain_size with zeros
        while a_col.len() < domain_size {
            a_col.push(F::zero());
            b_col.push(F::zero());
            c_col.push(F::zero());
            q_add.push(F::zero());
            q_mul.push(F::zero());
        }

        WitnessTable { a_col, b_col, c_col, q_add, q_mul }
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
pub struct Circuit<F: Field + FftField> {
    pub builder: CircuitBuilder<F>,
    pub witness: WitnessTable<F>,
    pub permutation: PermutationLayout,
    pub permutation_argument: Option<PermutationArgument<F>>,
    pub domain: GeneralEvaluationDomain<F>,
}

impl<F: Field + ark_ff::FftField> Circuit<F> {

    /// Builds circuit using `CircuitBuilder`.
    pub fn from_builder(builder: CircuitBuilder<F>, domain: GeneralEvaluationDomain<F>) -> Self {
        let witness = builder.generate_witness_table(domain.size());
        let permutation = builder.compute_permutation_layout();
        Self { builder, witness, permutation, permutation_argument: None, domain }
    }

     /// Builds gate constraint polynomial.
    pub fn build_gate_constraint(&self) -> DensePolynomial<F> {
        println!("a_col: {:?}", self.witness.a_col);
        println!("b_col: {:?}", self.witness.b_col);
        println!("c_col: {:?}", self.witness.c_col);
        println!("q_add: {:?}", self.witness.q_add);
        println!("q_mul: {:?}", self.witness.q_mul);
        println!("a_col len: {}", self.witness.a_col.len());
        println!("b_col len: {}", self.witness.b_col.len());
        println!("c_col len: {}", self.witness.c_col.len());
        println!("q_add len: {}", self.witness.q_add.len());
        println!("q_mul len: {}", self.witness.q_mul.len());
        println!("a_col first 5: {:?}", &self.witness.a_col[..5.min(self.witness.a_col.len())]);
        println!("b_col first 5: {:?}", &self.witness.b_col[..5.min(self.witness.b_col.len())]);
        println!("c_col first 5: {:?}", &self.witness.c_col[..5.min(self.witness.c_col.len())]);
        println!("q_add first 5: {:?}", &self.witness.q_add[..5.min(self.witness.q_add.len())]);
        println!("q_mul first 5: {:?}", &self.witness.q_mul[..5.min(self.witness.q_mul.len())]);

        let n = self.witness.a_col.len();
        let mut gate_vals = vec![F::zero(); n];
        for i in 0..n {
            let a = self.witness.a_col[i];
            let b = self.witness.b_col[i];
            let c = self.witness.c_col[i];
            let q_add = self.witness.q_add[i];
            let q_mul = self.witness.q_mul[i];
            gate_vals[i] = q_add * (a + b - c) + q_mul * (a * b - c);
        }
        let gate_poly = DensePolynomial::from_coefficients_vec(self.domain.ifft(&gate_vals));
        gate_poly
    }

    /// Builds permutation constraint polynomial.
    pub fn build_permutation_constraint(
        &self,
        a_col_flat: &[F],
        b_col_flat: &[F],
        c_col_flat: &[F],
        sigma: &[usize],
    ) -> DensePolynomial<F> {
        let pa = self.permutation_argument.as_ref().expect("Permutation argument not set");
        let n = self.domain.size();
    
        let mut constraint_vals = vec![F::zero(); n];
    
        for i in 0..n {
            let a = a_col_flat[i];
            let b = b_col_flat[i];
            let c = c_col_flat[i];

            // Use explicit wire indices for the identity permutation
            let a_term = a + pa.beta * F::from(3 * i as u64) + pa.gamma;
            let b_term = b + pa.beta * F::from(3 * i as u64 + 1) + pa.gamma;
            let c_term = c + pa.beta * F::from(3 * i as u64 + 2) + pa.gamma;

            // For sigma, use the s_id_vals at the sigma indices
            let a_term_sigma = a + pa.beta * pa.s_id_vals[sigma[3 * i]] + pa.gamma;
            let b_term_sigma = b + pa.beta * pa.s_id_vals[sigma[3 * i + 1]] + pa.gamma;
            let c_term_sigma = c + pa.beta * pa.s_id_vals[sigma[3 * i + 2]] + pa.gamma;

            if i < n - 1 {
                // Standard permutation constraint
                let lhs = pa.z_vals[i] * a_term * b_term * c_term;
                let rhs = pa.z_vals[i + 1] * a_term_sigma * b_term_sigma * c_term_sigma;
                constraint_vals[i] = lhs - rhs;
            } else {
                // Boundary condition: z[n] = 1, so we check that the last product equals 1
                let product = a_term * b_term * c_term * (a_term_sigma * b_term_sigma * c_term_sigma).inverse().unwrap();
                constraint_vals[i] = pa.z_vals[i] * product - F::one();
            }

            if i < 4 {
                println!("Row {} debug:", i);
                println!("  a: {}  b: {}  c: {}", a, b, c);
                println!("  a_term: {} (index {})", a_term, 3 * i);
                println!("  b_term: {} (index {})", b_term, 3 * i + 1);
                println!("  c_term: {} (index {})", c_term, 3 * i + 2);
                println!("  sigma indices: [{}, {}, {}]", sigma[3 * i], sigma[3 * i + 1], sigma[3 * i + 2]);
                println!("  a_term_sigma: {} (index {})", a_term_sigma, sigma[3 * i]);
                println!("  b_term_sigma: {} (index {})", b_term_sigma, sigma[3 * i + 1]);
                println!("  c_term_sigma: {} (index {})", c_term_sigma, sigma[3 * i + 2]);
                if i < n - 1 {
                    let lhs = pa.z_vals[i] * a_term * b_term * c_term;
                    let rhs = pa.z_vals[i + 1] * a_term_sigma * b_term_sigma * c_term_sigma;
                    println!("  lhs: {}", lhs);
                    println!("  rhs: {}", rhs);
                    println!("  lhs - rhs = {}", lhs - rhs);
                } else {
                    let product = a_term * b_term * c_term * (a_term_sigma * b_term_sigma * c_term_sigma).inverse().unwrap();
                    println!("  boundary constraint: {} * product - 1 = {}", pa.z_vals[i], pa.z_vals[i] * product - F::one());
                }
            }
        }
    
        DensePolynomial::from_coefficients_vec(self.domain.ifft(&constraint_vals))
    }
    
    

    /// Builds public input constraint polynomial: alpha * (wire_value - public_input)
    pub fn build_public_input_poly(&self) -> DensePolynomial<F> {
        let pa = self.permutation_argument.as_ref().expect("Permutation argument not set");
        let a_vals = &self.witness.a_col;
        let public_inputs: Vec<F> = self.builder.public_inputs
            .iter()
            .map(|v| self.builder.variables[v.0].unwrap())
            .collect();
        println!("public_inputs: {:?}", public_inputs);
        println!("public_input indices: {:?}", self.builder.public_inputs.iter().map(|v| v.0).collect::<Vec<_>>());
        let mut constraint = vec![F::zero(); a_vals.len()];
        for (pi_idx, var) in self.builder.public_inputs.iter().enumerate() {
            let pi_value = self.builder.variables[var.0].unwrap();
            // Find the first row where this variable is used as an input or output to a gate
            let mut found_row = None;
            for (row, gate) in self.builder.gates.iter().enumerate() {
                if gate.inputs[0].0 == var.0 || gate.inputs[1].0 == var.0 || gate.output.0 == var.0 {
                    found_row = Some(row);
                    break;
                }
            }
            if let Some(row) = found_row {
                println!("Public input variable {} (value {}) mapped to row {}", var.0, pi_value, row);
                constraint[row] = pa.alpha * (a_vals[row] - pi_value);
            } else {
                println!("Warning: public input variable {} not found in any gate input or output", var.0);
            }
        }
        println!("public input constraint vector: {:?}", constraint);
        DensePolynomial::from_coefficients_vec(self.domain.ifft(&constraint))
    }

    /// Builds the quotient polynomial.
    pub fn build_quotient_polynomial(&self, sigma: &[usize]) -> DensePolynomial<F> {
        // Compute gate constraint values at each domain point (still use self.witness.a_col, etc. for gate domain)
        let n = self.witness.a_col.len();
        let mut gate_vals = vec![F::zero(); n];
        for i in 0..n {
            let a = self.witness.a_col[i];
            let b = self.witness.b_col[i];
            let c = self.witness.c_col[i];
            let q_add = self.witness.q_add[i];
            let q_mul = self.witness.q_mul[i];
            gate_vals[i] = q_add * (a + b - c) + q_mul * (a * b - c);
        }
        let gate_poly = DensePolynomial::from_coefficients_vec(self.domain.ifft(&gate_vals));

        // Compute permutation and public input constraint values at each domain point
        let perm_poly = self.build_permutation_constraint(&self.witness.a_col, &self.witness.b_col, &self.witness.c_col, &sigma); // DensePolynomial<F>
        let pub_poly = self.build_public_input_poly();       // DensePolynomial<F>
        let perm_vals = self.domain.fft(&perm_poly.coeffs);
        let pub_vals = self.domain.fft(&pub_poly.coeffs);

        // Print all values and their sum at each domain point
        println!("i | gate | perm | pub | sum");
        for i in 0..n {
            let sum = gate_vals[i] + perm_vals[i] + pub_vals[i];
            println!("{:2} | {} | {} | {} | {}", i, gate_vals[i], perm_vals[i], pub_vals[i], sum);
        }

        let t_num = gate_poly.clone() + &perm_poly + &pub_poly;
        println!("gate_poly degree: {}", gate_poly.degree());
        println!("perm_poly degree: {}", perm_poly.degree());
        println!("pub_poly degree: {}", pub_poly.degree());
        println!("gate_poly first 5 coeffs: {:?}", &gate_poly.coeffs[..5.min(gate_poly.coeffs.len())]);
        println!("perm_poly first 5 coeffs: {:?}", &perm_poly.coeffs[..5.min(perm_poly.coeffs.len())]);
        println!("pub_poly first 5 coeffs: {:?}", &pub_poly.coeffs[..5.min(pub_poly.coeffs.len())]);
        println!("t_num degree: {}", t_num.degree());
        println!("t_num first 5 coeffs: {:?}", &t_num.coeffs[..5.min(t_num.coeffs.len())]);

        let z_h: DensePolynomial<F> = self.domain.vanishing_polynomial().into();
        let (t_quotient, remainder) = t_num.divide_by_vanishing_poly(self.domain);
        if !remainder.is_zero() {
            println!("Nonzero remainder: {:?}", remainder);
            for (i, coeff) in remainder.coeffs.iter().enumerate() {
                println!("  remainder coeff[{}] = {}", i, coeff);
            }
        }
        assert!(remainder.is_zero(), "t(X) not divisible by Z_H(X)");
        t_quotient
    }

    pub fn build_grand_product(
        witness_flat: &[F],
        sigma: &[usize],
        domain: GeneralEvaluationDomain<F>,
        beta: F,
        gamma: F,
        s_id_vals: &[F],
    ) -> DensePolynomial<F> {
        let n = domain.size();
        let mut z = vec![F::one(); n + 1];
    
        for i in 0..n {
            // Flattened witness layout: [A_0, B_0, C_0, A_1, B_1, C_1, ...]
            let a = witness_flat[3 * i];
            let b = witness_flat[3 * i + 1];
            let c = witness_flat[3 * i + 2];
    
            let a_sigma = witness_flat[sigma[3 * i]];
            let b_sigma = witness_flat[sigma[3 * i + 1]];
            let c_sigma = witness_flat[sigma[3 * i + 2]];
    
            let a_term = a + beta * s_id_vals[3 * i] + gamma;
            let b_term = b + beta * s_id_vals[3 * i + 1] + gamma;
            let c_term = c + beta * s_id_vals[3 * i + 2] + gamma;
    
            let a_term_sigma = a_sigma + beta * s_id_vals[sigma[3 * i]] + gamma;
            let b_term_sigma = b_sigma + beta * s_id_vals[sigma[3 * i + 1]] + gamma;
            let c_term_sigma = c_sigma + beta * s_id_vals[sigma[3 * i + 2]] + gamma;
    
            let numerator = a_term * b_term * c_term;
            let denominator = a_term_sigma * b_term_sigma * c_term_sigma;
    
            z[i + 1] = z[i] * numerator * denominator.inverse().unwrap();

            if i < 4 {
                println!("Row {} debug:", i);
                println!("  a: {}  b: {}  c: {}", a, b, c);
                println!("  s_id_vals: [{}, {}, {}]", s_id_vals[3 * i], s_id_vals[3 * i + 1], s_id_vals[3 * i + 2]);
                println!("  a_term: {} (s_id: {})", a_term, s_id_vals[3 * i]);
                println!("  b_term: {} (s_id: {})", b_term, s_id_vals[3 * i + 1]);
                println!("  c_term: {} (s_id: {})", c_term, s_id_vals[3 * i + 2]);
                println!("  sigma indices: [{}, {}, {}]", sigma[3 * i], sigma[3 * i + 1], sigma[3 * i + 2]);
                println!("  a_term_sigma: {} (sigma: {})", a_term_sigma, sigma[3 * i]);
                println!("  b_term_sigma: {} (sigma: {})", b_term_sigma, sigma[3 * i + 1]);
                println!("  c_term_sigma: {} (sigma: {})", c_term_sigma, sigma[3 * i + 2]);
                println!("  lhs (z[i+1]): {}", z[i + 1]);
                let rhs = z[i] * numerator * denominator.inverse().unwrap();
                println!("  rhs (should match lhs): {}", rhs);
                println!("  lhs - rhs = {}", z[i + 1] - rhs);
            }
        }
    
        DensePolynomial::from_coefficients_vec(domain.ifft(&z[..n].to_vec()))
    }
    
}