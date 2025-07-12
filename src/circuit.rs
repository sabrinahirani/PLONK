use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use std::collections::HashMap;

use crate::prover::{WireColumn, WirePosition, PermutationLayout};

use crate::poly_utils::interpolate_selector;

#[derive(Copy, Clone, Debug)]
pub struct Variable(pub usize);

#[derive(Debug)]
pub enum GateType {
    Add,
    Mul
}

#[derive(Debug)]
pub struct Gate<F: Field> {
    pub gate_type: GateType,
    pub inputs: [Variable; 2],
    pub output: Variable,
    pub selector_row: usize, 
    pub constant: F,
}

pub struct CircuitBuilder<F: Field> {
    pub next_var_index: usize,
    pub variables: Vec<Option<F>>, // value assigned to each variable (None until witness gen)
    pub public_inputs: Vec<Variable>,
    pub gates: Vec<Gate<F>>,       // constraint list
}

pub struct WitnessTable<F: Field> {
    pub a_col: Vec<F>,
    pub b_col: Vec<F>,
    pub c_col: Vec<F>,
    pub q_add: Vec<F>,
    pub q_mul: Vec<F>,
    pub q_add_poly: DensePolynomial<F>,
    pub q_mul_poly: DensePolynomial<F>,
}

impl<F: Field> WitnessTable<F> {
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


impl<F: Field + ark_ff::FftField> CircuitBuilder<F> {
    pub fn new() -> Self {
        Self {
            next_var_index: 0,
            variables: vec![],
            public_inputs: vec![],
            gates: vec![],
        }
    }

    pub fn new_variable(&mut self, value: Option<F>) -> Variable {
        let var = Variable(self.next_var_index);
        self.next_var_index += 1;
        self.variables.push(value);
        var
    }

    pub fn add_gate(&mut self, gate_type: GateType, a: Variable, b: Variable) -> Variable {
        let out_val = match gate_type {
            GateType::Add => {
                self.variables[a.0].unwrap() + self.variables[b.0].unwrap()
            }
            GateType::Mul => {
                self.variables[a.0].unwrap() * self.variables[b.0].unwrap()
            }
        };

        let out = self.new_variable(Some(out_val));
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

    pub fn mark_public(&mut self, var: Variable) {
        self.public_inputs.push(var);
    }

    pub fn get_public_input_values(&self) -> Vec<F> {
        self.public_inputs.iter().map(|v| self.variables[v.0].unwrap()).collect()
    }

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

        let q_add_poly = interpolate_selector::<F>(&q_add);
        let q_mul_poly = interpolate_selector::<F>(&q_mul);

        WitnessTable { a_col, b_col, c_col, q_add, q_mul, q_add_poly, q_mul_poly}
    }

    fn evaluate_gate(&self, gate: &Gate<F>) -> F {
        let a = self.variables[gate.inputs[0].0].unwrap();
        let b = self.variables[gate.inputs[1].0].unwrap();
        match gate.gate_type {
            GateType::Add => a + b,
            GateType::Mul => a * b,
        }
    }

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


