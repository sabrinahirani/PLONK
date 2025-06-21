use ark_bn254::Fr;

#[derive(Clone, Debug)]
pub enum GateType {
    Input,
    Add,
    Mul,
    ConstMul(Fr),
    Output,
}

#[derive(Clone, Debug)]
pub struct Gate {
    pub id: usize,
    pub gate_type: GateType,
    pub left: Option<usize>,
    pub right: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct Circuit {
    pub gates: Vec<Gate>,

    pub witness: Vec<Fr>,

    pub a_vec: Vec<Fr>,
    pub b_vec: Vec<Fr>,
    pub c_vec: Vec<Fr>,
}

impl Circuit {
    pub fn new() -> Self {
        Self {
            gates: Vec::new(),
            witness: Vec::new(),
            a_vec: Vec::new(),
            b_vec: Vec::new(),
            c_vec: Vec::new(),
        }
    }

    pub fn add_gate(
        &mut self,
        gate_type: GateType,
        left: Option<usize>,
        right: Option<usize>,
        owner: Option<usize>,
    ) -> usize {
        let id = self.gates.len();
        self.gates.push(Gate {
            id,
            gate_type,
            left,
            right,
            owner,
        });
        id
    }

    pub fn output_wires(&self) -> Vec<usize> {
        self.gates
            .iter()
            .filter(|g| matches!(g.gate_type, GateType::Output))
            .map(|g| g.id)
            .collect()
    }

    pub fn assign_witness(&mut self, witness: Vec<Fr>) {
        assert_eq!(
            witness.len(),
            self.gates.len(),
            "Witness length must match gate count"
        );
        self.witness = witness;
    }

    pub fn synthesize_r1cs(&mut self) {
        let n = self.gates.len();
        self.a_vec = vec![Fr::zero(); n];
        self.b_vec = vec![Fr::zero(); n];
        self.c_vec = vec![Fr::zero(); n];

        for gate in &self.gates {
            match &gate.gate_type {
                GateType::Input | GateType::Output => {
                    // No R1CS constraint for input/output markers
                    continue;
                }

                GateType::Add => {
                    // Enforce: left + right = output
                    self.a_vec[gate.id] = Fr::one();  // a · left
                    self.b_vec[gate.id] = Fr::one();  // b · right
                    self.c_vec[gate.id] = -Fr::one(); // c · output
                }

                GateType::Mul => {
                    // Enforce: left * right = output
                    let a_val = self.get_witness(gate.left);
                    let b_val = self.get_witness(gate.right);
                    let out_val = self.get_witness(Some(gate.id));

                    self.a_vec[gate.id] = a_val;
                    self.b_vec[gate.id] = b_val;
                    self.c_vec[gate.id] = -out_val;
                }

                GateType::ConstMul(constant) => {
                    // Enforce: left * constant = output
                    let left_val = self.get_witness(gate.left);
                    let out_val = self.get_witness(Some(gate.id));

                    self.a_vec[gate.id] = *constant;
                    self.b_vec[gate.id] = left_val;
                    self.c_vec[gate.id] = -out_val;
                }
            }
        }
    }

    /// Helper to safely get a witness value
    fn get_witness(&self, idx: Option<usize>) -> Fr {
        idx.map(|i| self.witness[i]).unwrap_or(Fr::zero())
    }

    /// Converts a_vec, b_vec, c_vec into DensePolynomials using Lagrange interpolation
    pub fn r1cs_to_polynomials(&self, domain: Radix2EvaluationDomain<Fr>) -> (DensePolynomial<Fr>, DensePolynomial<Fr>, DensePolynomial<Fr>) {
        assert_eq!(self.a_vec.len(), domain.size());
        assert_eq!(self.b_vec.len(), domain.size());
        assert_eq!(self.c_vec.len(), domain.size());

        let a_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&self.a_vec));
        let b_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&self.b_vec));
        let c_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&self.c_vec));

        (a_poly, b_poly, c_poly)
    }
}
