use std::{borrow::Borrow, iter, usize};

use ark_bn254::Fr;
use ark_ff::Field;
use util::mul_group::Radix2Group;

use crate::poly::{MultilinearPoly, UniPolyEvals, UniVarPoly};

pub struct Proofs {
    poly: MultilinearPoly,
    replicas: Vec<UniVarPoly>,
    proofs: Vec<Vec<UniPolyEvals>>,
    final_poly: UniVarPoly,
}

impl Proofs {
    pub fn n_th_replica(&self, n: usize) -> Symbol {
        let proofs = self
            .proofs
            .iter()
            .map(|x| {
                let len = x.len();
                x[n & (len - 1)].clone()
            })
            .collect::<Vec<_>>();
        Symbol {
            poly: self.poly.clone(),
            replica: self.replicas[n].clone(),
            proofs,
            final_poly: self.final_poly.clone(),
        }
    }
}

pub struct VeriRsProver {
    fft_group: Radix2Group,
    evaluate_groups: Vec<Radix2Group>,
    root_inv: Fr,
    code_rate: usize,
    log_row_length: usize,
    log_symbol_number: usize,
    layer_number: usize,
    step: usize,
}

pub struct Symbol {
    poly: MultilinearPoly,
    replica: UniVarPoly,
    proofs: Vec<UniPolyEvals>,
    final_poly: UniVarPoly,
}

impl VeriRsProver {
    pub fn setup(
        log_blob_size: usize,
        log_layer_num: usize,
        code_rate: usize,
        step: usize,
    ) -> Self {
        let fft_group = Radix2Group::new(log_blob_size - log_layer_num + code_rate);
        let log_row_length = log_blob_size - log_layer_num;
        let evaluate_groups = iter::successors(Some(fft_group.clone()), |x| Some(x.exp(1 << step)))
            .take(log_row_length / step)
            .collect::<Vec<_>>();
        let root_inv = Radix2Group::new(step).element_inv_at(1);
        VeriRsProver {
            fft_group,
            root_inv,
            code_rate,
            evaluate_groups,
            log_row_length,
            log_symbol_number: log_blob_size - log_layer_num + code_rate,
            layer_number: 1 << log_layer_num,
            step,
        }
    }

    pub fn prove(&self, data: Vec<Fr>) -> Proofs {
        let codewords = data
            .chunks(1 << self.log_row_length)
            .map(|x| self.fft_group.fft(x.to_vec()))
            .collect::<Vec<_>>();

        let replicas = (0..(1 << self.log_symbol_number))
            .map(|i| UniVarPoly::new((0..self.layer_number).map(|j| codewords[j][i]).collect()))
            .collect::<Vec<_>>();
        let mut claims = replicas
            .iter()
            .map(|x| x.eval(&19260817.into()))
            .collect::<Vec<_>>();
        let inv_2 = <Fr as Field>::inverse(&2.into()).unwrap();
        let mut proofs = vec![];
        for group in &self.evaluate_groups {
            let log_replica = claims.len().ilog2() as usize - self.step;
            let this_proof = (0..(claims.len() >> self.step))
                .map(|i| {
                    UniPolyEvals::new(
                        (0..(1 << self.step))
                            .map(|j| claims[i + (j << log_replica)])
                            .collect(),
                        group.element_inv_at(i),
                    )
                })
                .collect::<Vec<_>>();
            claims = this_proof
                .iter()
                .map(|x| x.clone().eval(19260817.into(), self.root_inv, inv_2))
                .collect();
            proofs.push(this_proof);
        }

        Proofs {
            poly: MultilinearPoly::new(data),
            replicas,
            proofs,
            final_poly: UniVarPoly::new({
                let log_order = claims.len().ilog2() as usize;
                let mut coeff = Radix2Group::new(log_order).ifft(claims);
                coeff.truncate(1 << (log_order - self.code_rate));
                coeff
            }),
        }
    }
}

pub struct VeriRsVerifier {
    locations: Vec<usize>,
    step: usize,
    inv_2: Fr,
    omega_inv: Fr,
    final_point: Fr,
    final_length: usize,
}

impl VeriRsVerifier {
    pub fn setup(
        index: usize,
        step: usize,
        mut log_symbol_number: usize,
        code_rate: usize,
    ) -> VeriRsVerifier {
        let nv = log_symbol_number - code_rate;
        let round = nv / step;
        let mut locations = vec![];
        for _ in 0..round {
            locations.push((index >> (log_symbol_number - step)) & ((1 << step) - 1));
            log_symbol_number -= step;
        }
        let inv_2 = <Fr as Field>::inverse(&2.into()).unwrap();
        let omega_inv = Radix2Group::new(step).element_inv_at(1);
        let final_point =
            Radix2Group::new(log_symbol_number).element_at(index & ((1 << log_symbol_number) - 1));
        VeriRsVerifier {
            locations,
            step,
            inv_2,
            omega_inv,
            final_point,
            final_length: log_symbol_number - code_rate,
        }
    }

    pub fn verify(&self, symbol: Symbol) -> bool {
        let Symbol {
            poly,
            replica,
            proofs,
            final_poly,
        } = symbol;
        let mut x = replica.eval(&19260817.into());
        let mut eval_point = vec![];

        for (poly, &location) in proofs.into_iter().zip(self.locations.iter()) {
            assert_eq!(x, poly.n_th_eval(location));
            x = poly.eval(19260817.into(), self.omega_inv, self.inv_2);
            eval_point.append(
                &mut iter::successors(Some(Fr::from(19260817)), |&x| Some(x * x))
                    .take(self.step)
                    .collect(),
            );
        }

        assert_eq!(x, final_poly.eval(&self.final_point));
        let x = final_poly.eval(&19260817.into());
        eval_point.append(
            &mut iter::successors(Some(Fr::from(19260817)), |&x| Some(x * x))
                .take(self.final_length)
                .collect(),
        );
        eval_point.append(
            &mut iter::successors(Some(Fr::from(19260817)), |&x| Some(x * x))
                .take(replica.len().ilog2() as usize)
                .collect::<Vec<_>>(),
        );
        assert_eq!(x, poly.eval(&eval_point));
        true
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::Fr;
    use ark_ff::UniformRand;
    use rand::thread_rng;

    use super::{VeriRsProver, VeriRsVerifier};

    #[test]
    fn consolidation_test() {
        let mut rng = thread_rng();
        let log_blob_size = 13;
        let step = 4;
        let log_layer_num = 3;
        let code_rate = 1;
        let prover = VeriRsProver::setup(log_blob_size, log_layer_num, code_rate, step);
        let data = (0..(1 << log_blob_size))
            .map(|_| <Fr as UniformRand>::rand(&mut rng))
            .collect::<Vec<_>>();
        let proofs = prover.prove(data);
        let symbol = proofs.n_th_replica(103);
        let verifier = VeriRsVerifier::setup(
            103,
            step,
            log_blob_size - log_layer_num + code_rate,
            code_rate,
        );
        assert!(verifier.verify(symbol));
    }
}
