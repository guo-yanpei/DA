use std::{iter, usize};

use ark_bn254::Fr;
use ark_ff::Field;
use util::mul_group::Radix2Group;

use crate::poly::{MultilinearPoly, UniPolyEvals, UniVarPoly};

pub struct VeriRsProver {
    poly: MultilinearPoly,
    replicas: Vec<UniVarPoly>,
    proofs: Vec<Vec<UniPolyEvals>>,
    final_poly: UniVarPoly,
}

pub struct Symbol {
    poly: MultilinearPoly,
    replica: UniVarPoly,
    proofs: Vec<UniPolyEvals>,
    final_poly: UniVarPoly,
}

impl VeriRsProver {
    fn new(data: Vec<Fr>, layer_num: usize, code_rate: usize, step: usize) -> VeriRsProver {
        let len = data.len();
        let group = Radix2Group::new((len / layer_num).ilog2() as usize + code_rate);
        let codewords = data
            .chunks(len / layer_num)
            .map(|x| group.fft(x.to_vec()))
            .collect::<Vec<_>>();

        let replicas = (0..(len / layer_num << code_rate))
            .map(|i| UniVarPoly::new((0..layer_num).map(|j| codewords[j][i]).collect()))
            .collect::<Vec<_>>();
        let mut claims = replicas
            .iter()
            .map(|x| x.eval(&19260817.into()))
            .collect::<Vec<_>>();
        let inv_2 = <Fr as Field>::inverse(&2.into()).unwrap();
        let nv = claims.len().ilog2() as usize - code_rate;
        let omega = Radix2Group::new(step).element_inv_at(1);
        let mut proofs = vec![];
        for k in 0..(nv / step) {
            let log_replica = claims.len().ilog2() as usize - step;
            let this_proof = (0..(claims.len() >> step))
                .map(|i| {
                    UniPolyEvals::new(
                        (0..(1 << step))
                            .map(|j| claims[i + (j << log_replica)])
                            .collect(),
                        group.element_inv_at(i << (step * k)),
                    )
                })
                .collect::<Vec<_>>();
            claims = this_proof
                .iter()
                .map(|x| x.clone().eval(19260817.into(), omega, inv_2))
                .collect();
            proofs.push(this_proof);
        }

        VeriRsProver {
            poly: MultilinearPoly::new(data),
            replicas,
            proofs,
            final_poly: UniVarPoly::new({
                let log_order = claims.len().ilog2() as usize;
                Radix2Group::new(log_order).ifft(claims)
            }),
        }
    }

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
        let nv = 13;
        let step = 4;
        let log_layer = 3;
        let code_rate = 1;
        let data = (0..(1 << nv))
            .map(|_| <Fr as UniformRand>::rand(&mut rng))
            .collect();
        let prover = VeriRsProver::new(data, 1 << log_layer, code_rate, step);
        let symbol = prover.n_th_replica(103);
        let verifier = VeriRsVerifier::setup(103, step, nv - log_layer + code_rate, code_rate);
        verifier.verify(symbol);
    }
}
