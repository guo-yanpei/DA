use std::collections::HashMap;

use ark_bn254::Fr;
use ark_ff::{Field, Zero};
use util::{
    merkle_tree::{MerkleTreeProver, MerkleTreeVerifier, MERKLE_ROOT_SIZE},
    mul_group::Radix2Group,
};

pub struct QueryResult {
    paths: Vec<u8>,
    values: HashMap<usize, Fr>,
}

impl QueryResult {
    pub fn verify_merkle_tree(
        &self,
        leaf_indices: &Vec<usize>,
        leaf_size: usize,
        merkle_verifier: &MerkleTreeVerifier,
    ) -> bool {
        let len = merkle_verifier.leave_number;

        let leaves: Vec<Vec<u8>> = leaf_indices
            .iter()
            .map(|x| {
                MerkleTreeProver::serialize_fields(
                    &(0..leaf_size)
                        .map(|j| self.values.get(&(x.clone() + j * len)).unwrap().clone())
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let res = merkle_verifier.verify(self.paths.clone(), leaf_indices, &leaves);
        assert!(res);
        res
    }
}

pub struct InterpolateValue {
    pub value: Vec<Fr>,
    leaf_size: usize,
    merkle_tree: MerkleTreeProver,
}

impl InterpolateValue {
    pub fn new(value: Vec<Fr>, leaf_size: usize) -> Self {
        let len = value.len() / leaf_size;
        let mt = MerkleTreeProver::new(
            &(0..len)
                .map(|i| {
                    MerkleTreeProver::serialize_fields(
                        &(0..leaf_size)
                            .map(|j| value[len * j + i])
                            .collect::<Vec<_>>(),
                    )
                })
                .collect(),
        );
        InterpolateValue {
            value,
            leaf_size,
            merkle_tree: mt,
        }
    }

    pub fn leave_num(&self) -> usize {
        self.merkle_tree.leave_num()
    }

    pub fn commit(&self) -> [u8; MERKLE_ROOT_SIZE] {
        self.merkle_tree.commit()
    }

    pub fn query(&self, leaf_indices: &Vec<usize>) -> QueryResult {
        let len = self.merkle_tree.leave_num();
        assert_eq!(len * self.leaf_size, self.value.len());
        let proof_values = (0..self.leaf_size)
            .flat_map(|i| {
                leaf_indices
                    .iter()
                    .map(|j| (j.clone() + i * len, self.value[j.clone() + i * len]))
                    .collect::<Vec<_>>()
            })
            .collect();
        let proof_bytes = self.merkle_tree.open(&leaf_indices);
        QueryResult {
            paths: proof_bytes,
            values: proof_values,
        }
    }
}

pub struct Proof {
    interpolations: Vec<[u8; MERKLE_ROOT_SIZE]>,
    final_value: Fr,
    query_results: Vec<QueryResult>,
}

pub struct Prover {
    interpolation: InterpolateValue,
    poly_num: usize,
    log_degree: usize,
}

impl Prover {
    fn evaluate_next_domain(
        last_interpolation: &Vec<Fr>,
        group: &Radix2Group,
        inv_2: Fr,
        challenge: Fr,
    ) -> Vec<Fr> {
        let mut res = vec![];
        let len = group.size();
        for i in 0..(len / 2) {
            let x = last_interpolation[i];
            let nx = last_interpolation[i + len / 2];
            let sum = x + nx;
            let new_v = sum + challenge * ((x - nx) * group.element_inv_at(i) - sum);
            res.push(new_v * inv_2);
        }
        res
    }

    pub fn new(polies: &[Vec<Fr>], group: &Radix2Group) -> Self {
        let log_degree = polies[0].len().ilog2() as usize;
        let value = polies
            .iter()
            .flat_map(|x| group.fft(x.clone()))
            .collect::<Vec<_>>();
        Prover {
            interpolation: InterpolateValue::new(value, polies.len() * 2),
            poly_num: polies.len(),
            log_degree,
        }
    }

    fn commit(&self) -> [u8; MERKLE_ROOT_SIZE] {
        self.interpolation.commit()
    }

    fn prove(
        &self,
        groups: &Vec<Radix2Group>,
        challenges: &(Fr, Vec<Fr>),
        mut leaf_indices: Vec<usize>,
    ) -> Proof {
        let poly_interpolations = {
            let len = groups[0].size();
            let mut v = (0..len).map(|_| <Fr as Zero>::zero()).collect::<Vec<_>>();
            for i in 0..len {
                let mut j = i;
                for _ in 0..self.poly_num {
                    v[i] *= challenges.0;
                    v[i] += self.interpolation.value[j];
                    j += len;
                }
            }
            v
        };
        let mut interpolations: Vec<InterpolateValue> = vec![];
        let mut final_value = None;
        let inv_2 = <Fr as Field>::inverse(&2.into()).unwrap();
        for i in 0..self.log_degree {
            let next_evaluation = Self::evaluate_next_domain(
                if i == 0 {
                    &poly_interpolations
                } else {
                    &interpolations[i - 1].value
                },
                &groups[i],
                inv_2,
                challenges.1[i],
            );
            if i < self.log_degree - 1 {
                let new_interpolation = InterpolateValue::new(next_evaluation, 2);
                interpolations.push(new_interpolation);
            } else {
                final_value = Some(next_evaluation[0]);
            }
        }
        let mut query_results = vec![];
        for i in 0..self.log_degree {
            let len = groups[i].size();
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();
            if i == 0 {
                query_results.push(self.interpolation.query(&leaf_indices));
            } else {
                query_results.push(interpolations[i - 1].query(&leaf_indices));
            }
        }
        Proof {
            interpolations: interpolations.iter().map(|x| x.commit()).collect(),
            final_value: final_value.unwrap(),
            query_results,
        }
    }
}

pub struct Verifier {
    mt_verifier: MerkleTreeVerifier,
    poly_num: usize,
}

impl Verifier {
    pub fn new(merkle_root: [u8; MERKLE_ROOT_SIZE], poly_num: usize, leave_number: usize) -> Self {
        Verifier {
            mt_verifier: MerkleTreeVerifier::new(leave_number, &merkle_root),
            poly_num,
        }
    }

    pub fn verify(
        &self,
        groups: &Vec<Radix2Group>,
        challenges: &(Fr, Vec<Fr>),
        mut leaf_indices: Vec<usize>,
        proof: Proof,
    ) {
        let mt_verifiers = {
            let mut v = vec![];
            let mut leave_num = self.mt_verifier.leave_number;
            for hash in proof.interpolations.iter() {
                leave_num /= 2;
                v.push(MerkleTreeVerifier::new(leave_num, hash));
            }
            v
        };

        let log_degree = challenges.1.len();
        for i in 0..log_degree {
            let len = groups[i].size();
            leaf_indices = leaf_indices.iter_mut().map(|v| *v % (len >> 1)).collect();
            leaf_indices.sort();
            leaf_indices.dedup();

            proof.query_results[i].verify_merkle_tree(
                &leaf_indices,
                if i == 0 { self.poly_num * 2 } else { 2 },
                if i == 0 {
                    &self.mt_verifier
                } else {
                    &mt_verifiers[i - 1]
                },
            );

            for j in leaf_indices.iter() {
                let new_v = if i == 0 {
                    let mut res = Fr::from(0);
                    let mut k = j.clone();
                    for _ in 0..self.poly_num {
                        let x = proof.query_results[0].values.get(&k).unwrap().clone();
                        let nx = proof.query_results[0]
                            .values
                            .get(&(k + len / 2))
                            .unwrap()
                            .clone();
                        let sum = x + nx;
                        res *= challenges.0;
                        res +=
                            sum + challenges.1[0] * ((x - nx) * groups[0].element_inv_at(*j) - sum);
                        k += len;
                    }
                    res
                } else {
                    let x = proof.query_results[i].values.get(&j).unwrap().clone();
                    let nx = proof.query_results[i]
                        .values
                        .get(&(j + len / 2))
                        .unwrap()
                        .clone();
                    let sum = x + nx;
                    sum + challenges.1[i] * ((x - nx) * groups[i].element_inv_at(*j) - sum)
                };
                if i < log_degree - 1 {
                    assert_eq!(new_v, proof.query_results[i + 1].values[j].double());
                } else {
                    assert_eq!(new_v, proof.final_value.double());
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use rand::{thread_rng, RngCore};

    use super::*;

    #[test]
    fn it_works() {
        let mut rng = thread_rng();
        let poly_num = 16;
        let log_degree = 12;
        let polies = (0..poly_num)
            .map(|x| {
                (0..(1 << log_degree))
                    .map(|_| <Fr as UniformRand>::rand(&mut rng))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let coderate = 1;
        let groups = (0..log_degree)
            .rev()
            .map(|x| Radix2Group::new(x + 1 + coderate))
            .collect::<Vec<_>>();
        let prover = Prover::new(&polies, &groups[0]);
        let challenges = {
            (
                <Fr as UniformRand>::rand(&mut rng),
                (0..log_degree)
                    .map(|_| <Fr as UniformRand>::rand(&mut rng))
                    .collect::<Vec<_>>(),
            )
        };
        let leaf_indices = (0..30).map(|_| rng.next_u32() as usize).collect::<Vec<_>>();
        let proof = prover.prove(&groups, &challenges, leaf_indices.clone());
        let commit = prover.commit();
        let verifier = Verifier::new(commit, poly_num, 1 << (log_degree + coderate - 1));
        verifier.verify(&groups, &challenges, leaf_indices, proof);
    }
}
