use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use std::vec::Vec;

/// Proof
pub struct Proof {
    // A in 'G1'
    pub a: G1Affine,
    // B in 'G2'
    pub b: G2Affine,
    // C in 'G1'
    pub c: G1Affine,
}

impl Default for Proof {
    fn default() -> self {
        Self {
            a: G1Affine::default(),
            b: G2Affine::default(),
            c: G1Affine::default(),
        }
    }
}

/// ProvingKey
pub struct ProvingKey {
    pub vk: VerifyingKey,
    // beta * G
    pub beta_g1: G1Affine,
    // delta * G
    pub delta_g1: G1Affine,
    // a_i * G
    pub a_query: Vec<G1Affine>,
    // b_i * G
    pub b_g1_query: Vec<G1Affine>,
    // b_i * H
    pub b_g2_query: Vec<G2Affine>,
    // h_i * G
    pub h_query: Vec<G1Affine>,
    // l_i * G
    pub l_query: Vec<G1Affine>,
}
