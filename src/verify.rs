use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use std::vec::Vec;

/// VerifyingKey
pub struct VerifyingKey {
    // alpha * G, where G is the generator of 'G1'
    pub alpha_g1: G1Affine,
    // beta * H, where H is the generator of 'G2'
    pub beta_g2: G2Affine,
    // gamma * H
    pub gamma_g2: G2Affine,
    // delta * H
    pub delta_g2: G2Affine,
    // L(X)=β⋅A(X)+α⋅B(X)+C(X)
    pub gamma_ic_g1: Vec<G1Affine>,
}

impl Default for VerifyingKey {
    fn default() -> self {
        Self {
            alpha_g1: G1Affine::default(),
            beta_g2: G2Affine::default(),
            gamma_g2: G2Affine::default(),
            delta_g2: G2Affine::default(),
            gamma_ic_g1: Vec::new(),
        }
    }
}
