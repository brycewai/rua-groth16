use crate::setup::CommonReferenceString;
use crate::{common::polynomial::Polynomial, r1cs_to_qap::Qap};
use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};

/// Proof
#[derive(Debug)]
pub struct Proof {
    // A in 'G1'
    pub a: G1Affine,
    // B in 'G2'
    pub b: G2Affine,
    // C in 'G1'
    pub c: G1Affine,
}

impl Default for Proof {
    fn default() -> Self {
        Self {
            a: G1Affine::default(),
            b: G2Affine::default(),
            c: G1Affine::default(),
        }
    }
}

/// 使用CRS、公共输入和QAP计算生成 Proof(A, B, C)，其中：
/// A = alpha + sum(ai*ui(x)) + r*delta
/// B = beta + sum(ai*vi(x)) + s*delta
/// C = 1/delta * sum(ai*(beta*ui(x) + alpha*vi(x) + wi(x))) + h(x)*t(x) + A*s + B*r - r*s*delta
pub fn prove<
    const N: usize,
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
>(
    crs: &CommonReferenceString<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>,
    qap: Qap<PUBLIC_WITNESS, PRIVATE_WITNESS>,
    public_input: [Scalar; PUBLIC_WITNESS],
    private_input: [Scalar; PRIVATE_WITNESS],
) -> Proof {
    dbg!(public_input);
    let r = Scalar::one();
    let s = Scalar::one();

    let u_x = public_input
        .iter()
        .zip(qap.public_u.iter())
        .map(|(a, u)| *a * u.clone())
        .sum::<Polynomial>()
        + private_input
            .iter()
            .zip(qap.private_u.iter())
            .map(|(a, u)| *a * u.clone())
            .sum();

    let v_x = public_input
        .iter()
        .zip(qap.public_v.iter())
        .map(|(a, v)| *a * v.clone())
        .sum::<Polynomial>()
        + private_input
            .iter()
            .zip(qap.private_v.iter())
            .map(|(a, v)| *a * v.clone())
            .sum();

    let w_x = public_input
        .iter()
        .zip(qap.public_w.iter())
        .map(|(a, w)| *a * w.clone())
        .sum::<Polynomial>()
        + private_input
            .iter()
            .zip(qap.private_w.iter())
            .map(|(a, w)| *a * w.clone())
            .sum();

    let p_x = u_x.clone() * v_x.clone() - w_x;
    let h_x = &p_x / &Polynomial::t_x(N); // 多项式除法
    assert!(h_x.remainder == Polynomial::zero()); // 余数不为0

    let a = (crs.alpha_g1
        + G1Projective::from(u_x.evaluation_secret_power(&crs.x_power_g1))
        + r * crs.delta_g1)
        .into();
    let b = (crs.beta_g2
        + G2Projective::from(v_x.evaluation_secret_power(&crs.x_power_g2))
        + s * crs.delta_g2)
        .into();
    let b_g1: G1Affine = (crs.beta_g1
        + G1Projective::from(v_x.evaluation_secret_power(&crs.x_power_g1))
        + s * crs.delta_g1)
        .into();
    let c = G1Affine::from(
        private_input
            .iter()
            .zip(crs.private_contribs_g1.iter())
            .map(|(a, contrib)| *a * *contrib)
            .sum::<G1Projective>()
            + h_x.quotient.evaluation_secret_power(&crs.x_power_t_x_g1)
            + a * s
            + r * b_g1
            - r * s * crs.delta_g1,
    );
    Proof { a, b, c }
}
