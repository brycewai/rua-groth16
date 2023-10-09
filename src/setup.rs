use crate::common::{polynomial::Polynomial, powers_of_x_on_curve};
use crate::r1cs_to_qap::Qap;
use bls12_381::{G1Affine, G2Affine, Scalar};

/////////////////////////////////////////////////////////////////////
/// CRS
///
// #[derive(Debug)]
pub struct CommonReferenceString<
    const N: usize,
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
> {
    pub alpha_g1: G1Affine,
    pub beta_g1: G1Affine,
    pub delta_g1: G1Affine,
    pub x_power_g1: [G1Affine; N],
    pub public_contribs_g1: [G1Affine; PUBLIC_WITNESS],
    pub private_contribs_g1: [G1Affine; PRIVATE_WITNESS],
    pub x_power_t_x_g1: [G1Affine; N_MINUS_ONE],

    pub beta_g2: G2Affine,
    pub gamma_g2: G2Affine,
    pub delta_g2: G2Affine,
    pub x_power_g2: [G2Affine; N],
}

impl<
        const N: usize,
        const N_MINUS_ONE: usize,
        const PUBLIC_WITNESS: usize,
        const PRIVATE_WITNESS: usize,
    > CommonReferenceString<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>
{
    /// 通过给定参数创建一个新 CRS
    /// 计算公式参考： https://learnblockchain.cn/2019/05/27/groth16/
    pub fn new(
        alpha: Scalar,
        beta: Scalar,
        gamma: Scalar,
        delta: Scalar,
        x: Scalar,
        qap: &Qap<PUBLIC_WITNESS, PRIVATE_WITNESS>,
    ) -> Self {
        let alpha_g1 = (G1Affine::generator() * alpha).into();
        let beta_g1 = (G1Affine::generator() * beta).into();
        let delta_g1: G1Affine = (G1Affine::generator() * delta).into();
        let x_power_g1 = powers_of_x_on_curve(x, G1Affine::generator());
        // public_contribs_g1
        let mut public_contribs_g1 = [G1Affine::generator(); PUBLIC_WITNESS];
        // 创建一个迭代器，该迭代器给出当前迭代次数以及下一个值: (i, val)
        for (i, public_contrib) in public_contribs_g1.iter_mut().enumerate() {
            *public_contrib = ((beta * qap.public_u[i].evaluation(&x)
                + alpha * qap.public_v[i].evaluation(&x)
                + qap.public_w[i].evaluation(&x))
                * gamma.invert().unwrap()
                * G1Affine::generator())
            .into();
        }

        // private_contribs_g1
        let mut private_contribs_g1 = [G1Affine::generator(); PRIVATE_WITNESS];
        for (i, private_contrib) in private_contribs_g1.iter_mut().enumerate() {
            *private_contrib = ((beta * qap.private_u[i].evaluation(&x)
                + alpha * qap.private_v[i].evaluation(&x)
                + qap.private_w[i].evaluation(&x))
                * gamma.invert().unwrap()
                * G1Affine::generator())
            .into();
        }

        // t(x) —— 标量
        let t_x = Polynomial::t_x(N).evaluation(&x);
        // t(x) / delta  —— G1上的点
        let t_x_div_delta_g1 = t_x * delta.invert().unwrap() * G1Affine::generator();
        // t(x) * x^i / delta  —— G1上的点
        let x_power_t_x_g1: [G1Affine; N_MINUS_ONE] =
            powers_of_x_on_curve(x, t_x_div_delta_g1.into());

        let beta_g2 = (G2Affine::generator() * beta).into();
        let gamma_g2 = (G2Affine::generator() * gamma).into();
        let delta_g2 = (G2Affine::generator() * delta).into();
        let x_power_g2 = powers_of_x_on_curve(x, G2Affine::generator());

        CommonReferenceString {
            alpha_g1,
            beta_g1,
            delta_g1,
            x_power_g1,
            public_contribs_g1,
            private_contribs_g1,
            x_power_t_x_g1,

            beta_g2,
            gamma_g2,
            delta_g2,
            x_power_g2,
        }
    }
}

/////////////////////////////////////////////////////////////////////
/// VerifyingKey
/// 本repo中不再实现VerifyingKey，直接使用CRS中的值进行证明计算
#[derive(Clone, Debug, PartialEq)]
pub struct VerifyingKey {
    // alpha * G
    pub alpha_g1: G1Affine,
    // beta * H
    pub beta_g2: G2Affine,
    // gamma * H
    pub gamma_g2: G2Affine,
    // delta * H
    pub delta_g2: G2Affine,
    // gamma^{-1} * (beta * a_i + alpha * b_i + c_i) * G
    pub gamma_abc_g1: Vec<G1Affine>,
}

impl Default for VerifyingKey {
    fn default() -> Self {
        Self {
            alpha_g1: G1Affine::default(),
            beta_g2: G2Affine::default(),
            gamma_g2: G2Affine::default(),
            delta_g2: G2Affine::default(),
            gamma_abc_g1: Vec::new(),
        }
    }
}

/////////////////////////////////////////////////////////////////////
/// ProvingKey
/// 本repo中不再实现provingKey，直接使用CRS中的值进行证明计算
#[derive(Clone, Debug, PartialEq)]
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
