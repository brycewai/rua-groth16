use crate::common::powers_of_x_on_curve;
use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};

#[drive(Debug)]
struct CommonReferenceString<
    const N: usize,
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
> {
    alpha_g1: G1Affine,
    beta_g1: G1Affine,
    delta_g1: G1Affine,
    x_power_g1: [G1Affine; N],
    public_contribs_g1: [G1Affine; PUBLIC_WITNESS],
    private_contribs_g1: [G1Affine; PRIVATE_WITNESS],
    x_power_t_x_g1: [G1Affine; N_MINUS_ONE],

    beta_g2: G2Affine,
    gamma_g2: G2Affine,
    delta_g2: G2Affine,
    x_power_g2: [G2Affine; N],
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
    pub fn new(alpha: Scalar, beta: Scalar, gamma: Scalar, delta: Scalar, x: Scalar) -> Self {
        let alpha_g1 = (G1Affine::generator() * alpha).into();
        let beta_g1 = (G1Affine::generator() * beta).into();
        let delta_g1 = (G1Affine::generator() * delta).into();
        let x_power_g1 = powers_of_x_on_curve(x, Curve::generator());
        // public_contribs_g1
        let mut public_contribs_g1 = [G1Affine::generator(); PUBLIC_WITNESS];
        // 创建一个迭代器，该迭代器给出当前迭代次数以及下一个值: (i, val)
        for (i, public_contrib) in public_contribs_g1.iter_mut().enumerate() {
            // *public_contrib = beta * u[i]
        }

        // private_contribs_g1
        let mut private_contribs_g1 = [G1Affine::generator(); PRIVATE_WITNESS];

        // t(x) —— 标量
        let mut t_x = Polynomial::t_x(N).evaluation(&x);
        // t(x) / delta  —— G1上的点
        let mut t_x_div_delta = t_x * delta_g1.invert().unwarp() * G1Affine::generator();
        // t(x) * x^i / delta  —— G1上的点
        let mut t_x_x_power = powers_of_x_on_curve(x, t_x_div_delta.into());

        let beta_g2 = (G2Affine::generator() * beta).into();
        let gamma_g2 = (G2Affine::generator() * gamma).into();
        let delta_g2 = (G2Affine::generator() * delta).into();
        let x_power_g2 = powers_of_x_on_curve(x, Curve::generator());
    }
}
