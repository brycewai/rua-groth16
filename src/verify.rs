use crate::{prove::Proof, setup::CommonReferenceString};
use bls12_381::{G1Affine, G1Projective, Scalar};
use pairing::PairingCurveAffine;

/////////////////////////////////////////////////////////////////////
/// verify
/// 通过crs、public input、proof来计算验证proof是否正确
pub fn verify<
    const N: usize,
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
>(
    crs: &CommonReferenceString<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>,
    public_input: &[Scalar; PUBLIC_WITNESS],
    proof: &Proof,
) -> bool {
    dbg!(public_input);
    let left = proof.a.pairing_with(&proof.b);
    dbg!(left);
    let right = crs.alpha_g1.pairing_with(&crs.beta_g2)
        + G1Affine::from(
            public_input
                .iter()
                .zip(crs.public_contribs_g1.iter())
                .map(|(a, contrib)| a * contrib)
                .sum::<G1Projective>(),
        )
        .pairing_with(&crs.gamma_g2)
        + proof.c.pairing_with(&crs.delta_g2);
    dbg!(right);
    left == right
}
