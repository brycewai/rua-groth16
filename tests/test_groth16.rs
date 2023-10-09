use bls12_381::Scalar;
use groth16::prove::prove;
use groth16::r1cs_to_qap::{ConstraintMatrix, Qap, R1cs};
use groth16::setup::CommonReferenceString;
use groth16::verify::verify;

#[test]
fn groth16() {
    pub const R1CS: R1cs<1, 1, 2> = R1cs::new(
        ConstraintMatrix {
            private_entries: [[Scalar::one(), Scalar::one()]],
            public_entries: [[Scalar::zero(), Scalar::zero()]],
        },
        ConstraintMatrix {
            private_entries: [[Scalar::one(), Scalar::zero()]],
            public_entries: [[Scalar::zero(), Scalar::one()]],
        },
        ConstraintMatrix {
            private_entries: [[Scalar::zero(), Scalar::one()]],
            public_entries: [[Scalar::one(), Scalar::zero()]],
        },
    );

    let qap: Qap<1, 1> = R1CS.into();
    let crs = CommonReferenceString::<2, 1, 1, 1>::new(
        123.into(),
        456.into(),
        789.into(),
        777.into(),
        4756.into(),
        &qap,
    );
    let public_input: [Scalar; 1] = [Scalar::one()];
    let private_input: [Scalar; 1] = [Scalar::one()];
    let proof = prove(&crs, qap, public_input, private_input);
    dbg!(&proof);
    let res = verify(crs, public_input, proof);

    assert!(res);
}
