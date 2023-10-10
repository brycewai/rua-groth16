/// Reference: <https://github.com/ThomasdenH/groth16>
use bls12_381::Scalar;
use groth16::prove::prove;
use groth16::r1cs_to_qap::{ConstraintMatrix, Qap, R1cs};
use groth16::setup::CommonReferenceString;
use groth16::verify::verify;
fn main() {
    // 1. compile circuit (...)
    // 2. setup (generate CRS)
    // 3. compute witness (private/public input)
    pub const R1CS: R1cs<2, 1, 2> = R1cs::new(
        ConstraintMatrix {
            private_entries: [[Scalar::one(), Scalar::one()]],
            public_entries: [
                [Scalar::zero(), Scalar::zero()],
                [Scalar::one(), Scalar::one()],
            ],
        },
        ConstraintMatrix {
            private_entries: [[Scalar::one(), Scalar::zero()]],
            public_entries: [
                [Scalar::zero(), Scalar::zero()],
                [Scalar::one(), Scalar::one()],
            ],
        },
        ConstraintMatrix {
            private_entries: [[Scalar::zero(), Scalar::one()]],
            public_entries: [
                [Scalar::zero(), Scalar::zero()],
                [Scalar::one(), Scalar::one()],
            ],
        },
    );

    let qap: Qap<2, 1> = R1CS.into();
    let crs = CommonReferenceString::<3, 2, 2, 1>::new(
        12345.into(),
        3456.into(),
        6789.into(),
        77.into(),
        4756.into(),
        &qap,
    );

    // 4. generate proof
    let public_input = [Scalar::from(12), Scalar::from(34)];
    dbg!(&public_input);
    let private_input = [Scalar::from(234)];
    dbg!(&private_input);
    let proof = prove(&crs, qap, public_input, private_input);
    dbg!(&proof);

    // 5. verify proof
    let res1 = verify(&crs, &public_input, &proof);
    dbg!(res1);

    // let public_input_2 = [Scalar::from(111111111111), Scalar::from(999)];
    // dbg!(public_input_2);
    // let res2 = verify(&crs, &public_input_2, &proof);
    // dbg!(res2);
}
