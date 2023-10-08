use crate::common::{self, Polynomial};
use bls12_381::Scalar;

/// 定义约束矩阵
pub struct ConstraintMatrix<
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
    const CONSTRAINTS: usize,
> {
    public_entries: [[Scalar; CONSTRAINTS]; PUBLIC_WITNESS], // 2维数组
    private_entries: [[Scalar; CONSTRAINTS]; PRIVATE_WITNESS], // 2维数组
}

/// R1CS
/// R1CS是约束系统，是约束的合集。
/// 每个逻辑门约束的形式为： (a DOT x) * (b DOT x) = (c DOT x)
/// 其中 a、b 和 c 是系数向量（由以某个素数 p 为模的整数组成的素数域的元素），x 是不同“伪变量”的向量（解向量）。
/// 每个伪变量要么是一个变量，表示域的一个元素，要么是特殊符号 1，表示域元素 1。
/// DOT 表示取两个向量的点积，除了所有加法和乘法都是以 p 为模进行的， * 表示两个标量模 p 的乘积。使用伪变量 1 允许用点积表示常数加数。
pub struct R1cs<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
{
    a: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    b: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    c: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
}

impl<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
    R1cs<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>
{
    pub fn new(
        a: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
        b: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
        c: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    ) -> Self {
        R1cs { a, b, c }
    }

    // 验证给定的r1cs中的public_witness和private_witness是否满足等式
    pub fn is_correct(
        &self,
        public_witness: [Scalar; PUBLIC_WITNESS],
        private_witness: [Scalar; PRIVATE_WITNESS],
    ) -> bool {
        // 矩阵运算，验证每个约束是否满足： (a DOT x) * (b DOT x) = (c DOT x)
        // CONSTRAINTS常量个约束就需要验证CONSTRAINTS次
        for constraint_index in 0..CONSTRAINTS {
            let mut a_sum = Scalar::zero();
            let mut b_sum = Scalar::zero();
            let mut c_sum = Scalar::zero();
            // 创建一个迭代器，该迭代器给出当前迭代次数以及下一个值: (i, val)
            for (private_index, private_witness_entry) in private_witness.iter().enumerate() {
                a_sum +=
                    private_witness_entry * self.a.private_entries[private_index][constraint_index];
                b_sum +=
                    private_witness_entry * self.b.private_entries[private_index][constraint_index];
                c_sum +=
                    private_witness_entry * self.c.private_entries[private_index][constraint_index];
            }

            for (public_index, public_witness_entry) in public_witness.iter().enumerate() {
                a_sum +=
                    public_witness_entry * self.a.public_entries[public_index][constraint_index];
                b_sum +=
                    public_witness_entry * self.b.public_entries[public_index][constraint_index];
                c_sum +=
                    public_witness_entry * self.c.public_entries[public_index][constraint_index];
            }

            if a_sum * b_sum != c_sum {
                return false;
            }
        }
        return true;
    }
}

/// QAP
#[derive(Clone)]
pub struct Qap<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize> {
    pub public_u: [Polynomial; PUBLIC_WITNESS],
    pub public_v: [Polynomial; PUBLIC_WITNESS],
    pub public_w: [Polynomial; PUBLIC_WITNESS],

    pub private_u: [Polynomial; PRIVATE_WITNESS],
    pub private_v: [Polynomial; PRIVATE_WITNESS],
    pub private_w: [Polynomial; PRIVATE_WITNESS],
}

/// 通过 R1CS 生成 QAP （From trait 允许一种类型定义 “怎么根据另一种类型生成自己”，因此它提供了一种类型转换的简单机制。）
impl<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
    from<R1cs<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>>
    for Qap<PUBLIC_WITNESS, PRIVATE_WITNESS>
{
    fn from(r1cs: R1cs<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>) -> Self {
        let mut qap = Qap {
            public_u: [Polynomial::zero(); PUBLIC_WITNESS],
            public_v: [Polynomial::zero(); PUBLIC_WITNESS],
            public_w: [Polynomial::zero(); PUBLIC_WITNESS],

            private_u: [Polynomial::zero(); PRIVATE_WITNESS],
            private_v: [Polynomial::zero(); PRIVATE_WITNESS],
            private_w: [Polynomial::zero(); PRIVATE_WITNESS],
        };
        // 拉格朗日
    }
}
