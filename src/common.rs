#[warn(unused_imports)]
use std::{
    cmp::Ordering,
    fmt::Debug,
    iter::{self, Sum},
    ops::{Add, Mul, MulAssign, Sub},
};

use bls12_381::Scalar;

use pairing::PairingCurveAffine;

/// Get an array with powers of x on a curve, i.e. obtain G, x*G, x^2*G, etc.
fn powers_of_x_on_curve<
    const N: usize,
    Curve: PairingCurveAffine
        + std::ops::Mul<bls12_381::Scalar>
        + std::convert::From<<Curve as std::ops::Mul<bls12_381::Scalar>>::Output>,
>(
    x: Scalar,
    generator: Curve,
) -> [Curve; N] {
    // let generator = Curve::generator();
    let mut powers = [generator; N];
    for i in 1..N {
        powers[i] = (powers[i - 1] * x).into();
    }

    powers
}

// fn power_of_x<
//     const N: usize,
//     Curve: PairingCurveAffine
//         + std::ops::Mul<bls12_381::Scalar>
//         + std::convert::From<<Curve as std::ops::Mul<bls12_381>>::Output>,
// >(
//     x: Scalar,
// ) -> [Curve; N] {
//     powers_of_x_on_curve(x, Curve::generator())
// }

/// 定义多项式
pub struct Polynomial {
    terms: Vec<Scalar>,
}

/// 定义多项式的阶
/// derive(PartialEq): 派生宏生成 trait PartialEq 的一个 impl.
#[derive(Copy, Clone, Eq, PartialEq)]
pub struct PolynomialDegree {
    degree: usize,
}

/// 参考 https://rustwiki.org/zh-CN/std/cmp/trait.Ord.html
impl Ord for PolynomialDegree {
    fn cmp(&self, other: &Self) -> Ordering {
        self.degree.cmp(&other.degree)
    }
}

/// 参考 https://rustwiki.org/zh-CN/std/cmp/trait.Ord.html
impl PartialOrd for PolynomialDegree {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// 实现多项式阶包含的方法
impl PolynomialDegree {
    pub fn is_zero_polynomial(&self) -> bool {
        self.degree == 0 || self.degree == usize::MAX
    }

    pub fn degree_expect_zero_polynomial(&self) -> Option<usize> {
        Some(self.degree).filter(|deg| *deg != 0 && *deg != usize::MAX)
    }
}

/// 实现多项式包含的方法
impl Polynomial {
    /// setup 中使用到的 t(x) = （x-r1)(x-r2)(x-r3)...(x-rn)
    /// repeat() 函数一次又一次地重复单个值。无限迭代器 (如 repeat()) 通常与适配器 (如 Iterator::take()) 一起使用，以使其具有有限性。
    /// 多项式的首位加上 -Scalar::one 和 Scalar::one 是为什么？
    pub fn t_x(degree: usize) -> Self {
        Polynomial {
            terms: iter::once(-Scalar::one())
                .chain(iter::repeat(Scalar::zero()).take(degree - 1))
                .chain(iter::once(Scalar::one()))
                .collect(),
        }
    }

    /// 零多项式
    pub fn zero() -> Self {
        Polynomial { terms: Vec::new() }
    }

    /// 恒等多项式，恒等于1
    pub fn identity() -> Self {
        Polynomial {
            terms: vec![Scalar::one()],
        }
    }

    /// 多项式：x - const
    pub fn x_minus_const(_const: Scalar) -> Self {
        Polynomial {
            terms: vec![-_const, Scalar::one()],
        }
    }

    /// 获取指定多项式的阶
    pub fn get_degree(&self) -> PolynomialDegree {
        PolynomialDegree {
            degree: self.terms.len().wrapping_sub(1),
        }
    }

    /// 多项式在指定点的求值计算
    /// rev 反转迭代器的方向
    /// 通过应用操作将每个元素 fold 到一个累加器中，返回最终结果
    pub fn evaluation(&self, x: &Scalar) -> Scalar {
        self.terms
            .iter()
            .rev()
            .fold(Scalar::zero(), |value, term| value * x + term)
    }
}
