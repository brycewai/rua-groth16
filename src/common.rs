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

/////////////////////////////////////////////////////////////////////
/// 定义多项式
///
pub mod polynomial {
    use bls12_381::Scalar;
    use std::{
        cmp::Ordering,
        fmt::Debug,
        iter::{self, Sum},
        ops::{Add, Mul, MulAssign, Sub},
    };

    #[derive(Clone, Eq, PartialEq, Default)]
    pub struct Polynomial {
        pub(crate) terms: Vec<Scalar>,
    }

    /// 定义多项式的阶
    /// derive(PartialEq): 派生宏生成 trait PartialEq 的一个 impl.
    #[derive(Copy, Clone, Eq, PartialEq)]
    pub struct PolynomialDegree {
        degree: usize,
    }

    pub const ZERO: Polynomial = Polynomial { terms: Vec::new() };

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

        pub fn is_zero(&self) -> bool {
            self.terms.is_empty()
        }
    }

    /// 实现多项式之间的乘法、多项式*标量，标量*多项式，以及*=运算
    impl Mul for Polynomial {
        type Output = Polynomial;
        fn mul(self, rhs: Self) -> Self::Output {
            if self.is_zero() {
                return Polynomial::zero();
            }
            // Polynomial::zero()
            let degree_1 = self.get_degree().degree_expect_zero_polynomial().unwrap();
            let degree_2 = rhs.get_degree().degree_expect_zero_polynomial().unwrap();
            let degree = degree_1 + degree_2;
            let mut p = Polynomial {
                terms: vec![Scalar::zero(); degree + 1],
            };
            for i in 0..=degree {
                for l_index in i.saturating_sub(degree_2)..=degree_1.min(i) {
                    let r_index = i - l_index;
                    p.terms[i] += self.terms[l_index] * rhs.terms[r_index];
                }
            }
            p
        }
    }

    impl Mul<Scalar> for Polynomial {
        type Output = Polynomial;
        fn mul(mut self, rhs: Scalar) -> Self::Output {
            for i in self.terms.iter_mut() {
                *i *= rhs;
            }
            self
        }
    }

    impl MulAssign<Polynomial> for Polynomial {
        fn mul_assign(&mut self, rhs: Polynomial) {
            *self = self.clone() * rhs
        }
    }

    impl MulAssign<Scalar> for Polynomial {
        fn mul_assign(&mut self, rhs: Scalar) {
            *self = self.clone() * rhs
        }
    }

    // 标量 * 多项式 = 多项式 * 标量
    impl Mul<Polynomial> for Scalar {
        type Output = Polynomial;
        fn mul(self, rhs: Polynomial) -> Self::Output {
            // 反过来就是死循环，实现错误
            rhs * self
        }
    }

    /// 实现多项式加法、减法、求和
    impl Sum for Polynomial {
        fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
            let mut res = Polynomial::zero();
            for i in iter {
                res = res + i;
            }
            res
        }
    }

    impl Add for Polynomial {
        type Output = Polynomial;
        fn add(mut self, rhs: Self) -> Self::Output {
            for (coef, other) in self.terms.iter_mut().zip(rhs.terms) {
                *coef += other;
            }
            self
        }
    }

    impl Sub for Polynomial {
        type Output = Polynomial;
        fn sub(mut self, rhs: Self) -> Self::Output {
            for (coef, other) in self.terms.iter_mut().zip(rhs.terms) {
                *coef -= other;
            }
            self
        }
    }

    /// 实现多项式的Debug trait
    impl Debug for Polynomial {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            // enumerate: 创建一个迭代器，该迭代器给出当前迭代次数以及下一个值: (i, val)
            // rev: 通常，迭代器从左到右进行迭代。 使用 rev() 之后，迭代器将改为从右向左进行迭代。
            for (i, coef) in self.terms.iter().enumerate().rev() {
                match i {
                    0 => {
                        write!(f, "{coef}")?;
                    }
                    1 => {
                        write!(f, "{coef} * x + ")?;
                    }
                    i => {
                        write!(f, "{coef} * x^{i} + ")?;
                    }
                }
            }
            Ok(())
        }
    }
}

/////////////////////////////////////////////////////////////////////
/// 拉格朗日插值法
///

pub mod lagrange {
    use super::polynomial::{Polynomial, ZERO};
    use bls12_381::Scalar;

    /// 创建一个N阶的基础多项式，使得第i个多项式在i点为1，其他点为0
    pub fn basis_polynomials<const N: usize>(evaluation_points: [Scalar; N]) -> [Polynomial; N] {
        let mut result = [ZERO; N];
        for p in result.iter_mut() {
            *p = Polynomial::identity();
        }
        for i in 0..N {
            let numerator_polynomial = &mut result[i];
            let mut one = Scalar::one();
            for j in 0..N {
                if i != j {
                    *numerator_polynomial *= Polynomial::x_minus_const(evaluation_points[j]);
                    one *= evaluation_points[i] - evaluation_points[j]
                }
            }
            *numerator_polynomial *= one.invert().unwrap();
        }
        result
    }
    /// 拉格朗日插值法
    pub fn interpolate<const N: usize>(
        basis_polynomials: &[Polynomial; N],
        evaluations: &[Scalar; N],
    ) -> Polynomial {
        basis_polynomials
            .iter()
            .zip(evaluations.iter())
            .map(|(p, v)| *v * p.clone())
            .sum()
    }
}

/////////////////////////////////////////////////////////////////////
/// 单位根
///
pub mod roots_of_unity {
    use bls12_381::Scalar;
    use ff::PrimeField;
    // (2m)^n
    pub fn primitive_root_of_unity(mut n: usize) -> Scalar {
        assert!(n.is_power_of_two());
        assert!(n <= (1 << 32));
        let mut root_of_unity = Scalar::ROOT_OF_UNITY;
        while n < (1 << 32) {
            n <<= 1;
            root_of_unity = root_of_unity.square();
        }
        root_of_unity
    }

    pub fn root_of_unity<const N: usize>() -> [Scalar; N] {
        let root_of_unity = primitive_root_of_unity(N);
        let mut all_roots = [Scalar::one(); N];
        all_roots[0] = root_of_unity;
        for i in 0..N {
            all_roots[i] = all_roots[i - 1] * root_of_unity;
        }
        all_roots
    }
}

/////////////////////////////////////////////////////////////////////
/// 单元测试
///
#[cfg(test)]
mod tests {
    use bls12_381::Scalar;

    use crate::common::{lagrange, polynomial::Polynomial, roots_of_unity::root_of_unity};

    #[test]
    fn test_evaluate() {
        // f(x) = x - 1
        let p = Polynomial::x_minus_const(Scalar::one());
        assert_eq!(p.evaluation(&Scalar::zero()), -Scalar::one());
        assert_eq!(p.evaluation(&Scalar::one()), Scalar::zero());
    }

    #[test]
    fn test_mul_zero() {
        let p1 = Polynomial::zero();
        let p2 = Polynomial {
            terms: vec![0.into(), 3.into(), 5.into()],
        };
        assert_eq!(p1.clone() * p2.clone(), Polynomial::zero());
    }

    #[test]
    fn test_lagrange_basis_polynomials() {
        const N: usize = 2 << 4;
        let roots = root_of_unity::<N>();
        let lagrange_basis_polynomials = lagrange::basis_polynomials(roots);

        for (i, p) in lagrange_basis_polynomials.iter().enumerate() {
            for (j, evaluation_point) in roots.iter().enumerate() {
                let val = p.evaluation(evaluation_point);
                if i == j {
                    assert_eq!(val, Scalar::one());
                } else {
                    assert_eq!(val, Scalar::zero());
                }
            }
        }
    }
}
