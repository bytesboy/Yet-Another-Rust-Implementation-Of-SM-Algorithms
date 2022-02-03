use std::cmp::Ordering;
use num_bigint::BigUint;
use num_integer::Integer;


pub trait EllipticProvider {
    fn blueprint(&self) -> &Elliptic;
    /// 标量乘法
    fn scalar_multiply(&self, x: BigUint, y: BigUint, scalar: BigUint) -> (BigUint, BigUint);
    /// 基点标量乘法
    fn scalar_base_multiply(&self, scalar: BigUint) -> (BigUint, BigUint);
    /// 标量缩小：小于循环子群的阶
    fn scalar_reduce(&self, scalar: BigUint) -> BigUint {
        let elliptic = self.blueprint();
        // compare scalar and order, n = (scalar mod order) if scalar > order else scalar
        if let Ordering::Greater = scalar.cmp(&elliptic.n) {
            scalar.mod_floor(&elliptic.n)
        } else {
            scalar
        }
    }
}

/// 使用SM2椭圆曲线公钥密码算法推荐曲线参数
///
/// y^2 = x^3 + ax + b
#[derive(Clone, Debug)]
pub struct Elliptic {
    pub p: BigUint,
    pub a: BigUint,
    pub b: BigUint,
    pub gx: BigUint,
    pub gy: BigUint,
    pub n: BigUint,
    pub bits: usize,
}





