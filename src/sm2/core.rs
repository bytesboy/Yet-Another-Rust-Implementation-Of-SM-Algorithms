use num_bigint::BigUint;


pub trait EllipticProvider {
    fn blueprint(&self) -> &Elliptic;
    /// 标量乘法
    fn scalar_multiply(&self, x: BigUint, y: BigUint, k: BigUint) -> (BigUint, BigUint);
    /// 基点标量乘法
    fn scalar_base_multiply(&self, k: BigUint) -> (BigUint, BigUint);
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





