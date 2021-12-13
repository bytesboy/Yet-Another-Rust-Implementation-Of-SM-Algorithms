use num_bigint::BigUint;


#[derive(Clone, Debug)]
pub struct CurveParams {
    // 大素数，有限域{0,1,2,...,p-1}
    pub p: BigUint,
    // 椭圆曲线参数a
    pub a: BigUint,
    // 椭圆曲线参数b
    pub b: BigUint,
    // 循环子域的阶
    pub n: BigUint,
    // 椭圆曲线基点G
    pub g: (BigUint, BigUint),
}