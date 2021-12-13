use std::mem;
use std::ops::{Add, Sub};
use std::sync::Once;
use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use crate::sm2::core::CurveParams;
use crate::sm2::p256::{BIT_SIZE, EC_A, EC_B, EC_GX, EC_GY, EC_N, EC_P};

/// 公钥 64bytes
#[derive(Clone, Debug)]
struct PublicKey(BigUint, BigUint);

impl PublicKey {
    pub fn new(private: &PrivateKey, base: &(BigUint, BigUint)) -> Self {
        PublicKey((*base).0.clone(), (*base).1.clone())
    }
}

/// 私钥 32bytes
#[derive(Clone, Debug)]
struct PrivateKey(BigUint);

impl PrivateKey {
    pub fn new(n: &BigUint) -> Self {
        // private key: 32 bytes
        let bytes: Vec<u8> = (0..BIT_SIZE / 8 + 8).map(|_| { rand::random::<u8>() }).collect();
        let mut k = BigUint::from_bytes_be(&bytes);
        // n-2
        let n = BigUint::sub((*n).clone(), BigUint::from_u64(2).unwrap());
        // k % n  ∈ [0, n-1]  => k % (n-2) + 1  ∈ [1, n-2]
        let key = k.mod_floor(&n).add(BigUint::from_u64(1).unwrap());
        PrivateKey(key)
    }
}

/// 秘钥对（d, P）d:私钥 P:公钥
#[derive(Debug)]
struct KeyPair(PrivateKey, PublicKey);


/// 秘钥生成器
///
/// 使用SM2椭圆曲线公钥密码算法推荐曲线参数
///
/// 256位素数域
///
/// y^2 = x^3 + ax + b
#[derive(Clone, Debug)]
pub struct KeyGenerator {
    curve: CurveParams,
}

impl KeyGenerator {
    pub fn new() -> Self {
        static mut GENERATOR: *const KeyGenerator = 0 as *const KeyGenerator;
        static ONCE: Once = Once::new();
        unsafe {
            ONCE.call_once(|| {
                let mut generator = KeyGenerator {
                    curve: CurveParams {
                        p: BigUint::from_bytes_be(&EC_P),
                        a: BigUint::from_bytes_be(&EC_A),
                        b: BigUint::from_bytes_be(&EC_B),
                        n: BigUint::from_bytes_be(&EC_N),
                        g: (BigUint::from_bytes_be(&EC_GX), BigUint::from_bytes_be(&EC_GY)),
                    }
                };
                GENERATOR = mem::transmute(Box::new(generator));
            });
            (*GENERATOR).clone()
        }
    }

    /// 输入:一个有效的Fq(q = p且p为大于3的素数，或q = 2m)上椭圆曲线系统参数的集合。
    ///
    /// 输出:与椭圆曲线系统参数相关的一个密钥对(d , P)。
    ///
    /// a) 用随机数发生器产生整数d ∈ \[1, n − 2];
    ///
    /// b) G为基点，计算点P = (xP,yP) = \[d]G;
    ///
    /// c) 密钥对是(d,P)，其中d为私钥，P为公钥。
    pub fn gen_key_pair(&self) -> KeyPair {
        let private_key = self.gen_private_key();
        KeyPair(private_key.clone(), self.gen_public_key(&private_key).clone())
    }

    /// 生成私钥
    ///
    /// 随机数生成整数d ∈ \[1, n − 2], n为椭圆曲线循环子群的阶
    fn gen_private_key(&self) -> PrivateKey {
        PrivateKey::new(&self.curve.n)
    }

    /// 生成公钥
    ///
    /// P = (xP,yP) = \[d]G, G为基点，d为私钥
    fn gen_public_key(&self, private_key: &PrivateKey) -> PublicKey {
        PublicKey::new(private_key, &self.curve.g)
    }
}


#[cfg(test)]
mod tests {
    use crate::sm2::key::KeyGenerator;

    #[test]
    fn main() {
        let prk = KeyGenerator::new().gen_private_key();
        println!("prk = {:?}", prk);
    }
}