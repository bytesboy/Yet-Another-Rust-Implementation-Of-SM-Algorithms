use std::ops::{Add, Sub};

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::FromPrimitive;

use crate::sm2::core::{BasePoint, Elliptic, Multiplication, Point};


/// 公钥 非压缩：65bytes 压缩：33bytes
#[derive(Clone, Debug)]
struct PublicKey(Point);

impl PublicKey {
    fn new(key: Point) -> Self {
        PublicKey(key)
    }

    fn from(private: &PrivateKey, base: &BasePoint) -> Self {
        PublicKey(base.multiply(private.0.clone()))
    }
}


/// 私钥 32bytes
#[derive(Clone, Debug)]
struct PrivateKey(BigUint);

impl PrivateKey {
    fn new(key: &BigUint) -> Self {
        PrivateKey(key.clone())
    }
}


/// 秘钥对（d, P）d:私钥 P:公钥
#[derive(Debug)]
struct KeyPair(PrivateKey, PublicKey);

impl KeyPair {
    pub fn private_key(&self) -> &PrivateKey {
        &self.0
    }

    pub fn public_key(&self) -> &PublicKey {
        &self.1
    }
}


/// 秘钥生成器
///
/// 使用SM2椭圆曲线公钥密码算法推荐曲线参数
///
/// 256位素数域
///
/// y^2 = x^3 + ax + b
#[derive(Clone, Debug)]
pub struct KeyGenerator {
    elliptic: Elliptic,
}

impl KeyGenerator {
    pub fn new(elliptic: &Elliptic) -> Self {
        KeyGenerator { elliptic: elliptic.clone() }
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
        // private key: 32 bytes
        let n = self.elliptic.base().order();
        let bytes: Vec<u8> = (0..self.elliptic.bits() / 8 + 8).map(|_| { rand::random::<u8>() }).collect();
        let mut k = BigUint::from_bytes_be(&bytes);
        // n-2
        let n = BigUint::sub((*n).clone(), BigUint::from_u64(2).unwrap());
        // k % n  ∈ [0, n-1]  => k % (n-2) + 1  ∈ [1, n-2]
        let key = k.mod_floor(&n).add(BigUint::from_u64(1).unwrap());
        PrivateKey::new(&key)
    }

    /// 生成公钥
    ///
    /// P = (xP,yP) = \[d]G, G为基点，d为私钥
    fn gen_public_key(&self, private_key: &PrivateKey) -> PublicKey {
        PublicKey::from(private_key, self.elliptic.base())
    }
}


#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::Num;

    use crate::sm2::core::{BasePoint, Elliptic, EllipticProvider, Point};
    use crate::sm2::key::{KeyGenerator, PrivateKey, PublicKey};
    use crate::sm2::p256::P256Elliptic;

    #[test]
    fn main() {
        let elliptic = P256Elliptic::initialize();
        let pair = KeyGenerator::new(elliptic.elliptic()).gen_key_pair();
        println!("prk = {:?}", pair.private_key());
        println!("puk = {:?}", pair.public_key());


        let elliptic = P256Elliptic::initialize();
        let pair = KeyGenerator::new(elliptic.elliptic()).gen_key_pair();
        println!("prk = {:?}", pair.private_key());
        println!("puk = {:?}", pair.public_key());
    }

    #[test]
    fn public_key() {
        // d: 48358803002808206747871163666773640956067045543241775523137833706911222329998
        // x: 76298453107918256108319614943154283626396976993715724710320433578462434588530
        // y: 22016840577845663905050918262284081863871275223913804750000840645022838962798

        let elliptic = P256Elliptic::initialize();
        let prk = "48358803002808206747871163666773640956067045543241775523137833706911222329998";
        let prk = BigUint::from_str_radix(prk, 10).unwrap();

        let private_key = PrivateKey(prk);
        let public_key = PublicKey::from(&private_key, elliptic.elliptic().base());

        println!("PriKey = {:?}", private_key);
        println!("PubKey = {:?}", public_key);
    }
}