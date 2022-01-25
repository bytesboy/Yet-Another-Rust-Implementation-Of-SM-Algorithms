use std::ops::{Add, Sub};
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::FromPrimitive;
use crate::sm2::core::EllipticProvider;

/// 公钥 非压缩：65bytes 压缩：33bytes
#[derive(Clone, Debug)]
pub struct PublicKey(BigUint, BigUint);

/// 私钥 32bytes
#[derive(Clone, Debug)]
pub struct PrivateKey(BigUint);

/// 秘钥对（d, P）d:私钥 P:公钥
#[derive(Debug)]
pub struct KeyPair(PrivateKey, PublicKey);

impl KeyPair {
    pub fn private_key(&self) -> &PrivateKey {
        &self.0
    }
    pub fn public_key(&self) -> &PublicKey {
        &self.1
    }
}

/// 秘钥生成器
pub struct KeyGenerator {
    elliptic: Box<dyn EllipticProvider>,
}

impl KeyGenerator {
    pub fn new(elliptic: Box<dyn EllipticProvider>) -> Self {
        KeyGenerator { elliptic }
    }

    pub fn gen_key_pair(&self) -> KeyPair {
        let private_key = self.gen_private_key();
        let public_key = self.gen_public_key(&private_key);
        KeyPair(private_key.clone(), public_key.clone())
    }

    /// 生成私钥
    ///
    /// 随机数生成整数d ∈ \[1, n − 2], n为椭圆曲线循环子群的阶
    fn gen_private_key(&self) -> PrivateKey {
        // private key: 32 bytes
        let elliptic = self.elliptic.blueprint();
        let n = &elliptic.n;
        let bytes: Vec<u8> = (0..elliptic.bits / 8 + 8).map(|_| { rand::random::<u8>() }).collect();
        let mut k = BigUint::from_bytes_be(&bytes);
        // n-2
        let n = BigUint::sub((*n).clone(), BigUint::from_u64(2).unwrap());
        // k % n  ∈ [0, n-1]  => k % (n-2) + 1  ∈ [1, n-2]
        let key = k.mod_floor(&n).add(BigUint::from_u64(1).unwrap());
        PrivateKey(key)
    }

    /// 生成公钥
    ///
    /// P = (xP,yP) = \[d]G, G为基点，d为私钥
    fn gen_public_key(&self, private_key: &PrivateKey) -> PublicKey {
        let key = self.elliptic.scalar_base_multiply(private_key.0.clone());
        PublicKey(key.0, key.1)
    }
}


#[cfg(test)]
mod tests {
    use num_traits::Num;
    use crate::sm2::p256::P256Elliptic;
    use super::*;

    #[test]
    fn main() {
        let generator = KeyGenerator::new(Box::new(P256Elliptic::init()));
        let pair = generator.gen_key_pair();
        println!("prk = {:?}", pair.private_key());
        println!("puk = {:?}", pair.public_key());
    }

    #[test]
    fn public_key() {
        // d: 48358803002808206747871163666773640956067045543241775523137833706911222329998
        // x: 76298453107918256108319614943154283626396976993715724710320433578462434588530
        // y: 22016840577845663905050918262284081863871275223913804750000840645022838962798

        let prk = "48358803002808206747871163666773640956067045543241775523137833706911222329998";
        let prk = BigUint::from_str_radix(prk, 10).unwrap();

        let private_key = PrivateKey(prk);
        let generator = KeyGenerator::new(Box::new(P256Elliptic::init()));
        let public_key = generator.gen_public_key(&private_key);

        println!("PriKey = {:?}", private_key);
        println!("PubKey = {:?}", public_key);
    }
}