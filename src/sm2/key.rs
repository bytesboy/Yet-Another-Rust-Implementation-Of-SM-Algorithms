use std::ops::Sub;

use num_bigint::BigUint;
use num_traits::{Num, One};

use crate::sm2::ecc::EllipticBuilder;

pub trait HexKey {
    fn encode(&self) -> String;
    fn decode(key: &str) -> Self;
}

/// 公钥
/// 非压缩公钥格式字节串长度为65字节，压缩格式长度为33字节;
/// 非压缩格式公钥首字节为0x04。
/// 压缩格式公钥，若公钥y坐标最后一位为0，则首字节为0x02，否则为0x03。
/// 签名长度：64字节。
#[derive(Clone, Debug)]
pub struct PublicKey(BigUint, BigUint);

impl PublicKey {
    pub fn value(&self) -> (BigUint, BigUint) {
        (self.0.clone(), self.1.clone())
    }
}

impl HexKey for PublicKey {
    fn encode(&self) -> String {
        let key = {
            let x = self.0.to_bytes_be();
            let y = self.1.to_bytes_be();
            [vec![0x04], to_32_bytes(x).to_vec(), to_32_bytes(y).to_vec()].concat()
        };
        hex::encode(key)
    }

    fn decode(key: &str) -> Self {
        if key.len() != 130 {
            panic!("The uncompressed public key's length must be 130.")
        }

        if !key.starts_with("04") {
            panic!("The compressed public key is invalid.")
        }

        let key = match hex::decode(key.trim_start_matches("04")) {
            Ok(data) => data,
            Err(_) => panic!("The public key must be composed of hex chars.")
        };

        PublicKey(
            BigUint::from_bytes_be(&key[..32]),
            BigUint::from_bytes_be(&key[32..]),
        )
    }
}


/// 私钥 32bytes
#[derive(Clone, Debug)]
pub struct PrivateKey(BigUint);

impl PrivateKey {
    pub fn value(&self) -> BigUint {
        self.0.clone()
    }
}

impl HexKey for PrivateKey {
    fn encode(&self) -> String {
        hex::encode(to_32_bytes(self.0.to_bytes_be()))
    }

    fn decode(key: &str) -> Self {
        if key.len() != 64 {
            panic!("The length of the private key must be 64.")
        }
        let key = match BigUint::from_str_radix(&*key, 16) {
            Ok(data) => data,
            Err(_) => panic!("The private key must be composed of hex chars.")
        };
        PrivateKey(key)
    }
}


/// 秘钥对（d, P）d:私钥 P:公钥
#[derive(Debug)]
pub struct KeyPair(PrivateKey, PublicKey);

impl KeyPair {
    pub fn new(prk: PrivateKey, puk: PublicKey) -> Self {
        KeyPair(prk, puk)
    }
    pub fn prk(&self) -> &PrivateKey {
        &self.0
    }
    pub fn puk(&self) -> &PublicKey {
        &self.1
    }
}

/// 秘钥生成器
pub struct KeyGenerator {
    builder: Box<dyn EllipticBuilder>,
}

impl KeyGenerator {
    pub fn init(builder: Box<dyn EllipticBuilder>) -> Self {
        KeyGenerator { builder }
    }

    pub fn gen_key_pair(&self) -> KeyPair {
        let private_key = self.gen_private_key();
        let public_key = self.gen_public_key(&private_key);
        KeyPair(private_key.clone(), public_key.clone())
    }

    /// 生成私钥
    ///
    /// d ∈ \[1, n − 2]
    fn gen_private_key(&self) -> PrivateKey {
        let e = self.builder.blueprint();
        let from = BigUint::one();
        let to = e.n.clone().sub(BigUint::from(2u8));
        PrivateKey(e.random(from, to))
    }

    /// 生成公钥
    ///
    /// P = (x,y) = dG, G为基点，d为私钥
    fn gen_public_key(&self, private_key: &PrivateKey) -> PublicKey {
        let key = self.builder.scalar_base_multiply(private_key.value());
        PublicKey(key.0, key.1)
    }
}

#[inline(always)]
pub fn copy_slice(dst: &mut [u8], src: &[u8]) {
    for (d, s) in dst.iter_mut().zip(src.iter()) {
        *d = *s;
    }
}

#[inline(always)]
pub fn to_32_bytes(data: Vec<u8>) -> [u8; 32] {
    let mut result = [0u8; 32];
    if data.len() > 32 {
        copy_slice(&mut result, &data[(data.len() - 32)..]);
    } else if data.len() < 32 {
        copy_slice(&mut result[(32 - data.len())..], &data);
    } else {
        copy_slice(&mut result, &data);
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::sm2::p256::P256Elliptic;

    use super::*;

    #[test]
    fn main() {
        let generator = KeyGenerator::init(Box::new(P256Elliptic::init()));
        let pair = generator.gen_key_pair();
        println!("prk = {:?}", pair.prk());
        println!("puk = {:?}", pair.puk());
    }

    #[test]
    fn generator() {
        // d: 48358803002808206747871163666773640956067045543241775523137833706911222329998
        // x: 76298453107918256108319614943154283626396976993715724710320433578462434588530
        // y: 22016840577845663905050918262284081863871275223913804750000840645022838962798

        let prk = "48358803002808206747871163666773640956067045543241775523137833706911222329998";
        let prk = BigUint::from_str_radix(prk, 10).unwrap();

        let private_key = PrivateKey(prk);
        let generator = KeyGenerator::init(Box::new(P256Elliptic::init()));
        let public_key = generator.gen_public_key(&private_key);

        assert_eq!(private_key.0.to_string(), "48358803002808206747871163666773640956067045543241775523137833706911222329998");
        assert_eq!(public_key.0.to_string(), "76298453107918256108319614943154283626396976993715724710320433578462434588530");
        assert_eq!(public_key.1.to_string(), "22016840577845663905050918262284081863871275223913804750000840645022838962798");

        assert_eq!(private_key.encode(), "6aea1ccf610488aaa7fddba3dd6d76d3bdfd50f957d847be3d453defb695f28e");
        assert_eq!(public_key.encode(), "04a8af64e38eea41c254df769b5b41fbaa2d77b226b301a2636d463c52b46c777230ad1714e686dd641b9e04596530b38f6a64215b0ed3b081f8641724c5443a6e");
    }

    #[test]
    fn key() {
        let prk = "6aea1ccf610488aaa7fddba3dd6d76d3bdfd50f957d847be3d453defb695f28e";
        let puk = "04a8af64e38eea41c254df769b5b41fbaa2d77b226b301a2636d463c52b46c777230ad1714e686dd641b9e04596530b38f6a64215b0ed3b081f8641724c5443a6e";

        let private_key = PrivateKey::decode(prk);
        assert_eq!(private_key.0.to_string(), "48358803002808206747871163666773640956067045543241775523137833706911222329998");

        let public_key = PublicKey::decode(puk);
        assert_eq!(public_key.0.to_string(), "76298453107918256108319614943154283626396976993715724710320433578462434588530");
        assert_eq!(public_key.1.to_string(), "22016840577845663905050918262284081863871275223913804750000840645022838962798");
    }
}