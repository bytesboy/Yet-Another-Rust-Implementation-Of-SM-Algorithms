use std::ops::{Add, Sub};

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{FromPrimitive, Num};

use crate::sm2::ecc::EllipticBuilder;

pub trait HexKey {
    fn encode(&self) -> String;
    fn decode(key: String) -> Self;
}

// pub trait ZipHexKey {
//     fn zip(&self) -> String;
//     fn unzip(&self, key: String, elliptic: Box<dyn EllipticProvider>) -> Self;
// }

/// 公钥
/// 非压缩公钥格式字节串长度为65字节，压缩格式长度为33字节;
/// 非压缩格式公钥首字节为0x04。
/// 压缩格式公钥，若公钥y坐标最后一位为0，则首字节为0x02，否则为0x03。
/// 签名长度：64字节。
#[derive(Clone, Debug)]
pub struct PublicKey(BigUint, BigUint);

impl HexKey for PublicKey {
    fn encode(&self) -> String {
        let (x, y) = (self.0.to_bytes_be(), self.1.to_bytes_be());
        if x.len() != 32 || y.len() != 32 {
            panic!("The public key is error, x and y coordinates must be 32 bytes.")
        }
        let key_bytes = {
            // key = 0x04 || x || y
            let mut key = vec![0u8; 65];
            key[0] = 0x04;
            for (dst, src) in key[1..33].iter_mut().zip(x.iter()) {
                *dst = *src;
            }
            for (dst, src) in key[33..].iter_mut().zip(y.iter()) {
                *dst = *src;
            }
            key
        };
        hex::encode(key_bytes)
    }

    fn decode(key: String) -> Self {
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

// impl ZipHexKey for PublicKey {
//     fn zip(&self) -> String {
//         let (x, y) = (self.0.to_bytes_be(), self.1.to_bytes_be());
//         if x.len() != 32 || y.len() != 32 {
//             panic!("The public key is error, x and y coordinates must be 32 bytes.")
//         }
//
//         let key_bytes = {
//             // key = 0x02 or 0x03 || x
//             let mut key = vec![0u8; 33];
//             if BigUint::zero() == self.1.clone().bitand(BigUint::from(1u8)) {
//                 key[0] = 0x02;
//             } else {
//                 key[1] = 0x03;
//             }
//             for (dst, src) in key[1..].iter_mut().zip(x.iter()) {
//                 *dst = *src;
//             }
//             key
//         };
//
//         hex::encode(key_bytes)
//     }
//
//     fn unzip(&self, key: String, elliptic: Box<dyn EllipticProvider>) -> Self {
//         if key.len() != 66 {
//             panic!("The compressed public key's length must be 66.")
//         }
//
//         if !(key.starts_with("02") && key.starts_with("03")) {
//             panic!("The compressed public key is invalid.")
//         }
//
//         let x = match hex::decode(key.split_at(2).1) {
//             Ok(data) => BigUint::from_bytes_be(data.as_slice()),
//             Err(_) => panic!("The public key must be composed of hex chars.")
//         };
//     }
// }

/// 私钥 32bytes
#[derive(Clone, Debug)]
pub struct PrivateKey(BigUint);

impl HexKey for PrivateKey {
    fn encode(&self) -> String {
        let key_bytes = self.0.to_bytes_be();
        if key_bytes.len() != 32 {
            panic!("The private key is error, it must be 32 bytes.")
        }
        hex::encode(key_bytes)
    }

    fn decode(key: String) -> Self {
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
    pub fn private_key(&self) -> &PrivateKey {
        &self.0
    }
    pub fn public_key(&self) -> &PublicKey {
        &self.1
    }
}

/// 秘钥生成器
pub struct KeyGenerator {
    builder: Box<dyn EllipticBuilder>,
}

impl KeyGenerator {
    pub fn new(builder: Box<dyn EllipticBuilder>) -> Self {
        KeyGenerator { builder }
    }

    pub fn gen_key_pair(&self) -> KeyPair {
        let private_key = self.gen_private_key();
        let public_key = self.gen_public_key(&private_key);
        KeyPair(private_key.clone(), public_key.clone())
    }

    /// 生成私钥
    fn gen_private_key(&self) -> PrivateKey {
        // private key: 32 bytes
        let elliptic = self.builder.blueprint();
        // subgroup order
        let n = &elliptic.n;
        let k = elliptic.random();
        // n-2
        let n = BigUint::sub((*n).clone(), BigUint::from_u64(2).unwrap());
        // k % n  ∈ [0, n-1]  => k % (n-2) + 1  ∈ [1, n-2] => key ∈ [1, n-1)
        let key = k.mod_floor(&n).add(BigUint::from_u64(1).unwrap());
        PrivateKey(key)
    }

    /// 生成公钥
    ///
    /// P = (x,y) = dG, G为基点，d为私钥
    fn gen_public_key(&self, private_key: &PrivateKey) -> PublicKey {
        let key = self.builder.scalar_base_multiply(private_key.0.clone());
        PublicKey(key.0, key.1)
    }
}


#[cfg(test)]
mod tests {
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
    fn generator() {
        // d: 48358803002808206747871163666773640956067045543241775523137833706911222329998
        // x: 76298453107918256108319614943154283626396976993715724710320433578462434588530
        // y: 22016840577845663905050918262284081863871275223913804750000840645022838962798

        let prk = "48358803002808206747871163666773640956067045543241775523137833706911222329998";
        let prk = BigUint::from_str_radix(prk, 10).unwrap();

        let private_key = PrivateKey(prk);
        let generator = KeyGenerator::new(Box::new(P256Elliptic::init()));
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

        let private_key = PrivateKey::decode(String::from(prk));
        assert_eq!(private_key.0.to_string(), "48358803002808206747871163666773640956067045543241775523137833706911222329998");

        let public_key = PublicKey::decode(String::from(puk));
        assert_eq!(public_key.0.to_string(), "76298453107918256108319614943154283626396976993715724710320433578462434588530");
        assert_eq!(public_key.1.to_string(), "22016840577845663905050918262284081863871275223913804750000840645022838962798");
    }
}