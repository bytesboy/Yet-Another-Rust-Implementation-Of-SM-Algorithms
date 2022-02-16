use std::cmp::Ordering;
use std::rc::Rc;

use num_bigint::BigUint;
use num_integer::Integer;

use crate::sm2::key::{PrivateKey, PublicKey};
use crate::sm3;

pub trait EllipticBuilder {
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

impl Elliptic {
    pub fn random(&self) -> BigUint {
        let bytes: Vec<u8> = (0..self.bits / 8 + 8).map(|_| { rand::random::<u8>() }).collect();
        BigUint::from_bytes_be(&bytes)
    }
}

#[derive(Debug, Copy, Clone)]
pub enum Mode {
    C1C2C3,
    C1C3C2,
}

pub struct Crypto {
    mode: Mode,
    builder: Rc<dyn EllipticBuilder>,
}

impl Crypto {
    pub fn init(mode: Mode, builder: Rc<dyn EllipticBuilder>) -> Self {
        Crypto { mode, builder }
    }

    pub fn encryptor(&self, key: PublicKey) -> Encryptor {
        Encryptor { key, mode: self.mode, builder: self.builder.clone() }
    }

    pub fn decryptor(&self, key: PrivateKey) -> Decryptor {
        Decryptor { key, mode: self.mode, builder: self.builder.clone() }
    }
}


pub trait Encryption {
    fn execute(&self, plain: &str) -> String;
}

pub trait Decryption {
    fn execute(&self, cipher: &str) -> String;
}

pub struct Encryptor {
    mode: Mode,
    key: PublicKey,
    builder: Rc<dyn EllipticBuilder>,
}


impl Encryption for Encryptor {
    fn execute(&self, plain: &str) -> String {
        let data = plain.as_bytes();
        let key = self.key.value();
        let cipher = loop {
            let k = self.builder.blueprint().random();
            // C1: [k]G
            let c1 = {
                let (x1, y1) = self.builder.scalar_base_multiply(k.clone());
                [vec![0x04], x1.to_bytes_be(), y1.to_bytes_be()].concat()
            };

            let (x2, y2) = self.builder.scalar_multiply(key.0.clone(), key.1.clone(), k.clone());
            let temp = [x2.to_bytes_be(), y2.to_bytes_be()].concat();
            let t = kdf(temp, data.len());

            if is_all_zero(t.clone()) {
                continue;
            }

            // C2: M ^ KDF(x2 ‖ γ2, len(M))
            let mut c2 = vec![];
            for i in 0..data.len() {
                c2.push(data[i] ^ t.clone()[i]);
            }

            // C3: hash(x2 ‖ M ‖ γ2)
            let c3 = {
                let data = [x2.to_bytes_be(), data.to_vec(), y2.to_bytes_be()].concat();
                sm3::hash(data.as_slice()).to_vec()
            };

            break match self.mode {
                Mode::C1C3C2 => [c1, c3, c2].concat(),
                Mode::C1C2C3 => [c1, c2, c3].concat()
            };
        };

        hex::encode(cipher)
    }
}

pub struct Decryptor {
    mode: Mode,
    key: PrivateKey,
    builder: Rc<dyn EllipticBuilder>,
}

impl Decryption for Decryptor {
    fn execute(&self, cipher: &str) -> String {
        let data = {
            if !cipher.starts_with("04") {
                panic!("The cipher data is invalid.")
            }
            match hex::decode(cipher) {
                Ok(data) => data[1..].to_vec(),
                Err(_) => panic!("The cipher data must be composed of hex chars.")
            }
        };
        let (c1, c2, c3) = {
            let len = data.len();
            match self.mode {
                Mode::C1C3C2 => {
                    (data.clone()[..64].to_vec(), data.clone()[96..].to_vec(), data.clone()[64..96].to_vec())
                }
                Mode::C1C2C3 => {
                    (data.clone()[..64].to_vec(), data.clone()[64..len - 32].to_vec(), data.clone()[len - 32..].to_vec())
                }
            }
        };


        let (x2, y2) = {
            let (x1, y1) = (
                BigUint::from_bytes_be(&c1.clone()[..32]),
                BigUint::from_bytes_be(&c1.clone()[32..])
            );
            self.builder.scalar_multiply(x1, y1, self.key.value())
        };


        let plain = {
            let temp = [x2.to_bytes_be(), y2.to_bytes_be()].concat();
            let t = kdf(temp, c2.len());

            if is_all_zero(t.clone()) {
                panic!("The cipher data is invalid.")
            }

            let mut plain = vec![];
            for i in 0..c2.len() {
                plain.push(c2.clone()[i] ^ t.clone()[i]);
            }
            plain
        };

        let hash = {
            let temp = [x2.to_bytes_be(), plain.clone(), y2.to_bytes_be()].concat();
            sm3::hash(&temp).to_vec()
        };

        if hash != c3 {
            panic!("The cipher data hash validation failed.");
        }

        String::from_utf8_lossy(plain.as_slice()).to_string()
    }
}


/// 秘钥派生函数
#[inline(always)]
fn kdf(data: Vec<u8>, len: usize) -> Vec<u8> {
    let mut counter: usize = 0x00000001;
    let mut result: Vec<u8> = vec![];
    let k = data.len() + 31 / 32;
    for i in 0..k {
        let temp = [data.as_slice(), to_bytes(counter).as_slice()].concat();
        let hash = sm3::hash(&temp);

        if (i + 1) == k && len % 32 != 0 {
            result = [result, hash[..(len % 32)].to_vec()].concat();
        } else {
            result = [result, hash.to_vec()].concat();
        }
        counter += 1;
    }
    result
}

#[inline(always)]
fn is_all_zero(data: Vec<u8>) -> bool {
    let mut flag = true;
    for i in 0..data.len() {
        if data[i] != 0 {
            flag = false;
            break;
        }
    }
    flag
}


#[inline(always)]
fn to_bytes(x: usize) -> [u8; 4] {
    let mut buf: [u8; 4] = [0; 4];

    buf[0] = (x >> 24) as u8;
    buf[1] = (x >> 16) as u8;
    buf[2] = (x >> 8) as u8;
    buf[3] = x as u8;

    buf
}



