use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Mul, Sub};
use std::rc::Rc;

use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::sm2::key::{KeyPair, PrivateKey, PublicKey, to_32_bytes};
use crate::sm2::p256::P256Elliptic;
use crate::sm3;

const UID: [u8; 16] = [
    0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38,
];

pub trait EllipticBuilder {
    fn blueprint(&self) -> &Elliptic;

    /// 点加
    fn point_add(&self, x1: BigUint, y1: BigUint, x2: BigUint, y2: BigUint) -> (BigUint, BigUint);
    /// 标量乘法
    fn scalar_multiply(&self, x: BigUint, y: BigUint, scalar: BigUint) -> (BigUint, BigUint);
    /// 基点标量乘法
    fn scalar_base_multiply(&self, scalar: BigUint) -> (BigUint, BigUint);
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
    /// 随机数 【from, to】
    pub fn random(&self, from: BigUint, to: BigUint) -> BigUint {
        let temp = match from.clone().cmp(&to) {
            Ordering::Greater => from.clone().sub(&to),
            Ordering::Less => to.clone().sub(&from),
            Ordering::Equal => BigUint::zero()
        };

        if temp.is_zero() {
            return from.clone();
        }

        let k = {
            let k: Vec<u8> = (0..self.bits / 8 + 8).map(|_| { rand::random::<u8>() }).collect();
            BigUint::from_bytes_be(&k)
        };

        // temp = to - from + 1
        let temp = temp.add(BigUint::one());
        // k % temp  ∈ [0, temp - 1] = [0, to - from]  =>  k % temp + from ∈ [from, to]
        k.mod_floor(&temp).add(&from)
    }

    pub fn scalar_reduce(&self, scalar: BigUint) -> BigUint {
        // compare scalar and order, n = (scalar mod order) if scalar > order else scalar
        if let Ordering::Greater = scalar.cmp(&self.n) {
            scalar.mod_floor(&self.n)
        } else {
            scalar
        }
    }
}

#[derive(Debug, Copy, Clone)]
enum Mode {
    C1C2C3,
    C1C3C2,
}

pub struct Crypto {
    mode: Mode,
    builder: Rc<dyn EllipticBuilder>,
}

impl Crypto {
    pub fn default() -> Self {
        Self::c1c3c2(Rc::new(P256Elliptic::init()))
    }

    pub fn c1c2c3(builder: Rc<dyn EllipticBuilder>) -> Self {
        Crypto { mode: Mode::C1C2C3, builder }
    }

    pub fn c1c3c2(builder: Rc<dyn EllipticBuilder>) -> Self {
        Crypto { mode: Mode::C1C3C2, builder }
    }

    pub fn encryptor(&self, key: PublicKey) -> Encryptor {
        Encryptor { key, mode: self.mode, builder: self.builder.clone() }
    }

    pub fn decryptor(&self, key: PrivateKey) -> Decryptor {
        Decryptor { key, mode: self.mode, builder: self.builder.clone() }
    }

    pub fn signer(&self, keypair: KeyPair) -> Signer {
        let za = self.digest(keypair.puk().clone());
        Signer { hash: za, keypair, builder: self.builder.clone() }
    }

    pub fn verifier(&self, key: PublicKey) -> Verifier {
        let za = self.digest(key.clone());
        Verifier { hash: za, key, builder: self.builder.clone() }
    }

    /// ZA=H256(ENTLA ∥ IDA ∥ a ∥ b ∥ xG ∥ yG ∥xA ∥yA)
    fn digest(&self, puk: PublicKey) -> Vec<u8> {
        let ent = {
            if UID.len() >= 8192 {
                panic!("UID is too large.");
            }
            let r = UID.len() * 8;
            [((r >> 8) & 0xFF) as u8, (r & 0xFF) as u8].to_vec()
        };

        let id = UID.to_vec();
        let e = self.builder.blueprint();
        let (a, b) = (e.a.to_bytes_be(), e.a.to_bytes_be());
        let (gx, gy) = (e.gx.to_bytes_be(), e.gy.to_bytes_be());

        let (px, py) = {
            let key = puk.value();
            let (x, y) = (key.0.to_bytes_be(), key.1.to_bytes_be());
            (to_32_bytes(x).to_vec(), to_32_bytes(y).to_vec())
        };

        sm3::hash([ent, id, a, b, gx, gy, px, py].concat().as_slice()).to_vec()
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
    /// 加密
    fn execute(&self, plain: &str) -> String {
        let data = plain.as_bytes();
        let cipher = loop {
            let k = {
                let elliptic = self.builder.blueprint();
                let from = BigUint::one();
                elliptic.random(from.clone(), elliptic.n.clone().sub(&from.clone()))
            };

            // C1: [k]G
            let c1 = {
                let (x1, y1) = self.builder.scalar_base_multiply(k.clone());
                [vec![0x04], x1.to_bytes_be(), y1.to_bytes_be()].concat()
            };

            let (x2, y2) = {
                let key = self.key.value();
                let (x, y) = (key.0.clone(), key.1.clone());
                self.builder.scalar_multiply(x, y, k.clone())
            };

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
    /// 解密
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


#[derive(Debug, Clone)]
pub struct Signature {
    r: BigUint,
    s: BigUint,
}

impl Display for Signature {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f, "Signature {{ r: {:0>64}, s: {:0>64} }}",
            self.r.to_str_radix(16),
            self.s.to_str_radix(16)
        )
    }
}

impl Signature {
    pub(crate) fn new(r: BigUint, s: BigUint) -> Self {
        Signature { r, s }
    }

    /// Encodes the signature to DER-encoded ASN.1 data.
    pub(crate) fn encode(&self) -> Vec<u8> {
        let data = yasna::construct_der(|writer| {
            writer.write_sequence(|writer| {
                writer.next().write_biguint(&self.r);
                writer.next().write_biguint(&self.s);
            })
        });
        data
    }

    /// Decodes the DER-encoded ASN.1 data to Signature.
    pub(crate) fn decode(signature: &[u8]) -> Self {
        let (r, s) = yasna::parse_der(signature, |reader| {
            reader.read_sequence(|reader| {
                let r = reader.next().read_biguint()?;
                let s = reader.next().read_biguint()?;
                Ok((r, s))
            })
        }).unwrap();

        Signature::new(r, s)
    }
}

pub struct Signer {
    hash: Vec<u8>,
    keypair: KeyPair,
    builder: Rc<dyn EllipticBuilder>,
}

impl Signer {
    /// 签名
    pub(crate) fn sign(&self, plain: &str) -> Signature {
        let m = [self.hash.clone(), plain.as_bytes().to_vec()].concat();
        let e = sm3::hash(m.as_slice());
        let elliptic = self.builder.blueprint();

        let key = self.keypair.prk();

        let (r, s) = loop {
            let k = {
                let from = BigUint::one();
                elliptic.random(from.clone(), elliptic.n.clone().sub(&from.clone()))
            };

            let r = {
                let (x, _) = self.builder.scalar_base_multiply(k.clone());
                BigUint::from_bytes_be(&e).add(&x).mod_floor(&elliptic.n)
            };

            if r == BigUint::zero() || r.clone().add(k.clone()) == elliptic.n {
                continue;
            }

            let s = {
                let n = elliptic.n.to_bigint().unwrap();
                let d = key.value().to_bigint().unwrap();
                let temp = d.clone().mul(&r.to_bigint().unwrap());
                // a = k - rd
                let a = k.to_bigint().unwrap().sub(&temp);
                let temp = d.clone().add(BigInt::one());
                // 1 / (1+d)
                let b = temp.extended_gcd(&n).x.mod_floor(&n);
                a.mul(b).mod_floor(&n).to_biguint().unwrap()
            };

            if s == BigUint::zero() {
                continue;
            }

            break (r, s);
        };

        Signature::new(r, s)
    }
}


pub struct Verifier {
    hash: Vec<u8>,
    key: PublicKey,
    builder: Rc<dyn EllipticBuilder>,
}

impl Verifier {
    /// 验签
    pub(crate) fn verify(&self, plain: &str, signature: &Signature) -> bool {
        let elliptic = self.builder.blueprint();
        let n1 = elliptic.n.clone().sub(BigUint::one());
        let (r, s) = (signature.r.clone(), signature.s.clone());

        if r < BigUint::one() || r > n1.clone() {
            return false;
        }

        if s < BigUint::one() || s > n1.clone() {
            return false;
        }

        let e = {
            let m = [self.hash.clone(), plain.as_bytes().to_vec()].concat();
            let h = sm3::hash(m.as_slice());
            BigUint::from_bytes_be(h.as_slice())
        };

        let t = r.clone().add(&s).mod_floor(&elliptic.n);

        if BigUint::zero().eq(&t) {
            return false;
        }

        let x = {
            let key = self.key.value();
            let p1 = self.builder.scalar_base_multiply(s.clone());
            let p2 = self.builder.scalar_multiply(key.0, key.1, t);
            let p3 = self.builder.point_add(p1.0, p1.1, p2.0, p2.1);
            p3.0
        };

        let rn = e.add(x).mod_floor(&elliptic.n);

        return if rn == r {
            true
        } else {
            false
        };
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo() {
        println!("BigUint::one() = {:?}", BigUint::one());
    }
}