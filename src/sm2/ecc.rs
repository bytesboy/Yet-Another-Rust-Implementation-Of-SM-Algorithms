use std::borrow::Borrow;
use std::cmp::Ordering;
use std::rc::Rc;
use std::thread::Builder;

use num_bigint::BigUint;
use num_integer::Integer;

use crate::sm2::key::{PrivateKey, PublicKey};

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

#[derive(PartialEq)]
pub enum Mode {
    C1C2C3,
    C1C3C2,
}


pub struct CryptoFactory {
    mode: Mode,
    builder: Rc<dyn EllipticBuilder>,
}

impl CryptoFactory {
    pub fn init(builder: Rc<dyn EllipticBuilder>) -> Self {
        CryptoFactory { mode: Mode::C1C2C3, builder }
    }

    pub fn mode(&mut self, mode: Mode) -> &mut Self {
        if self.mode != mode {
            self.mode = mode;
        }
        self
    }

    pub fn encryptor(&self, key: PublicKey) -> Box<dyn Encryption> {
        match self.mode {
            Mode::C1C2C3 => {
                Box::new(c1c2c3::Encryptor::new(key, self.builder.clone()))
            }
            Mode::C1C3C2 => {
                Box::new(c1c3c2::Encryptor::new(key, self.builder.clone()))
            }
        }
    }

    pub fn decryptor(&self, key: PrivateKey) -> Box<dyn Decryption> {
        match self.mode {
            Mode::C1C2C3 => {
                Box::new(c1c2c3::Decryptor::new(key, self.builder.clone()))
            }
            Mode::C1C3C2 => {
                Box::new(c1c3c2::Decryptor::new(key, self.builder.clone()))
            }
        }
    }
}


pub trait Encryption {
    fn execute(&self, plain: &str) -> String;
}

pub trait Decryption {
    fn execute(&self, cipher: &str) -> String;
}


mod c1c2c3 {
    use std::rc::Rc;
    use crate::sm2::ecc::{Decryption, EllipticBuilder, Encryption};
    use crate::sm2::key::{PrivateKey, PublicKey};

    pub struct Encryptor {
        key: PublicKey,
        builder: Rc<dyn EllipticBuilder>,
    }

    impl Encryptor {
        pub(crate) fn new(key: PublicKey, builder: Rc<dyn EllipticBuilder>) -> Self {
            Encryptor { key, builder }
        }
    }

    impl Encryption for Encryptor {
        fn execute(&self, plain: &str) -> String {
            println!("plain = {:?}", plain);
            String::from(plain)
        }
    }


    pub struct Decryptor {
        key: PrivateKey,
        builder: Rc<dyn EllipticBuilder>,
    }

    impl Decryptor {
        pub(crate) fn new(key: PrivateKey, builder: Rc<dyn EllipticBuilder>) -> Self {
            Decryptor { key, builder }
        }
    }

    impl Decryption for Decryptor {
        fn execute(&self, cipher: &str) -> String {
            println!("cipher = {:?}", cipher);
            String::from(cipher)
        }
    }
}


mod c1c3c2 {
    use std::rc::Rc;
    use crate::sm2::ecc::{Decryption, EllipticBuilder, Encryption};
    use crate::sm2::key::{PrivateKey, PublicKey};

    pub struct Encryptor {
        key: PublicKey,
        builder: Rc<dyn EllipticBuilder>,
    }

    impl Encryptor {
        pub(crate) fn new(key: PublicKey, builder: Rc<dyn EllipticBuilder>) -> Self {
            Encryptor { key, builder }
        }
    }

    impl Encryption for Encryptor {
        fn execute(&self, plain: &str) -> String {
            todo!()
        }
    }

    pub struct Decryptor {
        key: PrivateKey,
        builder: Rc<dyn EllipticBuilder>,
    }

    impl Decryptor {
        pub(crate) fn new(key: PrivateKey, builder: Rc<dyn EllipticBuilder>) -> Self {
            Decryptor { key, builder }
        }
    }

    impl Decryption for Decryptor {
        fn execute(&self, cipher: &str) -> String {
            todo!()
        }
    }
}


