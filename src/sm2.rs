use std::rc::Rc;

use crate::sm2::ecc::{Crypto, Decryption, Encryption};
use crate::sm2::key::{HexKey, KeyGenerator, PrivateKey, PublicKey};
use crate::sm2::p256::P256Elliptic;

mod key;
mod ecc;
mod p256;


pub fn generate_key() -> (String, String) {
    let p256 = P256Elliptic::init();
    let generator = KeyGenerator::init(Box::new(p256));
    let pair = generator.gen_key_pair();
    (pair.prk().encode(), pair.puk().encode())
}

pub fn encrypt(key: &str, plain: &str) -> String {
    let crypto = Crypto::default();
    crypto.encryptor(PublicKey::decode(key)).execute(plain)
}

pub fn decrypt(key: &str, cipher: &str) -> String {
    let crypto = Crypto::default();
    crypto.decryptor(PrivateKey::decode(key)).execute(cipher)
}