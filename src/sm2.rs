use std::rc::Rc;
use crate::sm2::ecc::{Crypto, Decryption, Encryption, Mode};
use crate::sm2::key::{HexKey, KeyGenerator, PrivateKey, PublicKey};
use crate::sm2::p256::P256Elliptic;
use crate::sm4::CryptoFactory;

mod key;
mod ecc;
mod p256;


pub fn generate_key() -> (String, String) {
    let p256 = P256Elliptic::init();
    let generator = KeyGenerator::init(Box::new(p256));
    let pair = generator.gen_key_pair();
    let (prk, puk) = (pair.private_key(), pair.public_key());
    (prk.encode(), puk.encode())
}

pub fn encrypt(key: &str, plain: &str) -> String {
    let puk = PublicKey::decode(key);
    let builder = Rc::new(P256Elliptic::init());
    let crypto = Crypto::init(Mode::C1C3C2, builder);
    crypto.encryptor(puk).execute(plain)
}

pub fn decrypt(key: &str, cipher: &str) -> String {
    let prk = PrivateKey::decode(key);
    let builder = Rc::new(P256Elliptic::init());
    let crypto = Crypto::init(Mode::C1C3C2, builder);
    crypto.decryptor(prk).execute(cipher)
}