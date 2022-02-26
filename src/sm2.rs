use std::rc::Rc;
use crate::sm2::ecc::{Crypto, Decryption, Encryption, Signature};
use crate::sm2::key::{HexKey, KeyGenerator, KeyPair, PrivateKey, PublicKey};
use crate::sm2::p256::P256Elliptic;

mod key;
mod ecc;
mod p256;


pub fn generate_keypair() -> (String, String) {
    let p256 = P256Elliptic::init();
    let generator = KeyGenerator::init(Box::new(p256));
    let pair = generator.gen_key_pair();
    (pair.prk().encode(), pair.puk().encode())
}

pub fn encrypt(public_key: &str, plain: &str) -> String {
    let crypto = Crypto::default();
    crypto.encryptor(PublicKey::decode(public_key)).execute(plain)
}

pub fn decrypt(private_key: &str, cipher: &str) -> String {
    let crypto = Crypto::default();
    crypto.decryptor(PrivateKey::decode(private_key)).execute(cipher)
}

pub fn encrypt_c1c2c3(public_key: &str, plain: &str) -> String {
    let crypto = Crypto::c1c2c3(Rc::new(P256Elliptic::init()));
    crypto.encryptor(PublicKey::decode(public_key)).execute(plain)
}

pub fn decrypt_c1c2c3(private_key: &str, cipher: &str) -> String {
    let crypto = Crypto::c1c2c3(Rc::new(P256Elliptic::init()));
    crypto.decryptor(PrivateKey::decode(private_key)).execute(cipher)
}

pub fn sign(private_key: &str, public_key: &str, plain: &str) -> String {
    let crypto = Crypto::default();
    let keypair = KeyPair::new(PrivateKey::decode(private_key), PublicKey::decode(public_key));
    hex::encode(crypto.signer(keypair).sign(&plain).encode())
}

pub fn verify(public_key: &str, plain: &str, signature: &str) -> bool {
    let crypto = Crypto::default();
    let s = Signature::decode(hex::decode(signature).unwrap().as_slice());
    crypto.verifier(PublicKey::decode(public_key)).verify(plain, &s)
}