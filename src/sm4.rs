mod core;
mod ecb;
mod cbc;
mod cfb;
mod ofb;
mod ctr;


/// 随机生成秘钥，返回由16进制字符组成的长度为32的字符串
pub fn generate_key() -> String {
    let key = uuid::Uuid::new_v4();
    key.to_simple().to_string()
}

/// 随机初始化向量，返回由16进制字符组成的长度为32的字符串，用于CBC、CFB、OFB、CTR分组模式
pub fn generate_iv() -> String {
    generate_key()
}

pub fn encrypt_ecb(key: String, plain: String) -> String {
    let mode = Mode::ECB { key };
    let crypto = CryptoFactory::new(mode);
    crypto.encrypt(plain)
}

pub fn decrypt_ecb(key: String, cipher: String) -> String {
    let mode = Mode::ECB { key };
    let crypto = CryptoFactory::new(mode);
    crypto.decrypt(cipher)
}

pub fn encrypt_cbc(key: String, iv: String, plain: String) -> String {
    let mode = Mode::CBC { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.encrypt(plain)
}

pub fn decrypt_cbc(key: String, iv: String, cipher: String) -> String {
    let mode = Mode::CBC { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.decrypt(cipher)
}

pub fn encrypt_cfb(key: String, iv: String, plain: String) -> String {
    let mode = Mode::CFB { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.encrypt(plain)
}

pub fn decrypt_cfb(key: String, iv: String, cipher: String) -> String {
    let mode = Mode::CFB { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.decrypt(cipher)
}

pub fn encrypt_ofb(key: String, iv: String, plain: String) -> String {
    let mode = Mode::OFB { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.encrypt(plain)
}

pub fn decrypt_ofb(key: String, iv: String, cipher: String) -> String {
    let mode = Mode::OFB { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.decrypt(cipher)
}

pub fn encrypt_ctr(key: String, iv: String, plain: String) -> String {
    let mode = Mode::CTR { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.encrypt(plain)
}

pub fn decrypt_ctr(key: String, iv: String, cipher: String) -> String {
    let mode = Mode::CTR { key, iv };
    let crypto = CryptoFactory::new(mode);
    crypto.decrypt(cipher)
}

pub enum Mode {
    ECB { key: String },
    CBC { key: String, iv: String },
    CFB { key: String, iv: String },
    OFB { key: String, iv: String },
    CTR { key: String, iv: String },
}

pub trait Cryptographer {
    fn encrypt_bytes(&self, plain: &[u8]) -> Vec<u8>;

    fn decrypt_bytes(&self, cipher: &[u8]) -> Vec<u8>;

    fn encrypt(&self, data: String) -> String {
        let cipher = self.encrypt_bytes(data.as_bytes());
        hex::encode(cipher)
    }

    fn decrypt(&self, data: String) -> String {
        let plain = self.decrypt_bytes(&hex::decode(data).unwrap());
        String::from_utf8_lossy(plain.as_ref()).to_string()
    }
}

pub struct CryptoFactory;


impl CryptoFactory {
    pub fn new(mode: Mode) -> Box<dyn Cryptographer> {
        match mode {
            Mode::ECB { key } => {
                Box::new(ecb::CryptoMode::new(&hex_decode_of_key(&key)))
            }
            Mode::CBC { key, iv } => {
                Box::new(cbc::CryptoMode::new(&hex_decode_of_key(&key), &hex_decode_of_iv(&iv)))
            }
            Mode::CFB { key, iv } => {
                Box::new(cfb::CryptoMode::new(&hex_decode_of_key(&key), &hex_decode_of_iv(&iv)))
            }
            Mode::OFB { key, iv } => {
                Box::new(ofb::CryptoMode::new(&hex_decode_of_key(&key), &hex_decode_of_iv(&iv)))
            }
            Mode::CTR { key, iv } => {
                Box::new(ctr::CryptoMode::new(&hex_decode_of_key(&key), &hex_decode_of_iv(&iv)))
            }
        }
    }
}


fn xor(a: &[u8], b: &[u8]) -> [u8; 16] {
    let mut out: [u8; 16] = [0; 16];
    for i in 0..16 {
        out[i] = a[i] ^ b[i];
    }
    out
}

fn hex_decode_of_key(key: &str) -> Vec<u8> {
    match hex::decode(key) {
        Ok(data) => data,
        Err(_) => panic!("The Key must be composed of hex chars with a length of 32.")
    }
}

fn hex_decode_of_iv(iv: &str) -> Vec<u8> {
    match hex::decode(iv) {
        Ok(data) => data,
        Err(_) => panic!("The IV must be composed of hex chars with a length of 32.")
    }
}



