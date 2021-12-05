use crate::sm4::core::Crypto;
use crate::sm4::{Cryptographer, xor};

/// OFB: Output FeedBack
///
/// 输出反馈模式
///
/// ### 推荐用CTR模式代替
///
/// 优点：
/// * 不需要填充
/// * 可事先进行加密、解密准备
/// * 加密、解密使用相同结构
/// * 对包含某些错误比特的密文进行解密时，只有明文中相应的比特会出错
///
/// 缺点：
/// * 不支持并行计算
/// * 主动攻击者反转密文分组中的某些比特时，明文分组中相应的比特也会被反转
pub struct CryptoMode {
    crypto: Crypto,
    iv: Vec<u8>,
}


impl CryptoMode {
    pub fn new(key: &[u8], iv: &[u8]) -> Self {
        crate::sm4::ofb::CryptoMode { crypto: Crypto::init(key), iv: iv.to_vec() }
    }
}

impl Cryptographer for CryptoMode {
    fn encrypt_bytes(&self, plain: &[u8]) -> Vec<u8> {
        let (quotients, remainder) = (plain.len() / 16, plain.len() % 16);

        let mut out: Vec<u8> = Vec::new();
        let mut buff = [0; 16];
        buff.copy_from_slice(&self.iv);

        for i in 0..quotients {
            let k = self.crypto.encrypt(&buff);
            let block = xor(&plain[i * 16..(i + 1) * 16], &k);
            block.iter().for_each(|e| out.push(*e));
            buff.copy_from_slice(&k);
        }
        let k = self.crypto.encrypt(&buff);
        // 最后一个分组，按位异或
        for i in 0..remainder {
            out.push(plain[quotients * 16 + i] ^ k[i]);
        }
        out
    }

    fn decrypt_bytes(&self, cipher: &[u8]) -> Vec<u8> {
        self.encrypt_bytes(cipher)
    }
}

#[cfg(test)]
mod tests {
    use crate::sm4::ofb::CryptoMode;
    use crate::sm4::Cryptographer;

    #[test]
    fn main() {
        let key = hex::decode("0123456789abcdeffedcba9876543210").unwrap();
        let iv = hex::decode("0123456789abcdeffedcba9876543210").unwrap();

        let plain = "Hello World, 哈哈，我的世界";

        let c = CryptoMode::new(&key, &iv);
        let cipher = c.encrypt(String::from(plain));
        let text = c.decrypt(cipher);

        assert_eq!(plain, text);
    }
}