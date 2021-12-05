use crate::sm4::core::Crypto;
use crate::sm4::{Cryptographer, xor};

/// CFB: Cipher FeedBack
///
/// 密文反馈模式
///
/// ### 现已不使用，推荐用CTR模式代替
///
/// 优点：
/// * 不需要填充
/// * 支持并行计算（仅解密）
/// * 能够解密任意密文分组
///
/// 缺点：
/// * 对包含某些错误比特的密文进行解密时，第一个分组的全部比特以及最后一个分组的相应比特会出错。
/// * 加密不支持并行计算
/// * 不能抵御重放攻击
pub struct CryptoMode {
    crypto: Crypto,
    iv: Vec<u8>,
}

impl CryptoMode {
    pub fn new(key: &[u8], iv: &[u8]) -> Self {
        crate::sm4::cfb::CryptoMode { crypto: Crypto::init(key), iv: iv.to_vec() }
    }
}

impl Cryptographer for CryptoMode {
    fn encrypt_bytes(&self, plain: &[u8]) -> Vec<u8> {
        let (quotients, remainder) = (plain.len() / 16, plain.len() % 16);

        let mut out: Vec<u8> = Vec::new();
        let mut buf = [0; 16];
        buf.copy_from_slice(&self.iv);

        for i in 0..quotients {
            let c = self.crypto.encrypt(&buf);
            let block = xor(&c, &plain[i * 16..(i + 1) * 16]);
            block.iter().for_each(|e| out.push(*e));
            buf.copy_from_slice(&block);
        }
        let c = self.crypto.encrypt(&buf);
        // 最后一个分组，按位异或
        for i in 0..remainder {
            out.push(plain[quotients * 16 + i] ^ c[i])
        }
        out
    }


    // 解密和加密的结构是一样的
    fn decrypt_bytes(&self, cipher: &[u8]) -> Vec<u8> {
        let (quotients, remainder) = (cipher.len() / 16, cipher.len() % 16);

        let mut out: Vec<u8> = Vec::new();
        let mut buf = [0; 16];
        buf.copy_from_slice(&self.iv);

        for i in 0..quotients {
            // 注意：是加密不是解密
            let c = self.crypto.encrypt(&buf);
            let plain = xor(&c, &cipher[i * 16..(i + 1) * 16]);
            plain.iter().for_each(|e| out.push(*e));
            buf.copy_from_slice(&cipher[i * 16..(i + 1) * 16])
        }

        let c = self.crypto.encrypt(&buf);
        for i in 0..remainder {
            out.push(cipher[quotients * 16 + i] ^ c[i])
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use crate::sm4::cfb::CryptoMode;
    use crate::sm4::Cryptographer;

    #[test]
    fn main() {
        let key = hex::decode("0123456789abcdeffedcba9876543210").unwrap();
        let iv = hex::decode("0123456789abcdeffedcba9876543210").unwrap();

        let plain = "Hello World, 哈罗，魔兽世界";

        let c = CryptoMode::new(&key, &iv);
        let cipher = c.encrypt(String::from(plain));
        let text = c.decrypt(cipher);

        assert_eq!(plain, text);
    }
}