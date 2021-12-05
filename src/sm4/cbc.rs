use crate::sm4::core::Crypto;
use crate::sm4::{Cryptographer, xor};

/// CBC: Cipher Block Chaining
///
/// 密码块链模式
///
/// ### 推荐使用
///
/// 优点：
/// * 明文的重复排列不会反映在密文中
/// * 支持并行计算（仅解密）
/// * 能够解密任意密文分组
///
/// 缺点：
/// * 对包含某些错误比特的密文进行解密时，第一个分组的全部比特以及最后一个分组的相应比特会出错。
/// * 加密不支持并行计算
///
/// 注意：加解密使用的初始化向量值必须相同
pub struct CryptoMode {
    crypto: Crypto,
    iv: Vec<u8>,
}

impl CryptoMode {
    pub fn new(key: &[u8], iv: &[u8]) -> Self {
        crate::sm4::cbc::CryptoMode { crypto: Crypto::init(key), iv: iv.to_vec() }
    }
}

impl Cryptographer for CryptoMode {
    fn encrypt_bytes(&self, plain: &[u8]) -> Vec<u8> {
        // 计算分组，每个分组应该是满16字节。最后一个分组要么是明文+填充总共满足16字节，要么是全填充16字节
        // 填充数据原则：(16-remainder)个(16-remainder)
        let (quotients, remainder) = (plain.len() / 16, plain.len() % 16);

        let mut out: Vec<u8> = Vec::new();
        let mut buf = [0; 16];
        buf.copy_from_slice(&self.iv);

        for i in 0..quotients {
            let block = xor(&buf, &plain[i * 16..(i + 1) * 16]);
            let cipher = self.crypto.encrypt(&block);

            out.extend_from_slice(&cipher);
            buf = cipher;
        }

        if remainder != 0 {
            // 如果数据长度除以16有余数，那就补充(16-余数)个(16-余数)
            let mut last = [(16 - remainder) as u8; 16];
            last[..remainder].copy_from_slice(&plain[quotients * 16..]);
            let block = xor(&buf, &last);
            let cipher = self.crypto.encrypt(&block);
            out.extend_from_slice(&cipher);
        } else {
            // 如果数据长度正好是16的倍数，那就补充16个字节,补充数据为0x10=16
            let block = xor(&buf, &[0x10; 16]);
            let cipher = self.crypto.encrypt(&block);
            out.extend_from_slice(&cipher);
        }
        out
    }

    fn decrypt_bytes(&self, cipher: &[u8]) -> Vec<u8> {
        let (quotients, remainder) = (cipher.len() / 16, cipher.len() % 16);
        if remainder != 0 {
            panic!("The cipher‘s length must be a multiple of 16 bytes.");
        }

        let mut out: Vec<u8> = Vec::new();
        let mut buf = [0; 16];
        buf.copy_from_slice(&self.iv);

        for i in 0..quotients {
            let block = self.crypto.decrypt(&cipher[i * 16..(i + 1) * 16]);
            let plain = xor(&buf, &block);
            plain.iter().for_each(|e| out.push(*e));
            buf.copy_from_slice(&cipher[i * 16..(i + 1) * 16])
        }

        let last_byte = out[cipher.len() - 1];
        // assert!(last_byte > 0 && last_byte <= 0x10);
        out.resize(cipher.len() - last_byte as usize, 0);
        out
    }
}


#[cfg(test)]
mod tests {
    use crate::sm4::cbc::CryptoMode;
    use crate::sm4::Cryptographer;

    #[test]
    fn main() {
        let key = hex::decode("0123456789abcdeffedcba9876543210").unwrap();
        let iv = hex::decode("0123456789abcdeffedcba9876543210").unwrap();

        let plain = "Hello World, 哈罗，世界";

        let c = CryptoMode::new(&key, &iv);
        let cipher = c.encrypt(String::from(plain));
        let text = c.decrypt(cipher);

        assert_eq!(plain, text);
    }
}