mod core;

/// 计算摘要信息：Hash值编码为Hex字符串
pub fn digest(data: &str) -> String {
    hex::encode(hash(data.as_bytes()))
}

/// 计算Hash值
pub fn hash(data: &[u8]) -> [u8; 32] {
    core::Crypto::new(data).hash()
}


