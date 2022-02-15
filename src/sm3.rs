mod core;

/// 计算摘要信息
pub fn digest(data: &str) -> String {
    hex::encode(hash(data.as_bytes()))
}

/// 计算Hash值
pub fn hash(data: vec<u8>) -> [u8; 32] {
    core::Crypto::new(data).hash()
}


