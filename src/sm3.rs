mod core;

/// 计算摘要信息
pub fn digest(data: &str) -> String {
    core::Crypto::init(data.as_bytes()).pad().block().iterate().hash()
}




