mod core;

pub fn digest(data: &str) -> String {
    core::Crypto::new(data.as_bytes()).pad().block().iterate().hash()
}




