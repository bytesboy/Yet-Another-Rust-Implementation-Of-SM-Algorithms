mod core;

pub fn digest(data: &str) -> String {
    core::Cryptographer::new(data.as_bytes()).pad().block().iterate().hash()
}




