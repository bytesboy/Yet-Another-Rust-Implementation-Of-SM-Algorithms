pub mod sm3;


#[cfg(test)]
mod tests {
    use crate::sm3;

    #[test]
    fn sm3_hash() {
        let hash = sm3::digest("abc");
        assert_eq!(hash, "66c7f0f462eeedd9d1f2d46bdc10e4e24167c4875cf2f7a2297da02b8f4ba8e0");
    }
}