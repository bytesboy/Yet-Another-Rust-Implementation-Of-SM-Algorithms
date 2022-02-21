use std::mem;
use std::sync::Once;

use num_bigint::{BigUint, ToBigInt};

use crate::sm2::ecc::{Elliptic, EllipticBuilder};
use crate::sm2::p256::params::{EC_A, EC_B, EC_GX, EC_GY, EC_N, EC_P, RI};
use crate::sm2::p256::payload::PayloadHelper;
use crate::sm2::p256::point::{Multiplication, P256AffinePoint, P256BasePoint};

mod point;
mod payload;
mod params;

#[derive(Clone, Debug)]
pub struct P256Elliptic {
    pub ec: Elliptic,
    pub ri: BigUint,
}

impl P256Elliptic {
    pub fn init() -> Self {
        static mut ELLIPTIC: *const P256Elliptic = std::ptr::null::<P256Elliptic>();
        static INITIALIZER: Once = Once::new();
        unsafe {
            INITIALIZER.call_once(|| {
                let p256 = P256Elliptic {
                    ec: Elliptic {
                        p: BigUint::from_bytes_be(&EC_P),
                        a: BigUint::from_bytes_be(&EC_A),
                        b: BigUint::from_bytes_be(&EC_B),
                        gx: BigUint::from_bytes_be(&EC_GX),
                        gy: BigUint::from_bytes_be(&EC_GY),
                        n: BigUint::from_bytes_be(&EC_N),
                        bits: 256,
                    },
                    ri: BigUint::from_bytes_be(&RI),
                };
                ELLIPTIC = mem::transmute(Box::new(p256));
            });
            (*ELLIPTIC).clone()
        }
    }
}

impl EllipticBuilder for P256Elliptic {
    fn blueprint(&self) -> &Elliptic {
        &self.ec
    }

    fn scalar_multiply(&self, x: BigUint, y: BigUint, scalar: BigUint) -> (BigUint, BigUint) {
        let elliptic = self.blueprint();
        let point = P256AffinePoint::new(
            PayloadHelper::transform(&x.to_bigint().unwrap()),
            PayloadHelper::transform(&y.to_bigint().unwrap()),
        );
        point.multiply(elliptic.scalar_reduce(scalar)).restore()
    }

    fn scalar_base_multiply(&self, scalar: BigUint) -> (BigUint, BigUint) {
        let elliptic = self.blueprint();
        let base = P256BasePoint::new(
            P256AffinePoint::new(
                PayloadHelper::transform(&elliptic.gx.to_bigint().unwrap()),
                PayloadHelper::transform(&elliptic.gy.to_bigint().unwrap()),
            ),
            elliptic.n.clone(),
        );
        base.multiply(elliptic.scalar_reduce(scalar)).restore()
    }
}


/// 0xffffffff for 0 < x <= 2^31  0xffffffff = 4294967295 = u32::MAX = 2^31 - 1
/// 0 for x == 0 or x > 2^31.
#[inline(always)]
fn mask(x: u32) -> u32 {
    x.wrapping_sub(1).wrapping_shr(31).wrapping_sub(1)
}


#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use crate::sm2::ecc::{Crypto, Decryption, Encryption, Signature};
    use crate::sm2::key::{HexKey, KeyPair, PrivateKey, PublicKey};

    use super::*;

    #[test]
    fn main() {
        let elliptic = P256Elliptic::init();

        let prk = "6aea1ccf610488aaa7fddba3dd6d76d3bdfd50f957d847be3d453defb695f28e";
        let puk = "04a8af64e38eea41c254df769b5b41fbaa2d77b226b301a2636d463c52b46c777230ad1714e686dd641b9e04596530b38f6a64215b0ed3b081f8641724c5443a6e";
        let private_key = PrivateKey::decode(prk);
        let public_key = PublicKey::decode(puk);

        let crypto = Crypto::c1c2c3(Rc::new(elliptic.clone()));
        let encryptor = crypto.encryptor(public_key.clone());
        let decryptor = crypto.decryptor(private_key.clone());
        let text = "兽人永不为奴，我们终将成王。——加尔鲁什·地狱咆哮";
        let cipher = encryptor.execute(text);
        let plain = decryptor.execute(&cipher);
        assert_eq!(plain, text);

        let crypto = Crypto::c1c3c2(Rc::new(elliptic.clone()));
        let encryptor = crypto.encryptor(public_key.clone());
        let decryptor = crypto.decryptor(private_key.clone());
        let text = "圣光会抛弃你的，英雄，就像抛弃我那样。——巫妖王";
        let cipher = encryptor.execute(text);
        let plain = decryptor.execute(&cipher);
        assert_eq!(plain, text);
    }


    #[test]
    fn signature() {
        let elliptic = P256Elliptic::init();

        let prk = "6aea1ccf610488aaa7fddba3dd6d76d3bdfd50f957d847be3d453defb695f28e";
        let puk = "04a8af64e38eea41c254df769b5b41fbaa2d77b226b301a2636d463c52b46c777230ad1714e686dd641b9e04596530b38f6a64215b0ed3b081f8641724c5443a6e";
        let prk = PrivateKey::decode(prk);
        let puk = PublicKey::decode(puk);
        let keypair = KeyPair::new(prk.clone(), puk.clone());

        let crypto = Crypto::c1c3c2(Rc::new(elliptic.clone()));

        let signer = crypto.signer(keypair);
        let s = signer.sign("圣光会抛弃你的，英雄，就像抛弃我那样。——巫妖王");

        let ans1 = hex::encode(s.encode());
        println!("ans1 = {:?}", ans1);
        let signature = Signature::decode(hex::decode(ans1).unwrap().as_slice());
        println!("sign = {:?}", signature);
    }
}