use std::mem;
use std::sync::Once;

use num_bigint::{BigUint, ToBigInt};

use crate::sm2::ecc::{Elliptic, EllipticProvider};
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

impl EllipticProvider for P256Elliptic {
    fn blueprint(&self) -> &Elliptic {
        &self.ec
    }

    fn scalar_multiply(&self, x: BigUint, y: BigUint, scalar: BigUint) -> (BigUint, BigUint) {
        let point = P256AffinePoint::new(
            PayloadHelper::transform(&x.to_bigint().unwrap()),
            PayloadHelper::transform(&y.to_bigint().unwrap()),
        );
        let scalar = self.scalar_to_bytes(scalar);
        point.multiply(scalar).restore()
    }

    fn scalar_base_multiply(&self, scalar: BigUint) -> (BigUint, BigUint) {
        let elliptic = self.ec.clone();
        let base = P256BasePoint::new(
            P256AffinePoint::new(
                PayloadHelper::transform(&elliptic.gx.to_bigint().unwrap()),
                PayloadHelper::transform(&elliptic.gy.to_bigint().unwrap()),
            ),
            elliptic.n,
        );

        let scalar = self.scalar_to_bytes(scalar);
        base.multiply(scalar).restore()
    }
}


/// 0xffffffff for 0 < x <= 2^31  0xffffffff = 4294967295 = u32::MAX = 2^31 - 1
/// 0 for x == 0 or x > 2^31.
#[inline(always)]
fn mask(x: u32) -> u32 {
    x.wrapping_sub(1).wrapping_shr(31).wrapping_sub(1)
}
