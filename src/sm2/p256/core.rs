use std::cmp::Ordering;
use std::mem;
use std::sync::Once;
use num_bigint::BigUint;
use num_integer::Integer;
use crate::sm2::core::{Elliptic, EllipticProvider};
use crate::sm2::p256::params::{EC_A, EC_B, EC_GX, EC_GY, EC_N, EC_P, RI};

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

    fn scalar_multiply(&self, x: BigUint, y: BigUint, k: BigUint) -> (BigUint, BigUint) {
        let point = P256Point(x, y);
        let p = point.multiply(k);
        (p.0, p.1)
    }

    fn scalar_base_multiply(&self, k: BigUint) -> (BigUint, BigUint) {
        let elliptic = self.ec.clone();
        let base = P256BasePoint {
            point: P256Point(elliptic.gx.clone(), elliptic.gy.clone()),
            order: elliptic.n,
        };
        let p = base.multiply(k);
        (p.0, p.1)
    }
}

trait Multiplication {
    fn multiply(&self, scalar: BigUint) -> P256Point;
}

/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
/// Affine coordinates: (X = x/z^2, Y = y/z^3)  Y^2 = X^3 + aX +b
#[derive(Clone, Debug)]
struct P256Point(BigUint, BigUint);

/// 基点
#[derive(Clone, Debug)]
struct P256BasePoint {
    point: P256Point,
    order: BigUint,
}

impl Multiplication for P256Point {
    fn multiply(&self, scalar: BigUint) -> P256Point {
        todo!()
    }
}

impl Multiplication for P256BasePoint {
    fn multiply(&self, scalar: BigUint) -> P256Point {
        let scalar = {
            // compare scalar and order, n = (scalar mod order) if scalar > order else scalar
            if let Ordering::Greater = scalar.cmp(&self.order) {
                scalar.mod_floor(&self.order)
            } else {
                scalar
            }
        };

        let mut scalar_bytes = [0u8; 32];
        for (i, v) in scalar.to_bytes_le().iter().enumerate() {
            scalar_bytes[i] = *v;
        }


        self.point.clone()
    }
}