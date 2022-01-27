use std::cmp::Ordering;
use std::mem;
use std::sync::Once;

use num_bigint::{BigUint, ToBigInt};
use num_integer::Integer;

use crate::sm2::core::{Elliptic, EllipticProvider};
use crate::sm2::p256::params::{EC_A, EC_B, EC_GX, EC_GY, EC_N, EC_P, RI};
use crate::sm2::p256::payload::{Payload, PayloadHelper};

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

/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
#[derive(Debug)]
struct P256JacobianPoint(Payload, Payload, Payload);

impl P256JacobianPoint {
    // fn new(x: PayLoad, y: PayLoad, z: PayLoad) -> Self {
    //     P256JacobianPoint(x, y, z)
    // }

    /// (x, y, z) => 2 * (x, y, z)
    /// [Formulas](https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l)
    fn double(&self) -> Self {
        let a = PayloadHelper::transform(&P256Elliptic::init().ec.a.to_bigint().unwrap());
        let (x, y, z) = (&self.0, &self.1, &self.2);

        let (alpha, beta) = (z.square(), y.square());
        // delta = 4xy^2
        let delta = x.multiply(&beta).scalar_multiply(4);
        // t1 = az^4
        let t1 = a.multiply(&alpha.square());
        // t2 = 8y^4
        let t2 = beta.square().scalar_multiply(8);
        // gama = 3x^2 + az^4
        let gama = x.square().scalar_multiply(3).add(&t1);
        // rx = (3x^2 + az^4)^2 - 8xy^2
        let rx = gama.square().subtract(&delta).subtract(&delta);
        let ry = delta.subtract(&rx).multiply(&gama).subtract(&t2);
        // rz = (y+z)^2 - z^2 - y^2
        let rz = y.add(&z).square().subtract(&alpha).subtract(&beta);

        P256JacobianPoint(rx, ry, rz)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn double() {
        let p = P256JacobianPoint(
            Payload::new([142920515, 258221801, 612883394, 247790219, 102162616, 256181319, 368653124, 339147441, 485647861]),
            Payload::new([131716495, 257805590, 847457731, 9891469, 365916039, 10897717, 75399777, 345048710, 61672909]),
            Payload::new([91126934, 246575011, 35050116, 166561688, 126087236, 206595946, 25361097, 132288796, 249238939]),
        );

        let point = p.double();

        let dx: [u32; 9] = [63255407, 227631960, 723093165, 65361332, 349345715, 60584340, 225318870, 397671582, 2985142];
        let dy: [u32; 9] = [109858056, 93563162, 762162539, 50265907, 127330792, 104238630, 142585591, 352255388, 504506288];
        let dz: [u32; 9] = [33808385, 18870127, 959285037, 176378705, 331289063, 266887158, 195778472, 241280794, 433045898];

        assert_eq!(dx, point.0.data());
        assert_eq!(dy, point.1.data());
        assert_eq!(dz, point.2.data());
    }
}