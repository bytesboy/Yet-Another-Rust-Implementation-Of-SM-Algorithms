use std::cmp::Ordering;
use std::mem;
use std::sync::Once;

use num_bigint::{BigUint, ToBigInt};
use num_integer::Integer;

use crate::sm2::core::{Elliptic, EllipticProvider};
use crate::sm2::p256::params::{BASE_TABLE, EC_A, EC_B, EC_GX, EC_GY, EC_N, EC_P, P256FACTOR, RI};
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
        let point = P256AffinePoint(
            PayloadHelper::transform(&x.to_bigint().unwrap()),
            PayloadHelper::transform(&y.to_bigint().unwrap()),
        );
        point.multiply(k).restore()
    }

    fn scalar_base_multiply(&self, k: BigUint) -> (BigUint, BigUint) {
        let elliptic = self.ec.clone();
        let base = P256BasePoint {
            point: P256AffinePoint(
                PayloadHelper::transform(&elliptic.gx.to_bigint().unwrap()),
                PayloadHelper::transform(&elliptic.gy.to_bigint().unwrap()),
            ),
            order: elliptic.n,
        };
        base.multiply(k).restore()
    }
}


/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
/// Affine coordinates: (X = x/z^2, Y = y/z^3)  Y^2 = X^3 + aX +b
#[derive(Clone, Debug)]
struct P256AffinePoint(Payload, Payload);

/// 基点
#[derive(Clone, Debug)]
struct P256BasePoint {
    point: P256AffinePoint,
    order: BigUint,
}

impl P256BasePoint {
    ///  Scalar is a little-endian number. Note that the value of scalar must be less than the order of the group.
    fn bytes_of_scalar(&self, scalar: BigUint) -> [u8; 32] {
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
        scalar_bytes
    }
}

/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
#[derive(Debug)]
struct P256JacobianPoint(Payload, Payload, Payload);


impl P256AffinePoint {
    fn restore(&self) -> (BigUint, BigUint) {
        let x = PayloadHelper::restore(&self.0).to_biguint().unwrap();
        let y = PayloadHelper::restore(&self.1).to_biguint().unwrap();
        (x, y)
    }

    /// get the entry of table by index.
    /// On entry: index < 16, table[0] must be zero.
    fn select(index: u32, table: Vec<u32>) -> Self {
        let mut table = table.clone();

        let mut point = P256AffinePoint(
            Payload::init(), Payload::init(),
        );

        for i in 1..16 {
            let mut mask: u32 = i ^ index;
            mask |= mask >> 2;
            mask |= mask >> 1;
            mask &= 1;

            if mask == 0 {
                mask = u32::MAX;    // 4294967295
            } else {
                mask -= 1;
            }

            for j in 0..9 {
                point.0.data()[j] |= table[0] & mask;
                if table.len() > 0 {
                    table = table[1..].to_vec();
                }
            }

            for j in 0..9 {
                point.1.data()[j] |= table[0] & mask;
                if table.len() > 0 {
                    table = table[1..].to_vec();
                }
            }
        }
        point
    }
}

impl P256JacobianPoint {
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

    /// add_affine sets {xOut,yOut,zOut} = {x1,y1,z1} + {x2,y2,1}.
    /// (i.e. the second point is affine.)
    ///
    /// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    ///
    /// Note that this function does not handle P+P, infinity+P nor P+infinity correctly.
    fn add_affine(&self, affine: &P256AffinePoint) -> Self {
        let (x1, y1, z1) = (&self.0, &self.1, &self.2);
        let (x2, y2) = (&affine.0, &affine.1);

        let z1z1 = z1.square();
        let temp = z1.add(&z1);
        let u2 = x2.multiply(&z1z1);
        let z1z1z1 = z1.multiply(&z1z1);
        let s2 = y2.multiply(&z1z1z1);
        let h = u2.subtract(&x1);

        let i = h.add(&h).square();
        let j = h.multiply(&i);

        let r = s2.subtract(y1);
        let r = r.add(&r);

        let v = x1.multiply(&i);

        let z_out = temp.multiply(&h);
        let rr = r.square();

        let x_out = rr.subtract(&j).subtract(&v).subtract(&v);
        let temp = v.subtract(&x_out);

        let y_out = temp.multiply(&r);
        let temp = y1.multiply(&j);
        let y_out = y_out.subtract(&temp).subtract(&temp);

        P256JacobianPoint(x_out, y_out, z_out)
    }

    /// sets out=source if mask = 0xffffffff in constant time.
    /// On entry: mask is either 0 or 0xffffffff.
    fn copy_from(&self, source: P256JacobianPoint, mask: u32) -> Self {
        let mut point = P256JacobianPoint(
            Payload::init(), Payload::init(), Payload::init(),
        );

        for i in 0..9 {
            point.0.data()[i] = self.0.data()[i] ^ (mask & (source.0.data()[i] ^ self.0.data()[i]));
            point.1.data()[i] = self.1.data()[i] ^ (mask & (source.1.data()[i] ^ self.1.data()[i]));
            point.2.data()[i] = self.2.data()[i] ^ (mask & (source.2.data()[i] ^ self.2.data()[i]));
        }
        point
    }

    /// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
    /// Affine coordinates: (X = x/z^2, Y = y/z^3)  Y^2 = X^3 + aX +b
    fn to_affine(&self) -> P256AffinePoint {
        let elliptic = P256Elliptic::init();
        let z = PayloadHelper::restore(&self.2);
        let p = elliptic.ec.p.to_bigint().unwrap();
        let zi = z.extended_gcd(&p).x.mod_floor(&p);

        let alpha = PayloadHelper::transform(&zi);
        let beta = alpha.square();
        let gama = alpha.multiply(&beta);

        let x = self.0.multiply(&beta);
        let y = self.1.multiply(&gama);

        P256AffinePoint(x, y)
    }
}


trait Multiplication {
    fn multiply(&self, scalar: BigUint) -> P256AffinePoint;
}

impl Multiplication for P256AffinePoint {
    fn multiply(&self, scalar: BigUint) -> P256AffinePoint {
        todo!()
    }
}

impl Multiplication for P256BasePoint {
    /// multiply sets P256Point = scalar*G where scalar is a little-endian number.
    fn multiply(&self, scalar: BigUint) -> P256AffinePoint {
        let scalar = self.bytes_of_scalar(scalar);

        let mut jacobian_point = P256JacobianPoint(
            Payload::init(), Payload::init(), Payload::init(),
        );
        let mut jacobian_temp_point = P256JacobianPoint(
            Payload::init(), Payload::init(), Payload::init(),
        );

        let mut n_is_infinity_mask = !(0 as u32);   // u32::MAX
        // The loop adds bits at positions 0, 64, 128 and 192, followed by positions 32, 96, 160
        // and 224 and does this 32 times.
        for i in 0..32 {
            if i != 0 {
                jacobian_point = jacobian_point.double();
            }
            let mut offset = 0;
            let mut j = 0;
            while j <= 32 {
                let bit0 = bit_of_scalar(scalar, 31 - i + j);
                let bit1 = bit_of_scalar(scalar, 95 - i + j);
                let bit2 = bit_of_scalar(scalar, 159 - i + j);
                let bit3 = bit_of_scalar(scalar, 223 - i + j);
                let idx = bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3);

                let affine_point = P256AffinePoint::select(
                    idx,
                    Vec::from(&BASE_TABLE[offset..]),
                );

                offset += 30 * 9;

                jacobian_temp_point = jacobian_point.add_affine(&affine_point);
                jacobian_point = jacobian_point.copy_from(
                    P256JacobianPoint(
                        affine_point.0.clone(),
                        affine_point.1.clone(),
                        Payload::new(P256FACTOR[1]),
                    ),
                    n_is_infinity_mask,
                );

                let p_is_finite_mask = mask(idx);
                let mask = p_is_finite_mask & !n_is_infinity_mask;

                jacobian_point = jacobian_point.copy_from(jacobian_temp_point, mask);

                // If p was not zero, then n is now non-zero.
                n_is_infinity_mask = n_is_infinity_mask & !p_is_finite_mask;

                j += 32;
            }
        }
        jacobian_point.to_affine()
    }
}


/// 0xffffffff for 0 < x <= 2^31  0xffffffff = 4294967295 = u32::MAX = 2^31 - 1
/// 0 for x == 0 or x > 2^31.
#[inline(always)]
pub(crate) fn mask(x: u32) -> u32 {
    x.wrapping_sub(1).wrapping_shr(31).wrapping_sub(1)
}

#[inline(always)]
fn bit_of_scalar(scalar: [u8; 32], bit: usize) -> u32 {
    (((scalar[bit >> 3]) >> (bit & 7)) & 1) as u32
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn double() {
        let p = P256JacobianPoint(
            Payload::new([
                142920515, 258221801, 612883394, 247790219, 102162616, 256181319, 368653124, 339147441, 485647861
            ]),
            Payload::new([
                131716495, 257805590, 847457731, 9891469, 365916039, 10897717, 75399777, 345048710, 61672909
            ]),
            Payload::new([
                91126934, 246575011, 35050116, 166561688, 126087236, 206595946, 25361097, 132288796, 249238939
            ]),
        );

        let point = p.double();

        let dx: [u32; 9] = [
            63255407, 227631960, 723093165, 65361332, 349345715, 60584340, 225318870, 397671582, 2985142
        ];
        let dy: [u32; 9] = [
            109858056, 93563162, 762162539, 50265907, 127330792, 104238630, 142585591, 352255388, 504506288
        ];
        let dz: [u32; 9] = [
            33808385, 18870127, 959285037, 176378705, 331289063, 266887158, 195778472, 241280794, 433045898
        ];

        assert_eq!(dx, point.0.data());
        assert_eq!(dy, point.1.data());
        assert_eq!(dz, point.2.data());
    }

    #[test]
    fn to_affine() {
        let jacobian = P256JacobianPoint(
            Payload::new([
                302587400, 224711462, 627912361, 12505049, 498636470, 226242352, 402285030, 277184676, 216966475
            ]),
            Payload::new([
                192016430, 212978101, 582317843, 172876572, 311643684, 126400666, 241514474, 362965479, 507691953
            ]),
            Payload::new([
                186636191, 229928314, 430146881, 262724875, 500465416, 219885119, 175182585, 128499041, 217581763
            ]),
        );
        let x: [u32; 9] = [
            194013013, 230698553, 317844872, 128801727, 111436768, 164685344, 76578606, 217356592, 311205467
        ];
        let y: [u32; 9] = [
            26049626, 112805900, 275795042, 259495837, 289529507, 146296588, 220416178, 146512122, 266185762
        ];

        let p = jacobian.to_affine();
        assert_eq!(p.0.data(), x);
        assert_eq!(p.1.data(), y);
    }


    #[test]
    fn add_affine() {
        let p1 = P256JacobianPoint(
            Payload::new([434464579, 232242225, 833663495, 95183971, 197589781, 65481707, 285356080, 397523777, 297319517]),
            Payload::new([105546064, 115648734, 616445926, 160673803, 382296094, 254935631, 24241561, 306433971, 112469103]),
            Payload::new([181993035, 232241130, 971204483, 180652253, 65532229, 175247468, 61056085, 229359646, 398806318]),
        );
        let p2 = P256AffinePoint(
            Payload::new([202984782, 49108071, 232741480, 255396639, 514738327, 218206935, 297234813, 116067631, 179908071]),
            Payload::new([5218908, 153082273, 421504040, 11374625, 412716736, 202538972, 20283405, 71924911, 112328172]),
        );

        let p3 = P256JacobianPoint(
            Payload::new([167460039, 227362747, 1005076632, 178921945, 76659602, 171371270, 426799015, 160435985, 428642590]),
            Payload::new([464015293, 22901587, 945207532, 41039408, 413094493, 244768035, 503070920, 229068862, 132259568]),
            Payload::new([404366665, 62541307, 262912748, 158805496, 464033083, 30021392, 180319644, 142373381, 27655256]),
        );

        let p = p1.add_affine(&p2);
        assert_eq!(p.0.data(), p3.0.data());
        assert_eq!(p.1.data(), p3.1.data());
        assert_eq!(p.2.data(), p3.2.data());
    }
}