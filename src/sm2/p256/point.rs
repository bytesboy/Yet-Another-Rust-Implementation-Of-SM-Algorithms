use std::cmp::Ordering;
use std::ops::{BitAnd, Neg, Shr};

use num_bigint::{BigInt, BigUint, Sign, ToBigInt};
use num_integer::Integer;
use num_traits::{One, ToPrimitive};

use crate::sm2::p256::{mask, P256Elliptic};
use crate::sm2::p256::params::{BASE_TABLE, P256FACTOR};
use crate::sm2::p256::payload::{Payload, PayloadHelper};

pub(crate) trait Multiplication {
    fn multiply(&self, scalar: BigUint) -> P256AffinePoint;
}

/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
/// Affine coordinates: (X = x/z^2, Y = y/z^3)  Y^2 = X^3 + aX +b
#[derive(Clone, Debug)]
pub(crate) struct P256AffinePoint(Payload, Payload);

impl P256AffinePoint {
    pub(crate) fn new(x: Payload, y: Payload) -> Self {
        P256AffinePoint(x, y)
    }

    pub(crate) fn restore(&self) -> (BigUint, BigUint) {
        let x = PayloadHelper::restore(&self.0).to_biguint().unwrap();
        let y = PayloadHelper::restore(&self.1).to_biguint().unwrap();
        (x, y)
    }

    /// get the entry of table by index.
    /// On entry: index < 16, table[0] must be zero.
    fn select(index: u32, table: Vec<u32>) -> Self {
        let (mut x, mut y) = (Payload::init().data(), Payload::init().data());
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

            let offset = ((i - 1) * 18) as usize;
            for j in 0..9 {
                x[j] |= table[offset + j] & mask;
                y[j] |= table[offset + 9 + j] & mask;
            }
        }
        P256AffinePoint(Payload::new(x), Payload::new(y))
    }

    pub(crate) fn to_jacobian(&self) -> P256JacobianPoint {
        P256JacobianPoint(
            self.0.clone(),
            self.1.clone(),
            PayloadHelper::transform(&BigInt::one()),
        )
    }
}


impl Multiplication for P256AffinePoint {
    fn multiply(&self, scalar: BigUint) -> P256AffinePoint {
        let points = {
            let mut precomp: [[[u32; 9]; 3]; 16] = [[[0; 9]; 3]; 16];

            precomp[1][0] = self.0.data();
            precomp[1][1] = self.1.data();
            precomp[1][2] = P256FACTOR[1];

            let mut i = 2;
            while i < 8 {
                let p = P256JacobianPoint(
                    Payload::new(precomp[i / 2][0]),
                    Payload::new(precomp[i / 2][1]),
                    Payload::new(precomp[i / 2][2]),
                );
                let temp = p.double();
                precomp[i][0] = temp.0.data();
                precomp[i][1] = temp.1.data();
                precomp[i][2] = temp.2.data();

                let p = P256JacobianPoint(
                    Payload::new(precomp[i][0]),
                    Payload::new(precomp[i][1]),
                    Payload::new(precomp[i][2]),
                );
                let temp = p.add_affine(&self);
                precomp[i + 1][0] = temp.0.data();
                precomp[i + 1][1] = temp.1.data();
                precomp[i + 1][2] = temp.2.data();

                i += 2;
            }
            precomp
        };

        let scalar = w_naf(scalar);
        let mut n_is_infinity_mask = u32::MAX;
        let mut counter = 0u16;

        let mut p1 = P256JacobianPoint(
            Payload::init(), Payload::init(), Payload::init(),
        );
        // let mut p2 = P256JacobianPoint(
        //     Payload::init(), Payload::init(), Payload::init(),
        // );

        for i in 0..scalar.len() {
            if scalar[i] == 0 {
                counter += 1;
                continue;
            }
            while counter > 0 {
                p1 = p1.double();
                counter -= 1;
            }

            let idx = (scalar[i].abs()) as u32;
            p1 = p1.double();
            let p2 = P256JacobianPoint::select(idx, points);

            let p3 = {
                if scalar[i] > 0 {
                    p1.add(&p2)
                } else {
                    p1.subtract(&p2)
                }
            };

            p1 = p1.copy_from_with_conditional(p2, n_is_infinity_mask);
            let p_is_finite_mask = mask(idx);
            let msk = p_is_finite_mask & !(n_is_infinity_mask);
            p1 = p1.copy_from_with_conditional(p3, msk);
            n_is_infinity_mask &= !(p_is_finite_mask);
        }

        while counter > 0 {
            p1 = p1.double();
            counter -= 1;
        }

        p1.to_affine_point()
    }
}


/// 基点
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub(crate) struct P256BasePoint {
    point: P256AffinePoint,
    order: BigUint,
}

impl P256BasePoint {
    pub(crate) fn new(point: P256AffinePoint, order: BigUint) -> Self {
        P256BasePoint { point, order }
    }
}

impl Multiplication for P256BasePoint {
    /// multiply sets P256Point = scalar*G where scalar is a little-endian number.
    fn multiply(&self, scalar: BigUint) -> P256AffinePoint {
        let scalar = {
            let mut bytes = [0u8; 32];
            for (i, v) in scalar.to_bytes_le().iter().enumerate() {
                bytes[i] = *v;
            }
            bytes
        };

        let mut jacobian = P256JacobianPoint(
            Payload::init(), Payload::init(), Payload::init(),
        );

        let mut n_is_infinity_mask = !(0 as u32);   // u32::MAX
        // The loop adds bits at positions 0, 64, 128 and 192, followed by positions 32, 96, 160
        // and 224 and does this 32 times.
        for i in 0..32 {
            if i != 0 {
                jacobian = jacobian.double();
            }
            let mut offset = 0;
            let mut j = 0;
            while j <= 32 {
                let bit0 = bit_of_scalar(scalar, 31 - i + j);
                let bit1 = bit_of_scalar(scalar, 95 - i + j);
                let bit2 = bit_of_scalar(scalar, 159 - i + j);
                let bit3 = bit_of_scalar(scalar, 223 - i + j);
                let idx = bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3);

                let affine = P256AffinePoint::select(
                    idx,
                    Vec::from(&BASE_TABLE[offset..]),
                );

                offset += 30 * 9;

                let temp = jacobian.add_affine(&affine);
                jacobian = jacobian.copy_from_with_conditional(
                    P256JacobianPoint(
                        affine.0.clone(),
                        affine.1.clone(),
                        Payload::new(P256FACTOR[1]),
                    ),
                    n_is_infinity_mask,
                );

                let p_is_finite_mask = mask(idx);
                let mask = p_is_finite_mask & !n_is_infinity_mask;

                jacobian = jacobian.copy_from_with_conditional(temp, mask);

                // If p was not zero, then n is now non-zero.
                n_is_infinity_mask = n_is_infinity_mask & !p_is_finite_mask;

                j += 32;
            }
        }
        jacobian.to_affine_point()
    }
}

/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
#[derive(Copy, Clone, Debug)]
pub(crate) struct P256JacobianPoint(Payload, Payload, Payload);


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
    pub(crate) fn add_affine(&self, affine: &P256AffinePoint) -> Self {
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
    fn copy_from_with_conditional(&self, source: P256JacobianPoint, mask: u32) -> Self {
        let (mut x, mut y, mut z) = (
            Payload::init().data(), Payload::init().data(), Payload::init().data()
        );
        for i in 0..9 {
            x[i] = self.0.data()[i] ^ (mask & (source.0.data()[i] ^ self.0.data()[i]));
            y[i] = self.1.data()[i] ^ (mask & (source.1.data()[i] ^ self.1.data()[i]));
            z[i] = self.2.data()[i] ^ (mask & (source.2.data()[i] ^ self.2.data()[i]));
        }
        P256JacobianPoint(
            Payload::new(x),
            Payload::new(y),
            Payload::new(z),
        )
    }

    /// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
    /// Affine coordinates: (X = x/z^2, Y = y/z^3)  Y^2 = X^3 + aX +b
    pub(crate) fn to_affine_point(&self) -> P256AffinePoint {
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

    /// get the entry of table by index.
    /// On entry: index < 16, table[0] must be zero.
    fn select(index: u32, table: [[[u32; 9]; 3]; 16]) -> Self {
        let (mut x, mut y, mut z) = ([0u32; 9], [0u32; 9], [0u32; 9]);
        // The implicit value at index 0 is all zero.
        // We don't need to perform that iteration of the loop because we already set out_* to zero.
        for i in 0..16 {
            let mut mask = i ^ index;
            mask |= mask >> 2;
            mask |= mask >> 1;
            mask &= 1;
            mask = mask.wrapping_sub(1);

            for j in 0..9 {
                x[j] |= table[i as usize][0][j] & mask;
                y[j] |= table[i as usize][1][j] & mask;
                z[j] |= table[i as usize][2][j] & mask;
            }
        }

        let x = Payload::new(x);
        let y = Payload::new(y);
        let z = Payload::new(z);

        P256JacobianPoint(x, y, z)
    }

    /// (x3, y3, z3) = (x1, y1, z1) + (x2, y2, z2)
    ///
    /// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    fn add(&self, other: &P256JacobianPoint) -> Self {
        let (x1, y1, z1) = (&self.0, &self.1, &self.2);
        let (x2, y2, z2) = (&other.0, &other.1, &other.2);

        // z1 = 0
        if let Sign::NoSign = PayloadHelper::restore(z1).sign() {
            return P256JacobianPoint(x2.clone(), y2.clone(), z2.clone());
        }
        // z2 = 0
        if let Sign::NoSign = PayloadHelper::restore(z2).sign() {
            return P256JacobianPoint(x1.clone(), y1.clone(), z1.clone());
        }

        let z12 = z1.square();
        let z22 = z2.square();

        let z13 = z12.multiply(z1);
        let z23 = z22.multiply(z2);

        // u1 = x1 * z2^2  u2 = x2 * z1^2
        let u1 = x1.multiply(&z22);
        let u2 = x2.multiply(&z12);

        // s1 = y1 * z2^3  s2 = y2 * z1^3
        let s1 = y1.multiply(&z23);
        let s2 = y2.multiply(&z13);

        let h = u2.subtract(&u1);
        let r = s2.subtract(&s1);

        let r2 = r.square();
        let h2 = h.square();
        let h3 = h2.multiply(&h);

        let tmp = u1.multiply(&h2);

        let x3 = r2.subtract(&h2.multiply(&h)).subtract(&tmp.scalar_multiply(2));
        let y3 = r.multiply(&tmp.subtract(&x3)).subtract(&h3.multiply(&s1));
        let z3 = {
            let u1b = PayloadHelper::restore(&u1);
            let u2b = PayloadHelper::restore(&u2);
            let s1b = PayloadHelper::restore(&s1);
            let s2b = PayloadHelper::restore(&s2);

            if Ordering::Equal == u1b.cmp(&u2b) && Ordering::Equal == s1b.cmp(&s2b) {
                let p = self.double();
                p.2.multiply(&z2).multiply(&h)
            } else {
                z1.multiply(&z2).multiply(&h)
            }
        };

        P256JacobianPoint(x3, y3, z3)
    }

    /// (x3, y3, z3) = (x1, y1, z1) - (x2, y2, z2)
    fn subtract(&self, other: &P256JacobianPoint) -> Self {
        let another = P256JacobianPoint(
            other.0.clone(),
            PayloadHelper::transform(&PayloadHelper::restore(&other.1).neg()),
            other.2.clone(),
        );
        self.add(&another)
    }
}


#[inline(always)]
fn bit_of_scalar(scalar: [u8; 32], bit: usize) -> u32 {
    (((scalar[bit >> 3]) >> (bit & 7)) & 1) as u32
}

#[inline(always)]
fn w_naf(scalar: BigUint) -> Vec<i8> {
    let mut k = scalar;

    let bits = k.bits() as usize;
    let mut naf: Vec<i8> = vec![0; bits + 1];

    if let Sign::NoSign = k.to_bigint().unwrap().sign() {
        return naf;
    }

    let mut carry = false;
    let mut length: usize = 0;
    let mut pos: usize = 0;

    while pos <= bits {
        let s = k.clone().shr(pos).bitand(BigUint::from(1u64));
        if s.to_usize().unwrap() == (carry as usize) {
            pos += 1;
            continue;
        }
        k = k.shr(pos);
        let mask = BigUint::from(15usize);
        let mut digit: isize = k.clone().bitand(mask).to_isize().unwrap();
        if carry {
            digit += 1;
        }
        carry = (digit & 8) != 0;
        if carry {
            digit -= 16;
        }
        length += pos;
        naf[length] = digit as i8;
        pos = 4usize;
    }

    if naf.len() > length + 1 {
        let mut t = vec![0; length + 1];
        for (d, s) in t.iter_mut().zip(naf[0..(length + 1)].iter()) {
            *d = *s;
        }
        naf = t
    }
    naf.reverse();
    naf
}

#[cfg(test)]
mod tests {
    use num_traits::Num;

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

    #[test]
    fn to_affine() {
        let jacobian = P256JacobianPoint(
            Payload::new([302587400, 224711462, 627912361, 12505049, 498636470, 226242352, 402285030, 277184676, 216966475]),
            Payload::new([192016430, 212978101, 582317843, 172876572, 311643684, 126400666, 241514474, 362965479, 507691953]),
            Payload::new([186636191, 229928314, 430146881, 262724875, 500465416, 219885119, 175182585, 128499041, 217581763]),
        );
        let x: [u32; 9] = [194013013, 230698553, 317844872, 128801727, 111436768, 164685344, 76578606, 217356592, 311205467];
        let y: [u32; 9] = [26049626, 112805900, 275795042, 259495837, 289529507, 146296588, 220416178, 146512122, 266185762];

        let p = jacobian.to_affine_point();
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

    #[test]
    fn add_jacobian() {
        let p1 = P256JacobianPoint(
            Payload::new([536870910, 268435455, 536871167, 268433407, 536870911, 268435455, 536870911, 234881023, 536870911]),
            Payload::new([536870910, 268435455, 536871167, 268433407, 536870911, 268435455, 536870911, 234881023, 536870911]),
            Payload::new([536870910, 268435455, 536871167, 268433407, 536870911, 268435455, 536870911, 234881023, 536870911]),
        );

        let p2 = P256JacobianPoint(
            Payload::new([213941498, 21300983, 60022125, 97060820, 192974655, 35884974, 326765193, 113910449, 256521185]),
            Payload::new([57250121, 220765648, 315404192, 140781057, 276132260, 27646902, 354194608, 33763371, 49435241]),
            Payload::new([2, 0, 536870656, 2047, 0, 0, 0, 33554432, 0]),
        );

        let p3 = p1.add(&p2);

        assert_eq!(p3.0.data(), [213941498, 21300983, 60022125, 97060820, 192974655, 35884974, 326765193, 113910449, 256521185]);
        assert_eq!(p3.1.data(), [57250121, 220765648, 315404192, 140781057, 276132260, 27646902, 354194608, 33763371, 49435241]);
        assert_eq!(p3.2.data(), [2, 0, 536870656, 2047, 0, 0, 0, 33554432, 0]);
    }

    #[test]
    fn sub_jacobian() {
        let p1 = P256JacobianPoint(
            Payload::new([169820625, 110064376, 930792524, 208967388, 469280954, 48848243, 81769946, 205036988, 325110895]),
            Payload::new([118502522, 33386085, 620331139, 260737378, 417984658, 26729668, 432554301, 226329341, 529338833]),
            Payload::new([258069116, 208474624, 937274820, 247524354, 178408971, 252379317, 302659768, 243985694, 317539455]),
        );

        let p2 = P256JacobianPoint(
            Payload::new([509832793, 181703115, 787976024, 210478362, 520855061, 159289624, 348831704, 240602015, 509540099]),
            Payload::new([271814258, 176524652, 698366764, 34633944, 357938607, 172608771, 211149456, 123200817, 417183214]),
            Payload::new([397113909, 50459711, 307720741, 153235104, 7288360, 139155237, 213953529, 206075123, 230363768]),
        );

        let p3 = p1.subtract(&p2);

        assert_eq!(p3.0.data(), [295090358, 236992739, 800053525, 147234841, 281370475, 197897281, 305280418, 361835277, 162960459]);
        assert_eq!(p3.1.data(), [229663282, 133755872, 656501873, 17946166, 228212011, 56920858, 395700549, 125368282, 347100819]);
        assert_eq!(p3.2.data(), [234698535, 154439292, 363189331, 134307834, 513337116, 113297570, 189927841, 204178274, 333316045]);
    }

    #[test]
    fn point_multiply() {
        let scalar = BigUint::from_str_radix("52097475535247475123296179337062319910931289617245574116042610944477699996763", 10).unwrap();
        let p = P256AffinePoint::new(
            Payload::new([213941498, 21300983, 60022125, 97060820, 192974655, 35884974, 326765193, 113910449, 256521185]),
            Payload::new([57250121, 220765648, 315404192, 140781057, 276132260, 27646902, 354194608, 33763371, 49435241]),
        );

        let point = p.multiply(scalar);
        let (x, y) = (PayloadHelper::restore(&point.0), PayloadHelper::restore(&point.1));

        let rx = BigUint::from_str_radix("86785762605949049655217626947696016468061046596615849956268278837626813985762", 10).unwrap();
        let ry = BigUint::from_str_radix("43927620980048020323390743598648642157054712812706113912175335564413059446606", 10).unwrap();

        assert_eq!(x.to_biguint().unwrap(), rx);
        assert_eq!(y.to_biguint().unwrap(), ry);
    }
}