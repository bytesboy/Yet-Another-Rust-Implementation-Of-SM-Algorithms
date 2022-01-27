use std::ops::{Add, Mul, Shl, Shr};

use num_bigint::{BigInt, ToBigInt};
use num_integer::Integer;
use num_traits::FromPrimitive;

use crate::sm2::p256::core::P256Elliptic;
use crate::sm2::p256::params::{P256CARRY, P256FACTOR, P256ZERO31};

/// Field elements are represented as nine, unsigned 32-bit words. The value of a field element is:
///
/// ```
/// Value = (x8 * 2^228) + (x7 * 2^200) + (x6 * 2^171) + (x5 * 2^143) + (x4 * 2^114) + (x3 * 2^86) +
///         (x2 * 2^57)  + (x1 * 2^29)  + x0
/// ```
///
/// That is, each limb is alternately 29 or 28-bits wide in little-endian order.
///
/// This means that a field element hits 2^257, rather than 2^256 as we would like.
/// A 28, 29, ... pattern would cause us to hit 2^256, but that causes problems
/// when multiplying as terms end up one bit short of a limb
/// which would require much bit-shifting to correct.
///
/// ***
/// ### Pattern
/// * bits-257: |29bits|28bits|29bits|29bits|28bits|29bits|29bits|28bits|29bits|
///
/// Finally, the values stored in a field element are in Montgomery form.
/// So the value |y| is stored as (y*R) mod p, where p is the P-256 prime and R is 2^257.
/// ***
/// ### little-endian order
///
/// Example: payload = \[x0, x1, x2, x3, x4, x5, x6, x7, x8]
///
/// | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits |
/// | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
/// |   x8   |   x7   |   x6   |   x5   |   x4   |   x3   |   x2   |   x1   |   x0   |
///
/// 0xFFFFFFF  = 1111111111111111111111111111
/// 0x1FFFFFFF = 11111111111111111111111111111
enum LimbPattern {
    WIDTH28BITS = 0xFFFFFFF,
    WIDTH29BITS = 0x1FFFFFFF,
}

#[derive(Debug)]
pub(crate) struct Payload {
    data: [u32; 9],
}

impl Payload {
    pub(crate) fn init() -> Self {
        Payload { data: [0u32; 9] }
    }

    pub(crate) fn new(data: [u32; 9]) -> Self {
        Payload { data }
    }

    pub(crate) fn data(&self) -> [u32; 9] {
        self.data
    }

    /// payload3 = payload1 + payload2
    ///
    /// payload1 = \[x0, x1, x2, x3, x4, x5, x6, x7, x8]
    /// payload2 = \[y0, y1, y2, y3, y4, y5, y6, y7, y8]
    ///
    /// |capacity| 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits |
    /// |        |   x8   |   x7   |   x6   |   x5   |   x4   |   x3   |   x2   |   x1   |   x0   |
    /// |        |   y8   |   y7   |   y6   |   y5   |   y4   |   y3   |   y2   |   y1   |   y0   |
    ///
    /// |capacity| 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits |
    /// | length | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits |
    /// |   257  |   228  |   200  |   171  |   143  |   114  |   86   |   57   |   29   |    0   |
    /// | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
    /// | carry  |   r8   |   r7   |   r6   |   r5   |   r4   |   r3   |   r2   |   r1   |   r0   |
    ///
    /// On entry, payload1\[i] + payload2\[i] must not overflow a 32-bit word.
    /// On exit: payload3\[0,2,...] < 2^30, payload3\[1,3,...] < 2^29
    pub(crate) fn add(&self, other: &Payload) -> Payload {
        let mut result = Payload::init();
        let mut carry: u32 = 0;
        let mut i = 0;
        loop {
            let x = self.data[i].wrapping_add(other.data[i]).wrapping_add(carry);
            carry = x.shr(29);
            result.data[i] = x & (LimbPattern::WIDTH29BITS as u32);
            i += 1;
            if i == 9 {
                break;
            }
            let x = self.data[i].wrapping_add(other.data[i]).wrapping_add(carry);
            carry = x.shr(28);
            result.data[i] = x & (LimbPattern::WIDTH28BITS as u32);
            i += 1;
        }
        PayloadHelper::reduce_carry(&mut result, carry as usize);
        result
    }

    /// payload3 = payload1 - payload2
    ///
    /// On entry: payload1\[0,2,...] < 2^30, payload1\[1,3,...] < 2^29 and
    ///           payload2\[0,2,...] < 2^30, payload2\[1,3,...] < 2^29.
    /// On exit:  payload3\[0,2,...] < 2^30, payload3\[1,3,...] < 2^29.
    pub(crate) fn subtract(&self, other: &Payload) -> Payload {
        let mut result = Payload::init();
        let mut carry: u32 = 0;
        let mut i = 0;
        loop {
            let x = self.data[i].wrapping_sub(other.data[i]).wrapping_add(P256ZERO31[i]).wrapping_add(carry);
            carry = x.shr(29);
            result.data[i] = x & (LimbPattern::WIDTH29BITS as u32);
            i += 1;
            if i == 9 {
                break;
            }
            let x = self.data[i].wrapping_sub(other.data[i]).wrapping_add(P256ZERO31[i]).wrapping_add(carry);
            carry = x.shr(28);
            result.data[i] = x & (LimbPattern::WIDTH28BITS as u32);
            i += 1;
        }
        PayloadHelper::reduce_carry(&mut result, carry as usize);
        result
    }

    /// multiply sets payload3 = payload1 * payload2.
    ///
    /// On entry: payload1\[0,2,...] < 2^30, payload1\[1,3,...] < 2^29 and
    ///           payload2\[0,2,...] < 2^30, payload2\[1,3,...] < 2^29.
    /// On exit:  payload3\[0,2,...] < 2^30, payload3\[1,3,...] < 2^29.
    pub(crate) fn multiply(&self, other: &Payload) -> Payload {
        let mut result = Payload::init();
        let mut tmp: [u64; 17] = [0; 17];
        tmp[0] = (self.data[0] as u64) * (other.data[0] as u64);
        tmp[1] = (self.data[0] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[0] as u64) << 0);
        tmp[2] = (self.data[0] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[0] as u64) << 0);
        tmp[3] = (self.data[0] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[0] as u64) << 0);
        tmp[4] = (self.data[0] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[0] as u64) << 0);
        tmp[5] = (self.data[0] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[0] as u64) << 0);
        tmp[6] = (self.data[0] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[0] as u64) << 0);
        tmp[7] = (self.data[0] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[0] as u64) << 0);
        tmp[8] = (self.data[0] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[0] as u64) << 0);
        tmp[9] = (self.data[1] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[1] as u64) << 0);
        tmp[10] = (self.data[2] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[2] as u64) << 0);
        tmp[11] = (self.data[3] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[3] as u64) << 0);
        tmp[12] = (self.data[4] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[4] as u64) << 0);
        tmp[13] = (self.data[5] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[5] as u64) << 0);
        tmp[14] = (self.data[6] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[6] as u64) << 0);
        tmp[15] = (self.data[7] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[7] as u64) << 0);
        tmp[16] = (self.data[8] as u64) * ((other.data[8] as u64) << 0);

        PayloadHelper::reduce_degree(&mut result, &mut tmp);
        result
    }

    pub(crate) fn square(&self) -> Payload {
        let mut result = Payload::init();
        let mut tmp: [u64; 17] = [0; 17];
        tmp[0] = (self.data[0] as u64) * (self.data[0] as u64);
        tmp[1] = (self.data[0] as u64) * ((self.data[1] as u64) << 1);
        tmp[2] = (self.data[0] as u64) * ((self.data[2] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[1] as u64) << 1);
        tmp[3] = (self.data[0] as u64) * ((self.data[3] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[2] as u64) << 1);
        tmp[4] = (self.data[0] as u64) * ((self.data[4] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[3] as u64) << 2) +
            (self.data[2] as u64) * (self.data[2] as u64);
        tmp[5] = (self.data[0] as u64) * ((self.data[5] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[4] as u64) << 1) +
            (self.data[2] as u64) * ((self.data[3] as u64) << 1);
        tmp[6] = (self.data[0] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[5] as u64) << 2) +
            (self.data[2] as u64) * ((self.data[4] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[3] as u64) << 1);
        tmp[7] = (self.data[0] as u64) * ((self.data[7] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[2] as u64) * ((self.data[5] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[4] as u64) << 1);
        tmp[8] = (self.data[0] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[7] as u64) << 2) +
            (self.data[2] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[5] as u64) << 2) +
            (self.data[4] as u64) * (self.data[4] as u64);
        tmp[9] = (self.data[1] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[2] as u64) * ((self.data[7] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[4] as u64) * ((self.data[5] as u64) << 1);
        tmp[10] = (self.data[2] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[7] as u64) << 2) +
            (self.data[4] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[5] as u64) * ((self.data[5] as u64) << 1);
        tmp[11] = (self.data[3] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[4] as u64) * ((self.data[7] as u64) << 1) +
            (self.data[5] as u64) * ((self.data[6] as u64) << 1);
        tmp[12] = (self.data[4] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[5] as u64) * ((self.data[7] as u64) << 2) +
            (self.data[6] as u64) * (self.data[6] as u64);
        tmp[13] = (self.data[5] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[6] as u64) * ((self.data[7] as u64) << 1);
        tmp[14] = (self.data[6] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[7] as u64) * ((self.data[7] as u64) << 1);
        tmp[15] = (self.data[7] as u64) * ((self.data[8] as u64) << 1);
        tmp[16] = (self.data[8] as u64) * (self.data[8] as u64);

        PayloadHelper::reduce_degree(&mut result, &mut tmp);
        result
    }

    pub(crate) fn scalar_multiply(&self, n: usize) -> Payload {
        let p = Payload { data: P256FACTOR[n] };
        self.multiply(&p)
    }
}

pub(crate) struct PayloadHelper;

impl PayloadHelper {
    /// ### Example
    ///
    /// n: 115792089210356248756420345214020892766250353991924191454421193933289684991996
    ///
    /// * step 1 :
    ///     ```
    ///     x = (n * 2^257) % p;
    ///     x = 115792089048596568753516506446018802244132569949625955944202853485549017104377
    ///       = 1111111111111111111111111111100011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111001000000000000000000000000000001101111111111111111111111111111111111111111111111111111111111111001
    ///     ```
    /// * step 2: while loop, every time extract 29bits and 28 bits.
    ///     + step 2.1.1: change x: bigint to x: vec\<u32>
    ///         ```
    ///         x = [4294967289, 4294967295, 6, 4294967289, 4294967295, 4294967295, 4294967295, 4294967288]
    ///         ```
    ///     + step 2.1.2: extract 29 bits of x using operator &, `x[0] & 0x1FFFFFFF`
    ///         ```
    ///         11111111111111111111111111111001
    ///         &  11111111111111111111111111111
    ///         =  11111111111111111111111111001
    ///         ```
    ///     + step 2.1.3: right shift 29 bits on purpose to delete the extracted 29 bits.
    ///        ```
    ///         x = 11111111111111111111111111111000111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110010000000000000000000000000000011011111111111111111111111111111111111
    ///        ```
    ///     + step 2.2.1: some operation like 2.1.1
    ///     + step 2.2.2: extract 28 bits of x using operator &, `x[0] & 0xFFFFFFF`
    ///     + step 2.2.3: right shift 28 bits.
    /// * step 3: get the result
    ///     ```
    ///     data = [536870905, 268435455, 895, 268428288, 536870911, 268435455, 536870911, 150994943, 268435455]
    pub(crate) fn transform(n: &BigInt) -> Payload {
        let elliptic = P256Elliptic::init();
        let mut data: [u32; 9] = [0; 9];
        // data = n * R mod p = n * 2^257 % p
        let mut x: BigInt = BigInt::shl(n.clone(), 257);
        x = x.mod_floor(&elliptic.ec.p.to_bigint().unwrap());
        let mut i: usize = 0;
        while i < 9 {
            // x -> [u32]
            let bits = x.to_u32_digits().1;
            if !bits.is_empty() {
                // extract 29 bits using operator &
                data[i] = bits[0] & (LimbPattern::WIDTH29BITS as u32);
            } else {
                data[i] = 0
            }
            // right shift 29 bits
            x = BigInt::shr(x, 29);
            i += 1;
            if i == 9 {
                break;
            }
            let bits = x.to_u32_digits().1;
            if !bits.is_empty() {
                // extract 28 bits using operator &
                data[i] = bits[0] & (LimbPattern::WIDTH28BITS as u32);
            } else {
                data[i] = 0
            }
            // right shift 28 bits
            x = BigInt::shr(x, 28);
            i += 1;
        }
        Payload { data }
    }

    /// Example: payload = \[x0, x1, x2, x3, x4, x5, x6, x7, x8]
    ///
    ///           n = x8
    /// * i=7  => n = n * 2^28 + x7 = x8 * 2^28 + x7
    /// * i=6  => n = n * 2^29 + x6 = x8 * 2^57 + x7 * 2^29 + x6
    /// * i=5  => n = n * 2^28 + x5 = x8 * 2^85 + x7 * 2^57 + x6 * 2^28 + x5
    /// * i=4  => n = n * 2^29 + x4 = x8 * 2^114 + x7 * 2^86 + x6 * 2^57 + x5 * 2^29 + x4
    /// * i=3  => n = n * 2^28 + x3 = x8 * 2^142 + x7 * 2^114 + x6 * 2^85 + x5 * 2^57 + x4 * 2^28 + x3
    /// * i=2  => n = n * 2^29 + x2 = x8 * 2^171 + x7 * 2^143 + x6 * 2^114 + x5 * 2^86 + x4 * 2^57 + x3 * 2^29 + x2
    /// * i=1  => n = n * 2^28 + x1 = x8 * 2^199 + x7 * 2^171 + x6 * 2^142 + x5 * 2^114 + x4 * 2^85 + x3 * 2^57 + x2 * 2^28 + x1
    /// * i=0  => n = n * 2^29 + x0 = x8 * 2^228 + x7 * 2^200 + x6 * 2^171 + x5 * 2^143 + x4 * 2^114 + x3 * 2^86 + x2 * 2^57 + x1 * 2^29 + x0
    ///
    /// return n * RI mod p
    pub(crate) fn restore(payload: &Payload) -> BigInt {
        let elliptic = P256Elliptic::init();
        let mut n = BigInt::from_u32(payload.data[8]).unwrap();
        let mut temp: BigInt;
        let mut i: isize = 7;
        while i >= 0 {
            // i & 1 = 0 => i is even, else i is odd
            if (i & 1) == 0 {
                // even index, n * 2^29
                n = n.shl(29);
            } else {
                // odd index, n * 2^28
                n = n.shl(28);
            }
            temp = BigInt::from_u32(payload.data[i as usize]).unwrap();
            n = n.add(temp);
            i -= 1;
        }
        // formula: data = n * R mod P  => n = data * RI mod p
        n = n.mul(elliptic.ri.to_bigint().unwrap());
        n = n.mod_floor(&elliptic.ec.p.to_bigint().unwrap());
        n
    }

    /// 0xffffffff for 0 < x <= 2^31  0xffffffff = 4294967295 = u32::MAX = 2^31 - 1
    /// 0 for x == 0 or x > 2^31.
    fn mask(x: u32) -> u32 {
        x.wrapping_sub(1).wrapping_shr(31).wrapping_sub(1)
    }

    /// reduce_carry adds a multiple of p in order to cancel |carry|,which is a term at 2^257.
    ///
    /// payload = \[r0, r1, r2, r3, r4, r5, r6, r7, r8]
    ///
    /// we can count Res = carry * 2^257 + r8 * 2^228 + r7 * 2^200 + r6 * 2^171 + r5 * 2^143 + r4 * 2^114 + r3 * 2^86 + r2 * 2^57 + r1 * 2^29 + r0,
    /// and carry * 2^257 could transform to array P256CARRY, through with PayloadHelper::transform(carry), where carry is from 0 to 7.
    ///
    /// So we can mark ResArray = payload + P256CARRY, and Res = PayloadHelper::restore(ResArray)
    ///
    /// |capacity| 32bits |    32bits   | 32bits | 32bits | 32bits |     32bits   |    32bits   | 32bits |    32bits   |
    /// | length | 29bits |   <=29bits  | 29bits | 28bits | 29bits |   <=29bits   |   <=30bits  | 28bits |   <=30bits  |
    /// | ------ | ------ | ----------- | ------ | ------ | ------ | -----------  | ----------- | ------ | ----------- |
    /// |        |   r8   | r7+T[c*9+7] |   r6   |   r5   |   r4   |  r3+T[c*9+3] | r2+T[c*9+2] |   r1   | r0+T[c*9+2] |
    ///
    /// On entry: carry < 2^3, payload\[0,2,...] < 2^29, payload\[1,3,...] < 2^28.
    /// On exit: payload\[0,2,..] < 2^30, payload\[1,3,...] < 2^29.
    fn reduce_carry(payload: &mut Payload, carry: usize) {
        payload.data[0] += P256CARRY[carry * 9 + 0];
        payload.data[2] += P256CARRY[carry * 9 + 2];
        payload.data[3] += P256CARRY[carry * 9 + 3];
        payload.data[7] += P256CARRY[carry * 9 + 7];
    }

    /// reduce_degree sets a = b/R mod p where b contains 64-bit words with the same
    /// 29,28,... bit positions as a field element.
    ///
    /// The values in field elements are in Montgomery form: x*R mod p where R = 2^257.
    /// Since we just multiplied two Montgomery values together, the result is x * y * R * R mod p.
    /// We wish to divide by R in order for the result also to be in Montgomery form.
    ///
    /// On entry: tmp\[i] < 2^64
    /// On exit:  a\[0,2,...] < 2^30, a\[1,3,...] < 2^29
    ///
    /// Limb number:   0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10...
    /// Width (bits):  29| 28| 29| 28| 29| 28| 29| 28| 29| 28| 29
    /// Start bit:     0 | 29| 57| 86|114|143|171|200|228|257|285
    /// (odd phase):   0 | 28| 57| 85|114|142|171|199|228|256|285
    fn reduce_degree(a: &mut Payload, b: &mut [u64; 17]) {
        let mut tmp: [u32; 18] = [0; 18];
        let mut carry: u32;
        let mut x: u32;
        let mut x_mask: u32;

        tmp[0] = (b[0] as u32) & (LimbPattern::WIDTH29BITS as u32);
        tmp[1] = (b[0] as u32) >> 29;
        tmp[1] |= (((b[0] >> 32) as u32) << 3) & (LimbPattern::WIDTH28BITS as u32);
        tmp[1] += (b[1] as u32) & (LimbPattern::WIDTH28BITS as u32);
        carry = tmp[1] >> 28;
        tmp[1] &= LimbPattern::WIDTH28BITS as u32;

        let mut i = 2;
        while i < 17 {
            tmp[i] = ((b[i - 2] >> 32) as u32) >> 25;
            tmp[i] += ((b[i - 1]) as u32) >> 28;
            tmp[i] += (((b[i - 1] >> 32) as u32) << 4) & (LimbPattern::WIDTH29BITS as u32);
            tmp[i] += (b[i] as u32) & (LimbPattern::WIDTH29BITS as u32);
            tmp[i] += carry;
            carry = tmp[i] >> 29;
            tmp[i] &= LimbPattern::WIDTH29BITS as u32;

            i += 1;
            if i == 17 {
                break;
            }

            tmp[i] = ((b[i - 2] >> 32) as u32) >> 25;
            tmp[i] += (b[i - 1] as u32) >> 29;
            tmp[i] += (((b[i - 1] >> 32) as u32) << 3) & (LimbPattern::WIDTH28BITS as u32);
            tmp[i] += (b[i] as u32) & (LimbPattern::WIDTH28BITS as u32);
            tmp[i] += carry;
            carry = tmp[i] >> 28;
            tmp[i] &= LimbPattern::WIDTH28BITS as u32;

            i += 1
        }

        tmp[17] = ((b[15] >> 32) as u32) >> 25;
        tmp[17] += (b[16] as u32) >> 29;
        tmp[17] += ((b[16] >> 32) as u32) << 3;
        tmp[17] += carry;

        i = 0;
        loop {
            tmp[i + 1] += tmp[i] >> 29;
            x = tmp[i] & (LimbPattern::WIDTH29BITS as u32);
            tmp[i] = 0;

            if x > 0 {
                let mut set4: u32 = 0;
                let mut set7: u32 = 0;
                x_mask = Self::mask(x);
                tmp[i + 2] += (x << 7) & (LimbPattern::WIDTH29BITS as u32);
                tmp[i + 3] += x >> 22;

                if tmp[i + 3] < 0x10000000 {
                    set4 = 1;
                    tmp[i + 3] += 0x10000000 & x_mask;
                    tmp[i + 3] -= (x << 10) & (LimbPattern::WIDTH28BITS as u32);
                } else {
                    tmp[i + 3] -= (x << 10) & (LimbPattern::WIDTH28BITS as u32);
                }
                if tmp[i + 4] < 0x20000000 {
                    tmp[i + 4] += 0x20000000 & x_mask;
                    tmp[i + 4] -= set4;
                    tmp[i + 4] -= x >> 18;
                    if tmp[i + 5] < 0x10000000 {
                        tmp[i + 5] += 0x10000000 & x_mask;
                        tmp[i + 5] -= 1;
                        if tmp[i + 6] < 0x20000000 {
                            set7 = 1;
                            tmp[i + 6] += 0x20000000 & x_mask;
                            tmp[i + 6] -= 1;
                        } else {
                            tmp[i + 6] -= 1;
                        }
                    } else {
                        tmp[i + 5] -= 1;
                    }
                } else {
                    tmp[i + 4] -= set4;
                    tmp[i + 4] -= x >> 18;
                }

                if tmp[i + 7] < 0x10000000 {
                    tmp[i + 7] += 0x10000000 & x_mask;
                    tmp[i + 7] -= set7;
                    tmp[i + 7] -= (x << 24) & (LimbPattern::WIDTH28BITS as u32);
                    tmp[i + 8] += (x << 28) & (LimbPattern::WIDTH29BITS as u32);
                    if tmp[i + 8] < 0x20000000 {
                        tmp[i + 8] += 0x20000000 & x_mask;
                        tmp[i + 8] -= 1;
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += ((x >> 1) - 1) & x_mask;
                    } else {
                        tmp[i + 8] -= 1;
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += (x >> 1) & x_mask;
                    }
                } else {
                    tmp[i + 7] -= set7;
                    tmp[i + 7] -= (x << 24) & (LimbPattern::WIDTH28BITS as u32);
                    tmp[i + 8] += (x << 28) & (LimbPattern::WIDTH29BITS as u32);
                    if tmp[i + 8] < 0x20000000 {
                        tmp[i + 8] += 0x20000000 & x_mask;
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += ((x >> 1) - 1) & x_mask;
                    } else {
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += (x >> 1) & x_mask;
                    }
                }
            }

            if (i + 1) == 9 {
                break;
            }
            tmp[i + 2] += tmp[i + 1] >> 28;
            x = tmp[i + 1] & (LimbPattern::WIDTH28BITS as u32);
            tmp[i + 1] = 0;

            if x > 0 {
                let mut set5 = 0;
                let mut set8 = 0;
                let mut set9 = 0;
                x_mask = Self::mask(x);
                tmp[i + 3] += (x << 7) & (LimbPattern::WIDTH28BITS as u32);
                tmp[i + 4] += x >> 21;

                if tmp[i + 4] < 0x20000000 {
                    set5 = 1;
                    tmp[i + 4] += 0x20000000 & x_mask;
                    tmp[i + 4] -= (x << 11) & (LimbPattern::WIDTH29BITS as u32);
                } else {
                    tmp[i + 4] -= (x << 11) & (LimbPattern::WIDTH29BITS as u32);
                }
                if tmp[i + 5] < 0x10000000 {
                    tmp[i + 5] += 0x10000000 & x_mask;
                    tmp[i + 5] -= set5;
                    tmp[i + 5] -= x >> 18;
                    if tmp[i + 6] < 0x20000000 {
                        tmp[i + 6] += 0x20000000 & x_mask;
                        tmp[i + 6] -= 1;
                        if tmp[i + 7] < 0x10000000 {
                            set8 = 1;
                            tmp[i + 7] += 0x10000000 & x_mask;
                            tmp[i + 7] -= 1;
                        } else {
                            tmp[i + 7] -= 1;
                        }
                    } else {
                        tmp[i + 6] -= 1;
                    }
                } else {
                    tmp[i + 5] -= set5;
                    tmp[i + 5] -= x >> 18;
                }

                if tmp[i + 8] < 0x20000000 {
                    set9 = 1;
                    tmp[i + 8] += 0x20000000 & x_mask;
                    tmp[i + 8] -= set8;
                    tmp[i + 8] -= (x << 25) & (LimbPattern::WIDTH29BITS as u32);
                } else {
                    tmp[i + 8] -= set8;
                    tmp[i + 8] -= (x << 25) & (LimbPattern::WIDTH29BITS as u32);
                }
                if tmp[i + 9] < 0x10000000 {
                    tmp[i + 9] += 0x10000000 & x_mask;
                    tmp[i + 9] -= set9;
                    tmp[i + 9] -= x >> 4;
                    tmp[i + 10] += (x - 1) & x_mask;
                } else {
                    tmp[i + 9] -= set9;
                    tmp[i + 9] -= x >> 4;
                    tmp[i + 10] += x & x_mask;
                }
            }
            i += 2;
        }

        carry = 0;
        i = 0;
        while i < 8 {
            a.data[i] = tmp[i + 9];
            a.data[i] += carry;
            a.data[i] += (tmp[i + 10] << 28) & (LimbPattern::WIDTH29BITS as u32);
            carry = a.data[i] >> 29;
            a.data[i] &= LimbPattern::WIDTH29BITS as u32;

            i += 1;
            a.data[i] = tmp[i + 9] >> 1;
            a.data[i] += carry;
            carry = a.data[i] >> 28;
            a.data[i] &= LimbPattern::WIDTH28BITS as u32;

            i += 1;
        }

        a.data[8] = tmp[17];
        a.data[8] += carry;
        carry = a.data[8] >> 29;
        a.data[8] &= LimbPattern::WIDTH29BITS as u32;

        Self::reduce_carry(a, carry as usize)
    }
}


#[cfg(test)]
mod tests {
    use num_traits::Num;

    use super::*;

    #[test]
    fn payload() {
        let n = "115792089210356248756420345214020892766250353991924191454421193933289684991996";
        let n = BigInt::from_str_radix(n, 10).unwrap();
        let payload = PayloadHelper::transform(&n);
        assert_eq!(payload.data, [
            536870905, 268435455, 895,
            268428288, 536870911, 268435455,
            536870911, 150994943, 268435455
        ]);
        let m = PayloadHelper::restore(&payload);
        assert_eq!(m, n);
    }
}
