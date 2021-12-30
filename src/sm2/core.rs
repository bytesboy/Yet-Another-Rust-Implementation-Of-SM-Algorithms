use std::cmp::Ordering;
use std::mem;
use std::sync::Once;

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::Num;

pub trait EllipticProvider {
    fn initialize() -> Self;
}

pub trait Multiplication {
    fn multiply(&self, scalar: BigUint) -> Point;
}


#[derive(Clone, Debug)]
pub struct Elliptic {
    // the order of the finite field
    p: BigUint,
    // the constant a and b of the curve equation: y^2 = x^3 + ax +b
    a: BigUint,
    b: BigUint,
    // the base point g(gx, gy) and order n
    g: BasePoint,
    bits: usize,
}


impl Elliptic {
    pub fn new(p: BigUint, a: BigUint, b: BigUint, g: BasePoint, bits: usize) -> Self {
        Elliptic { p, a, b, g, bits }
    }

    pub fn base(&self) -> &BasePoint {
        &self.g
    }

    pub fn bits(&self) -> usize {
        self.bits
    }

    pub fn p(&self) -> &BigUint {
        &self.p
    }
}


#[derive(Clone, Debug)]
pub struct Point(BigUint, BigUint);

impl Point {
    pub fn new(x: BigUint, y: BigUint) -> Self { Point(x, y) }
}

impl Multiplication for Point {
    fn multiply(&self, scalar: BigUint) -> Point {
        self.clone()
    }
}

#[derive(Clone, Debug)]
pub struct BasePoint {
    point: Point,
    order: BigUint,
}


impl BasePoint {
    pub fn new(point: Point, order: BigUint) -> Self { BasePoint { point, order } }

    pub fn point(&self) -> &Point {
        &self.point
    }

    pub fn order(&self) -> &BigUint {
        &self.order
    }
}


impl Multiplication for BasePoint {
    fn multiply(&self, scalar: BigUint) -> Point {
        let scalar = {
            // compare scalar and order, n = (scalar mod order) if scalar > order else scalar
            if let Ordering::Greater = scalar.cmp(&self.order) {
                scalar.mod_floor(&self.order)
            } else {
                scalar
            }
        };


        self.point.clone()
    }
}

/// 椭圆曲线上的点集记为E(Fp) = {(x,y,z)|x,y,z ∈ Fp且满足曲线方程y^2 = x^3 + axz^4 + bz^6}，
/// 其中a,b ∈ Fp， 且4a^3 +27b^2 != 0 mod p。
pub struct JacobPoint(BigUint, BigUint, BigUint);

impl JacobPoint {
    pub fn new(x: BigUint, y: BigUint, z: BigUint) -> Self {
        JacobPoint(x, y, z)
    }
}






