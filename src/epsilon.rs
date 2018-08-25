use num::{Float, Signed, abs};

pub trait Epsilon {
    type RHS;
    type PrecisionType;

    fn close(&self, other: Self::RHS, precision: Self::PrecisionType) -> bool;
    fn near_zero(&self, precision: Self::PrecisionType) -> bool;
}

impl<T: Float + Signed> Epsilon for T {
    type RHS = T;
    type PrecisionType = T;

    fn close(&self, other: T, precision: T) -> bool {
        abs(other - *self) < abs(precision)
    }

    fn near_zero(&self, precision: T) -> bool {
        abs(*self) < abs(precision)
    }
}
