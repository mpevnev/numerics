use num::{Float, Signed, abs};

/// A trait for things that can be approximately equal.
pub trait Epsilon {
    type RHS;
    type PrecisionType;

    /// Return true if self and `other` differ no more than by a given amount.
    fn close(&self, other: Self::RHS, precision: Self::PrecisionType) -> bool;

    /// Return true if self is close to zero.
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
