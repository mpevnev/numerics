use num::{Float, Signed, abs};

/// A trait for things that can be approximately equal.
pub trait Epsilon {
    type RHS;
    type Precision;

    /// Return true if self and `other` differ no more than by a given amount.
    fn close(&self, other: Self::RHS, precision: Self::Precision) -> bool;

    /// Return true if self is close to zero.
    fn near_zero(&self, precision: Self::Precision) -> bool;
}

impl<T: Float + Signed> Epsilon for T {
    type RHS = T;
    type Precision = T;

    fn close(&self, other: T, precision: T) -> bool {
        abs(other - *self) < abs(precision)
    }

    fn near_zero(&self, precision: T) -> bool {
        abs(*self) < abs(precision)
    }
}
