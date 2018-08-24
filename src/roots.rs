use num::{Float, FromPrimitive};

/* ---------- Functions to deal with single-root intervals ---------- */

pub struct OneRootConfig<T> {
    precision: T,
    max_iters: Option<u32>  // pass None to remove the limit.
}

pub fn bisect_one<T, F>(config: OneRootConfig<T>,
                        left: T, 
                        right: T, 
                        target: F)
    -> Option<T>
    where T: Float + FromPrimitive,
          F: Fn(T) -> T
{
    let mut iter = 0;
    let mut left = left;
    let mut right = right;
    let mut left_val = target(left);
    let mut right_val = target(right);

    if left_val * right_val > T::zero() {
        return None;
    }

    let mut mid = (left + right) / T::from_i32(2).unwrap();
    let max = config.max_iters;
    while right - left > config.precision && max.map_or(true, |m| iter < m) {
        let mid_val = target(mid);
        if left_val * mid_val <= T::zero() {
            right = mid;
            right_val = mid_val;
        } else if mid_val * right_val <= T::zero() {
            left = mid;
            left_val = mid_val;
        } else {
            return None;
        }
        iter += 1;
        mid = (left + right) / T::from_i32(2).unwrap();
    }

    Some(mid)
}

/* ---------- unit tests ---------- */

#[cfg(test)]
mod tests {

}
