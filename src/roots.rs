use num::{Float, FromPrimitive};

/* ---------- Functions to deal with single-root intervals ---------- */

pub struct OneRootConfig<T> {
    precision: T,
    max_iters: Option<u32>  // pass None to remove the limit.
}

pub fn bisect_one<T, F>(config: &OneRootConfig<T>,
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
    test_suite! {
        name bisect_one_tests;

        use epsilon::Epsilon;
        use roots::*;

        test bisect_one_pos_1() {
            let target = |x| x;
            let prec = 1e-6;
            let cfg = OneRootConfig { precision: prec, max_iters: None };
            let root = bisect_one(&cfg, -1.0, 1.0, target);
            assert_that!(root.unwrap().close(0.0, prec));
        }

        test bisect_one_pos_2() {
            let target = |x| (x - 2.0) * (x + 2.0);
            let prec = 1e-9;
            let cfg = OneRootConfig { precision: prec, max_iters: None };
            let root1 = bisect_one(&cfg, 1.8, 2.1, target);
            let root2 = bisect_one(&cfg, -10.0, 0.0, target);
            assert_that!(root1.unwrap().close(2.0, prec));
            assert_that!(root2.unwrap().close(-2.0, prec));
        }

        test bisect_one_neg_1() {
            let target = |x| x;
            let prec = 1e-6;
            let cfg = OneRootConfig { precision: prec, max_iters: None };
            let root = bisect_one(&cfg, 1.0, 2.0, target);
            assert_that!(&root, is_variant!(None));
        }
    }
}
