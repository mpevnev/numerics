use num::{Float, FromPrimitive};

use epsilon::Epsilon;

/* ---------- bisection for intervals with a single root ---------- */

/// Configuration structure for the bisection method (one root version).
#[derive(Debug, Clone, Copy)]
pub struct OneRootBisectCfg<T> {
    /// The real root, if any, will be no further than this from the reported
    /// root.
    pub precision: T,
    /// A limit on the number of iterations to perform. Pass `None` if you
    /// don't want a limit.
    pub max_iters: Option<u32>  
}

/// Find a root for a given function in a given interval, assuming there is
/// only one root there.
pub fn bisect_one<T, F>(config: OneRootBisectCfg<T>,
                        left: T, 
                        right: T, 
                        target: &F)
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

/* ---------- bisection for intervals with several roots ---------- */

/// Configuration structure for the bisection method (multiple roots version).
#[derive(Debug, Clone, Copy)]
pub struct MultiRootBisectCfg<T> {
    /// Real roots will be no further than that from reported roots.
    pub precision: T,
    /// A limit on the number of iterations to perform. Pass `None` if you 
    /// don't want a limit.
    pub max_iters: Option<u32>,
    /// The requested interval will be split into this many chunks, and each
    /// chunk will be tried for a root.
    pub num_intervals: usize
}

#[derive(Debug, Clone, Copy)]
struct MultiRootBisectState<'a, T, F: 'a> {
    cfg: MultiRootBisectCfg<T>, 
    left: T,
    right: T,
    target: &'a F,
    cur_interval: usize,
    last_root: Option<T>
}

/*
pub fn bisect_multi<'a, T, F>(config: MultiRootBisectCfg<T>,
                          left: T,
                          right: T,
                          target: F)
    -> impl Iterator<Item = T>
    where T: Float + FromPrimitive,
          F: 'a + Fn(T) -> T
{
    MultiRootBisectState {
        cfg: config,
        left,
        right,
        target: &target,
        cur_interval: 0
    }
}
*/

impl<'a, T, F> Iterator for MultiRootBisectState<'a, T, F> 
    where T: Float + FromPrimitive + Epsilon<RHS=T, PrecisionType=T>,
          F: 'a + Fn(T) -> T
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        if self.cur_interval > self.cfg.num_intervals {
            return None
        }
        let num_ints = T::from_usize(self.cfg.num_intervals)
            .expect("Failed to convert the number of intervals into a float");
        let interval_width = (self.right - self.left) / num_ints;
        while self.cur_interval < self.cfg.num_intervals {
            let int = T::from_usize(self.cur_interval)
                .expect("Failed to convert an index into a float");
            let left = self.left + interval_width * int;
            let right = left + interval_width;
            let one_cfg = OneRootBisectCfg { 
                precision: self.cfg.precision,
                max_iters: self.cfg.max_iters
            };
            let res = bisect_one(one_cfg, left, right, &self.target);
            self.cur_interval += 1;
            if let Some(root) = res {
                let two = T::from_i32(2).unwrap();
                let double_prec = two * self.cfg.precision;
                let mapper = |old: T| old.close(root, double_prec);
                let duplicate = self.last_root.map_or(false, mapper);
                if duplicate {
                    continue
                }
                return Some(root)
            }
        }
        None
    }
}

/* ---------- unit tests ---------- */

#[cfg(test)]
mod tests {
    test_suite! {
        name bisect_one_tests;

        use galvanic_assert::matchers::*;

        use epsilon::Epsilon;
        use roots::*;

        test bisect_one_pos_1() {
            let target = |x| x;
            let prec = 1e-6;
            let cfg = OneRootBisectCfg { precision: prec, max_iters: None };
            let root = bisect_one(cfg, -1.0, 1.0, &target);
            assert_that!(root.unwrap().close(0.0, prec));
        }

        test bisect_one_pos_2() {
            let target = |x| (x - 2.0) * (x + 2.0);
            let prec = 1e-9;
            let cfg = OneRootBisectCfg { precision: prec, max_iters: None };
            let root1 = bisect_one(cfg, 1.8, 2.1, &target);
            let root2 = bisect_one(cfg, -10.0, 0.0, &target);
            assert_that!(root1.unwrap().close(2.0, prec));
            assert_that!(root2.unwrap().close(-2.0, prec));
        }

        test bisect_one_neg_1() {
            let target = |x| x;
            let prec = 1e-6;
            let cfg = OneRootBisectCfg { precision: prec, max_iters: None };
            let root = bisect_one(cfg, 1.0, 2.0, &target);
            assert_that!(&root, is_variant!(None));
        }

        test bisect_multi_pos_1() {
            let target = |x| (x - 2.0) * (x + 2.0);
            let prec = 1e-6;
            let cfg = MultiRootBisectCfg {
                precision: prec,
                max_iters: None,
                num_intervals: 20
            };
            let state = MultiRootBisectState {
                cfg,
                left: -3.0,
                right: 3.0,
                target: &target,
                cur_interval: 0
            };
            let roots: Vec<_> = state.collect();
            assert_that!(&roots.len(), eq(2));
            assert_that!(roots[0].close(-2.0, prec));
            assert_that!(roots[1].close(2.0, prec));
        }
    }
}
