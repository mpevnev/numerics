use num::{Float, FromPrimitive, Signed, abs};

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
    where T: Float + FromPrimitive + Signed,
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
    let mut mid_val = target(mid);
    let max = config.max_iters;
    while right - left > config.precision && max.map_or(true, |m| iter < m) {
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
        mid_val = target(mid);
    }

    if abs(left_val) < abs(mid_val) {
        Some(left)
    } else if abs(right_val) < abs(mid_val) {
        Some(right)
    } else {
        Some(mid)
    }
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

pub fn bisect_multi<'a, T: 'a, F>(config: MultiRootBisectCfg<T>,
                              left: T,
                              right: T,
                              target: &'a F)
    -> impl Iterator<Item = T> + 'a
    where T: Float + FromPrimitive + Epsilon<RHS=T, Precision=T> + Signed,
          F: Fn(T) -> T
{
    MultiRootBisectState {
        cfg: config,
        left,
        right,
        target,
        cur_interval: 0,
        last_root: None
    }
}

impl<'a, T, F> Iterator for MultiRootBisectState<'a, T, F>
    where T: Float + FromPrimitive + Signed + Epsilon<RHS=T, Precision=T>,
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
                self.last_root = Some(root);
                return Some(root)
            }
        }
        None
    }
}

/* ---------- Newton's method ---------- */

/// Configuration structure for the Newton's method (one root version).
#[derive(Debug, Clone, Copy)]
pub struct OneRootNewtonCfg<T> {
    /// The real root, if any, is most likely to be within this distance from
    /// the reported root, but this is not guaranteed.
    pub precision: T,
    /// A limit on the number of iterations to perform. Pass `None` if you
    /// don't want a limit.
    pub max_iters: Option<u32>
}

pub fn newton_one<T, F, D>(config: OneRootNewtonCfg<T>,
                           left: T,
                           right: T,
                           first_approx: T,
                           target: &F,
                           derivative: &D)
    -> Option<T>
    where T: Float + Epsilon<RHS=T, Precision=T>,
          F: Fn(T) -> T,
          D: Fn(T) -> T
{
    let mut left = left;
    let mut right = right;
    let mut root = first_approx;
    let mut prev_root = None;
    let mut iter = 0;
    while prev_root.map_or(true, |old| !root.close(old, config.precision))
        && config.max_iters.map_or(true, |max| iter < max) {
            iter += 1;
            let left_val = target(left);
            let right_val = target(right);
            if let Some(next) = next_newton_iter(config.precision,
                                                 left, 
                                                 right, 
                                                 root, 
                                                 target, 
                                                 derivative) {
                prev_root = Some(root);
                root = next;
            } else if let Some(fallback_root) 
                = linear_fallback(left, right, left_val, right_val) {
                    prev_root = Some(root);
                    root = fallback_root;
            } else {
                return None
            }
            let val_at_root = target(root);
            if left_val * val_at_root <= T::zero() {
                right = root;
            } else {
                left = root;
            }
    }
    Some(root)
}

fn next_newton_iter<T, F, D>(prec: T,
                             left: T,
                             right: T,
                             old: T,
                             target: &F,
                             derivative: &D)
    -> Option<T>
    where T: Float + Epsilon<RHS=T, Precision=T>,
          F: Fn(T) -> T,
          D: Fn(T) -> T
{
    let d = derivative(old);
    if d.near_zero(prec) {
        return None
    }
    let res = old - target(old) / d;
    if res < left {
        None
    } else if res > right {
        None
    } else {
        Some(res)
    }
}

fn linear_fallback<T: Float>(x1: T , x2: T, y1: T, y2: T) -> Option<T>
{
    let res = ((y2 - y1) * x1 - (x2 - x1) * y1) / (y2 - y1);
    if res < x1 {
        None
    } else if res > x2 {
        None
    } else {
        Some(res)
    }
}

/* ---------- unit tests ---------- */

#[cfg(test)]
mod tests {
    test_suite! {
        name bisect_one_tests;

        use num::pow::Pow;

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
            let roots: Vec<_> = bisect_multi(cfg, -3.0, 3.0, &target).collect();
            assert_that!(&roots.len(), eq(2));
            assert_that!(roots[0].close(-2.0, prec));
            assert_that!(roots[1].close(2.0, prec));
        }

        test bisect_multi_pos_2() {
            let target = |x| x;
            let prec = 1e-6;
            let cfg = MultiRootBisectCfg {
                precision: prec,
                max_iters: None,
                num_intervals: 2
            };
            let roots: Vec<_> = bisect_multi(cfg, -1.0, 1.0, &target).collect();
            assert_that!(&roots.len(), eq(1));
            assert_that!(roots[0].close(0.0, prec));
        }

        test bisect_multi_neg_1() {
            let target = |x| (x - 1.0) * (x - 2.0);
            let prec = 1e-6;
            let cfg = MultiRootBisectCfg {
                precision: prec,
                max_iters: None,
                num_intervals: 10
            };
            let roots = bisect_multi(cfg, 3.0, 4.0, &target).collect::<Vec<_>>();
            assert_that!(&roots, eq(vec![]));
        }

        test newton_one_pos_1() {
            let target = |x: f64| (x - 1.0) * (x - 2.0) * (x - 3.0);
            let der = |x: f64| 3.0 * x.pow(2) - 12.0 * x + 11.0;
            let prec = 1e-6;
            let cfg = OneRootNewtonCfg {
                precision: prec,
                max_iters: None
            };
            let root = newton_one(cfg, 0.5, 1.5, 0.55, &target, &der);
            assert_that!(root.unwrap().close(1.0, prec));
            let root = newton_one(cfg, 1.5, 2.5, 1.55, &target, &der);
            assert_that!(root.unwrap().close(2.0, prec));
            let root = newton_one(cfg, 2.5, 4.0, 3.15, &target, &der);
            assert_that!(root.unwrap().close(3.0, prec));
        }

        test newton_one_pos_2() {
            let target = |x: f64| x.pow(0.1) - 1.0;
            let der = |x: f64| 0.1 * x.pow(-0.9);
            let prec = 1e-6;
            let cfg = OneRootNewtonCfg {
                precision: prec,
                max_iters: None
            };
            let root = newton_one(cfg, 0.5, 1.5, 0.55, &target, &der);
            assert_that!(root.unwrap().close(1.0, prec));
        }

        test newton_one_neg_1() {
            let target = |x: f64| (x - 1.0) * (x - 2.0) * (x - 3.0);
            let der = |x: f64| 3.0 * x.pow(2) - 12.0 * x + 11.0;
            let prec = 1e-6;
            let cfg = OneRootNewtonCfg {
                precision: prec,
                max_iters: None
            };
            let root = newton_one(cfg, 5.0, 6.0, 5.5, &target, &der);
            assert_that!(root.is_none());
        }
    }
}
