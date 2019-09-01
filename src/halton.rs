use crate::{traits::QRng, utils::primes};

const MAX_LOG_N: usize = 48;

/// One-dimensional Halton sequence generator with a given base.
#[derive(Clone)]
struct HaltonSeq1D {
    base: u32,
    digits: Vec<u32>,
    remainders: Vec<f64>,
    next_power: u64,
}

impl HaltonSeq1D {
    #[inline]
    fn new(base: u32) -> Self {
        Self { base, digits: vec![0], remainders: vec![0.], next_power: 1 }
    }

    #[inline]
    fn reset(&mut self) {
        self.digits.clear();
        self.digits.push(0);
        self.remainders.clear();
        self.remainders.push(0.);
        self.next_power = 1;
    }

    #[inline]
    unsafe fn next(&mut self, index: u64) -> f64 {
        // In order to avoid pre-allocating too much memory for digits and remainders,
        // we only extend digit/remainders vectors when the new bit appears. For a given
        // base, this happens on indices base^0, base^1, base^2, ...
        if index == self.next_power {
            self.digits.push(0);
            self.remainders.push(0.);
            self.next_power *= u64::from(self.base);
        }

        let base_f = f64::from(self.base);
        let mut digit = self.digits.as_mut_ptr();
        let rem = self.remainders.as_mut_ptr();

        // increase the least significant bit and see what happens
        *digit += 1;
        let h = if *digit == self.base {
            // handle carry over if it occurs
            let mut k = 0;
            while {
                k += 1;
                *digit = 0;
                digit = digit.add(1);
                *digit += 1;
                *digit == self.base
            } {}
            *rem.add(k - 1) = (f64::from(*digit) + *rem.add(k)) / base_f;
            for i in (1..k).rev() {
                *rem.add(i - 1) = *rem.add(i) / base_f;
            }
            *rem
        } else {
            // simple case, no carry
            f64::from(*digit) + *rem
        };
        h / base_f
    }
}

/// Halton low-discrepancy sequence generator.
///
/// The implementation follows "Fast, Portable and Reliable Algorithm for the
/// Calculation of Halton Numbers" (Kolar and Shea, 1993), where successive
/// elements are generated based on the previous ones. This makes the algorithm
/// an order of magnitude faster than the naive one.
///
/// The first `ndim` prime numbers are used as bases.
#[derive(Clone)]
pub struct HaltonSeq {
    index: u64,
    seqs: Vec<HaltonSeq1D>,
}

impl HaltonSeq {
    /// Returns a new Halton sequence generator with dimensionality `ndim`.
    #[inline]
    pub fn new(ndim: usize) -> Self {
        Self { index: 0, seqs: primes().take(ndim).map(|x| HaltonSeq1D::new(x as _)).collect() }
    }
}

impl QRng for HaltonSeq {
    #[inline]
    fn ndim(&self) -> usize {
        self.seqs.len()
    }

    #[inline]
    unsafe fn gen_fill_unchecked(&mut self, out: &mut [f64]) {
        if self.index >= (1 << MAX_LOG_N) {
            self.index = 0;
            self.seqs.iter_mut().for_each(HaltonSeq1D::reset);
        }
        self.index += 1;
        for (i, s) in self.seqs.iter_mut().enumerate() {
            *out.get_unchecked_mut(i) = s.next(self.index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::HaltonSeq;
    use crate::{utils::primes, QRng};

    const TOL: f64 = 1e-15;

    #[test]
    fn test_halton_seq_mod2() {
        let halton2 = vec![
            0.50000, // first cycle (zero excluded)
            0.25000, 0.75000, // second cycle
            0.12500, 0.62500, 0.37500, 0.87500, // third cycle
            // fourth cycle
            0.06250, 0.56250, 0.31250, 0.81250, 0.18750, 0.68750, 0.43750, 0.93750,
            // fifth cycle
            0.03125, 0.53125, 0.28125, 0.78125, 0.15625, 0.65625, 0.40625, 0.90625, 0.09375,
            0.59375, 0.34375, 0.84375, 0.21875, 0.71875, 0.46875, 0.96875,
        ];
        let mut seq = HaltonSeq::new(1).with_buf();
        for &x in halton2.iter() {
            assert!((x - seq.gen()[0]).abs() < TOL);
        }
    }

    #[test]
    fn test_halton_seq_mod3() {
        let halton3 = vec![
            // first cycle (zero excluded)
            1. / 3.,
            2. / 3.,
            // second cycle
            1. / 9.,
            4. / 9.,
            7. / 9.,
            2. / 9.,
            5. / 9.,
            8. / 9.,
            // third cycle
            1. / 27.,
            10. / 27.,
            19. / 27.,
            4. / 27.,
            13. / 27.,
            22. / 27.,
            7. / 27.,
            16. / 27.,
            25. / 27.,
            2. / 27.,
            11. / 27.,
            20. / 27.,
            5. / 27.,
            14. / 27.,
            23. / 27.,
            8. / 27.,
            17. / 27.,
            26. / 27.,
        ];
        let mut seq = HaltonSeq::new(2).with_buf();
        for &x in halton3.iter() {
            assert!((x - seq.gen()[1]).abs() < TOL);
        }
    }

    #[test]
    fn test_halton_seq_mean() {
        const MAX_DIM: usize = 7;
        for ndim in 1..=MAX_DIM {
            let p = primes().nth(ndim - 1).unwrap();
            let (mut sum, mut count) = (0., 0);
            let mut seq = HaltonSeq::new(ndim).with_buf();
            let mut cycle_len = p - 1; // initial cycle length
            for _ in 0..4 {
                // run 4 full cycles for the selected dimension
                for _ in 0..cycle_len {
                    sum += seq.gen()[ndim - 1];
                    count += 1;
                }
                cycle_len *= p; // cycle length is (p - 1) * p^n
            }
            let mean = sum / (count as f64);
            assert!((mean - 0.5).abs() < TOL);
        }
    }
}
