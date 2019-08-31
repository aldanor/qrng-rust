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
