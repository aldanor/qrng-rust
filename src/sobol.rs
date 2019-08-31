use crate::traits::QRng;

#[cfg(not(feature = "sobol-high-dim"))]
mod assets {
    pub type DirNum = u16;
    pub const SOBOL_MAX_DIM: usize = 1111;
    pub const SOBOL_LOCS: &[DirNum] = include!("assets/sobol-lim-locs.rs");
    pub const SOBOL_COEF: &[DirNum] = include!("assets/sobol-lim-coef.rs");
    pub const SOBOL_DNUM: &[DirNum] = include!("assets/sobol-lim-dnum.rs");
}
#[cfg(feature = "sobol-high-dim")]
mod assets {
    pub type DirNum = u32;
    pub const SOBOL_MAX_DIM: usize = 21201;
    pub const SOBOL_LOCS: &[DirNum] = include!("assets/sobol-all-locs.rs");
    pub const SOBOL_COEF: &[DirNum] = include!("assets/sobol-all-coef.rs");
    pub const SOBOL_DNUM: &[DirNum] = include!("assets/sobol-all-dnum.rs");
}
use self::assets::*;

const MAX_LOG_N: usize = 48;

#[inline]
fn get_raw_data(index: usize) -> (DirNum, &'static [DirNum]) {
    if index + 2 > SOBOL_MAX_DIM {
        panic!("invalid Sobol sequence dimension: {}", index + 2);
    }
    let (start, end) = (SOBOL_LOCS[index] as usize, SOBOL_LOCS[index + 1] as usize);
    (SOBOL_COEF[index], &SOBOL_DNUM[start..end])
}

unsafe fn get_dirnums(axis: usize, out: &mut [u64], stride: usize) {
    const MAX: usize = MAX_LOG_N - 1;
    let mut dirnums = [0; MAX_LOG_N];
    if axis == 0 {
        for (i, x) in dirnums.iter_mut().enumerate() {
            *x = 1 << (MAX - i);
        }
    } else {
        let (coef, m) = get_raw_data(axis - 1);
        let coef = u64::from(coef);
        let s = m.len() + 1;
        *dirnums.get_unchecked_mut(0) = 1 << MAX;
        for i in 1..s.min(MAX_LOG_N) {
            *dirnums.get_unchecked_mut(i) =
                (u64::from(*m.get_unchecked(i - 1)) * 2 + 1) << (MAX - i);
        }
        for i in s..MAX_LOG_N {
            let mut k = i - s;
            let dk = *dirnums.get_unchecked(k);
            let mut x = dk ^ (dk >> s);
            let mut coef_s = coef;
            k += 1;
            for _ in 1..s {
                x ^= (coef_s & 1) * *dirnums.get_unchecked(k);
                coef_s >>= 1;
                k += 1;
            }
            *dirnums.get_unchecked_mut(i) = x;
        }
    }
    for (i, &x) in dirnums.iter().enumerate() {
        *out.get_unchecked_mut(i * stride) = x;
    }
}

/// Sobol low-discrepancy sequence generator.
///
/// The implementation relies on primitive polynomials module two suggested in
/// "Constructing Sobol Sequences with Better Two-Dimensional Projections" (Joe and
/// Kuo, 2008).
#[derive(Clone)]
pub struct SobolSeq {
    ndim: usize,
    dirnums: Vec<u64>,
    value: Vec<u64>,
    index: u64,
}

impl SobolSeq {
    /// Returns a new Sobol sequence generator with dimensionality `ndim`.
    ///
    /// Panics if `ndim` is greater than 1111 (or 21201 if `sobol-high-dim` feature
    /// is enabled).
    #[inline]
    pub fn new(ndim: usize) -> Self {
        let mut dirnums = vec![0; ndim * MAX_LOG_N];
        for i in 0..ndim {
            unsafe { get_dirnums(i, &mut dirnums[i..], ndim) };
        }
        Self { ndim, dirnums, value: vec![0; ndim], index: 0 }
    }
}

impl QRng for SobolSeq {
    #[inline]
    fn ndim(&self) -> usize {
        self.ndim
    }

    #[inline]
    unsafe fn gen_fill_unchecked(&mut self, out: &mut [f64]) {
        const MAX_N: u64 = 1u64 << MAX_LOG_N;
        const DENUM: f64 = MAX_N as f64;
        let c = (!self.index).trailing_zeros() as usize;
        let v = self.dirnums.get_unchecked(c * self.ndim..);
        for j in 0..self.ndim {
            let x = self.value.get_unchecked_mut(j);
            *x ^= *v.get_unchecked(j);
            *out.get_unchecked_mut(j) = (*x as f64) / DENUM;
        }
        self.index = (self.index + 1) % MAX_N;
    }
}
