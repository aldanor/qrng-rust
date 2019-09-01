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

#[cfg(test)]
mod tests {
    use super::{get_raw_data, SobolSeq};
    use crate::QRng;

    #[test]
    fn test_sobol_seq() {
        let expected = vec![
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            [0.75, 0.25, 0.25, 0.25, 0.75, 0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75, 0.75, 0.25, 0.25, 0.75, 0.25],
            [0.375, 0.375, 0.625, 0.875, 0.375, 0.125, 0.375, 0.875],
            [0.875, 0.875, 0.125, 0.375, 0.875, 0.625, 0.875, 0.375],
            [0.625, 0.125, 0.875, 0.625, 0.625, 0.875, 0.125, 0.125],
            [0.125, 0.625, 0.375, 0.125, 0.125, 0.375, 0.625, 0.625],
            [0.1875, 0.3125, 0.9375, 0.4375, 0.5625, 0.3125, 0.4375, 0.9375],
            [0.6875, 0.8125, 0.4375, 0.9375, 0.0625, 0.8125, 0.9375, 0.4375],
            [0.9375, 0.0625, 0.6875, 0.1875, 0.3125, 0.5625, 0.1875, 0.1875],
            [0.4375, 0.5625, 0.1875, 0.6875, 0.8125, 0.0625, 0.6875, 0.6875],
            [0.3125, 0.1875, 0.3125, 0.5625, 0.9375, 0.4375, 0.0625, 0.0625],
            [0.8125, 0.6875, 0.8125, 0.0625, 0.4375, 0.9375, 0.5625, 0.5625],
            [0.5625, 0.4375, 0.0625, 0.8125, 0.1875, 0.6875, 0.3125, 0.8125],
            [0.0625, 0.9375, 0.5625, 0.3125, 0.6875, 0.1875, 0.8125, 0.3125],
            [0.09375, 0.46875, 0.46875, 0.65625, 0.28125, 0.96875, 0.53125, 0.84375],
            [0.59375, 0.96875, 0.96875, 0.15625, 0.78125, 0.46875, 0.03125, 0.34375],
            [0.84375, 0.21875, 0.21875, 0.90625, 0.53125, 0.21875, 0.78125, 0.09375],
            [0.34375, 0.71875, 0.71875, 0.40625, 0.03125, 0.71875, 0.28125, 0.59375],
            [0.46875, 0.09375, 0.84375, 0.28125, 0.15625, 0.84375, 0.90625, 0.21875],
        ];
        let mut seq = SobolSeq::new(8).with_buf();
        for e in &expected {
            let x = seq.gen();
            assert_eq!(x, e.as_ref());
        }
    }

    #[test]
    fn test_sobol_seq_mean() {
        const LEN: usize = 100_000;
        const NDIM: usize = 20;
        let mut seq = SobolSeq::new(NDIM).with_buf();
        let mut sum = vec![0.; NDIM];
        for _ in 0..LEN {
            for (i, &x) in seq.gen().iter().enumerate() {
                sum[i] += x;
            }
        }
        for i in 0..NDIM {
            assert!(((sum[i] / (LEN as f64)) - 0.5).abs() < 1e-4);
        }
    }

    #[test]
    fn test_raw_data() {
        assert_eq!(get_raw_data(0), (0, vec![].as_slice()));
        assert_eq!(get_raw_data(1), (1, vec![1].as_slice()));
        assert_eq!(get_raw_data(2), (1, vec![1, 0].as_slice()));
        assert_eq!(get_raw_data(3), (2, vec![0, 0].as_slice()));
        assert_eq!(get_raw_data(4), (1, vec![0, 1, 1].as_slice()));
        assert_eq!(get_raw_data(5), (4, vec![1, 2, 6].as_slice()));
        assert_eq!(get_raw_data(6), (2, vec![0, 2, 2, 8].as_slice()));
        assert_eq!(get_raw_data(7), (4, vec![0, 2, 2, 2].as_slice()));
        assert_eq!(get_raw_data(100), (4, vec![1, 0, 7, 2, 2, 18, 113, 111, 229].as_slice()));
        assert_eq!(
            get_raw_data(1108),
            (4091, vec![1, 0, 1, 13, 1, 1, 86, 195, 106, 401, 1640, 1603].as_slice())
        );
        assert_eq!(
            get_raw_data(1109),
            (4094, vec![0, 2, 7, 9, 0, 3, 105, 78, 301, 201, 693, 791].as_slice())
        );
    }

    #[cfg(feature = "sobol-high-dim")]
    fn test_raw_data_high_dim() {
        assert_eq!(
            get_raw_data(1110),
            (21, vec![1, 2, 6, 8, 26, 62, 6, 169, 361, 260, 206, 2900, 5225].as_slice())
        );
        assert_eq!(
            get_raw_data(21198),
            (
                131020,
                vec![
                    0, 2, 0, 9, 0, 41, 1, 212, 436, 971, 1967, 2128, 7293, 5914, 27608, 10981,
                    19841
                ]
                .as_slice()
            )
        );
        assert_eq!(
            get_raw_data(21199),
            (
                131059,
                vec![
                    0, 3, 5, 7, 3, 18, 119, 168, 122, 778, 1840, 3678, 4819, 13683, 13434, 57301,
                    43158
                ]
                .as_slice()
            )
        );
    }

    #[test]
    #[should_panic]
    fn test_raw_data_panic() {
        let index = if cfg!(feature = "sobol-high-dim") { 21200 } else { 1110 };
        get_raw_data(index);
    }
}
