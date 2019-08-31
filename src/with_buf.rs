use crate::traits::QRng;

#[derive(Clone)]
pub struct QRngWithBuf<R: QRng> {
    qrng: R,
    buf: Vec<f64>,
}

impl<R: QRng> QRngWithBuf<R> {
    #[inline(always)]
    pub fn new(qrng: R) -> Self {
        let ndim = qrng.ndim();
        Self { qrng, buf: vec![0.; ndim] }
    }

    #[inline(always)]
    pub fn gen(&mut self) -> &[f64] {
        unsafe { self.qrng.gen_fill_unchecked(&mut self.buf) };
        &self.buf
    }
}

impl<R: QRng> QRng for QRngWithBuf<R> {
    #[inline(always)]
    fn ndim(&self) -> usize {
        self.qrng.ndim()
    }

    #[inline(always)]
    unsafe fn gen_fill_unchecked(&mut self, out: &mut [f64]) {
        self.qrng.gen_fill_unchecked(out);
    }
}
