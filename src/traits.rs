use crate::with_buf::QRngWithBuf;

/// Multi-dimensional quasi-random sequence generator.
///
/// Implementors of this trait are sequence generators that hold internal mutable
/// state and can generate the next multi-dimensional sample point on demand.
/// Each sample point is an `f64` slice containing values between 0 and 1.
///
/// This is the core trait implemented by multi-dimensional QRNGs.
/// For most of the use cases, it is meant to be imported in user's
/// scope since it encapsulates most of the high-level methods.
pub trait QRng: Clone {
    /// Returns the sequence dimensionality.
    ///
    /// On each step, this many coordinates of a sample point will be generated.
    fn ndim(&self) -> usize;

    /// Writes the next element of the sequence to `out` (no bounds checks).
    ///
    /// The output values are expected to be floating-point numbers between 0 and 1.
    ///
    /// Note: this method does **not** perform bound checks and is not meant to be
    /// called directly from outside of this crate. If called, it is the user's
    /// responsibility to provide a buffer of length `ndim()` or higher.
    unsafe fn gen_fill_unchecked(&mut self, out: &mut [f64]);

    /// Writes the next element of the sequence to `out` (with a bounds check).
    ///
    /// The output values are floating-point numbers between 0 and 1.
    ///
    /// # Examples
    ///
    /// ```
    /// # use qrng::*;
    /// let mut buf = vec![0.; 3];
    /// let mut seq = HaltonSeq::new(3);
    /// seq.gen_fill(&mut buf);
    /// ```
    #[inline]
    fn gen_fill(&mut self, out: &mut [f64]) {
        if out.len() < self.ndim() {
            panic!(
                "index out of bounds: the len is {} but the index is {}",
                out.len(),
                self.ndim()
            );
        }
        unsafe {
            self.gen_fill_unchecked(out);
        }
    }

    /// Returns a wrapper (TODO: ...).
    ///
    /// See [`QRngWithBuf`](struct.QRngWithBuf.html) (TODO: ...).
    ///
    /// # Examples
    ///
    /// ```
    /// # use qrng::*;
    /// let mut seq = SobolSeq::new(5).with_buf();
    /// let next = seq.gen();
    /// ```
    #[inline]
    fn with_buf(self) -> QRngWithBuf<Self> {
        QRngWithBuf::new(self)
    }
}
