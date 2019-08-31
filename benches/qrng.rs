#[macro_use]
extern crate criterion;

extern crate qrng;

use qrng::{HaltonSeq, QRng, SobolSeq};

use criterion::{black_box, Criterion};

const NDIM_SET: &[usize] = &[2, 20, 100, 1000];
const NDIM_LEN_SET: &[(usize, usize)] = &[(2, 20_000), (20, 5_000), (100, 1_000), (1000, 100)];

fn halton_seq_gen(c: &mut Criterion) {
    for &(ndim, len) in NDIM_LEN_SET {
        c.bench_function(&format!("HaltonSeq::gen (ndim={}, len={})", ndim, len), |b| {
            let seq = HaltonSeq::new(black_box(ndim)).with_buf();
            b.iter(|| {
                let mut seq = seq.clone();
                for _ in 0..black_box(len) {
                    seq.gen();
                }
            })
        });
    }
}

fn sobol_seq_new(c: &mut Criterion) {
    for &ndim in NDIM_SET {
        c.bench_function(&format!("SobolSeq::new (ndim={})", ndim), |b| {
            b.iter(|| SobolSeq::new(black_box(ndim)))
        });
    }
}

fn sobol_seq_gen(c: &mut Criterion) {
    for &(ndim, len) in NDIM_LEN_SET {
        c.bench_function(&format!("SobolSeq::gen (ndim={}, len={})", ndim, len), |b| {
            let seq = SobolSeq::new(black_box(ndim)).with_buf();
            b.iter(|| {
                let mut seq = seq.clone();
                for _ in 0..black_box(len) {
                    seq.gen();
                }
            })
        });
    }
}

criterion_group!(benches, halton_seq_gen, sobol_seq_new, sobol_seq_gen);
criterion_main!(benches);
