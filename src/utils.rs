#[derive(Clone)]
pub struct PrimeSeq {
    primes: Vec<u64>,
    index: usize,
}

impl PrimeSeq {
    #[inline]
    pub fn new() -> Self {
        Self { primes: vec![2, 3], index: 0 }
    }
}

impl Iterator for PrimeSeq {
    type Item = u64;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.primes.len() <= self.index {
            let mut x = self.primes[self.primes.len() - 1];
            'outer: loop {
                x += 2;
                for &p in &self.primes {
                    if x % p == 0 {
                        continue 'outer;
                    } else if p * p > x {
                        break;
                    }
                }
                self.primes.push(x);
                break;
            }
        }
        self.index += 1;
        Some(self.primes[self.index - 1])
    }
}

pub fn primes() -> PrimeSeq {
    PrimeSeq::new()
}

#[cfg(test)]
mod tests {
    use super::primes;

    #[test]
    fn test_prime_seq() {
        assert_eq!(primes().take(10).collect::<Vec<_>>(), vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }
}
