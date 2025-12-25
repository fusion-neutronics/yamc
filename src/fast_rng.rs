// Fast random number generator matching OpenMC's PCG-LCG implementation
// See: openmc-source-code/random_lcg.cpp
//
// This is significantly faster than rand_pcg::Pcg64 because:
// 1. No struct construction overhead - just a u64 seed
// 2. Fully inlineable - compiler can optimize completely
// 3. Minimal state - 8 bytes vs Pcg64's larger internal state

use rand::{RngCore, SeedableRng};

/// LCG multiplier (same as OpenMC)
const PRN_MULT: u64 = 6364136223846793005;
/// LCG additive constant (same as OpenMC)
const PRN_ADD: u64 = 1442695040888963407;

/// Fast RNG using OpenMC's PCG-LCG algorithm.
///
/// This is a PCG (Permuted Congruential Generator) variant that uses
/// an LCG as the base generator with output permutation for quality.
///
/// Reference: Melissa E. O'Neill, "PCG: A Family of Simple Fast Space-Efficient
/// Statistically Good Algorithms for Random Number Generation"
#[derive(Clone, Copy, Debug)]
pub struct FastRng {
    seed: u64,
}

impl FastRng {
    /// Create a new FastRng with the given seed
    #[inline]
    pub fn new(seed: u64) -> Self {
        Self { seed }
    }

    /// Generate a random f64 in [0, 1) - matches OpenMC's prn() function
    #[inline(always)]
    pub fn random(&mut self) -> f64 {
        // Advance the LCG
        self.seed = PRN_MULT.wrapping_mul(self.seed).wrapping_add(PRN_ADD);

        // PCG output permutation (RXS-M-XS variant)
        let word = ((self.seed >> ((self.seed >> 59) + 5)) ^ self.seed)
            .wrapping_mul(12605985483714917081);
        let result = (word >> 43) ^ word;

        // Convert to f64 in [0, 1) - equivalent to ldexp(result, -64)
        (result as f64) * 5.421010862427522e-20
    }

    /// Reseed the RNG (for reuse across particles)
    #[inline]
    pub fn reseed(&mut self, seed: u64) {
        self.seed = seed;
    }
}

impl SeedableRng for FastRng {
    type Seed = [u8; 8];

    fn from_seed(seed: Self::Seed) -> Self {
        Self {
            seed: u64::from_le_bytes(seed),
        }
    }
}

impl RngCore for FastRng {
    #[inline(always)]
    fn next_u32(&mut self) -> u32 {
        self.next_u64() as u32
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        // Advance the LCG
        self.seed = PRN_MULT.wrapping_mul(self.seed).wrapping_add(PRN_ADD);

        // PCG output permutation
        let word = ((self.seed >> ((self.seed >> 59) + 5)) ^ self.seed)
            .wrapping_mul(12605985483714917081);
        (word >> 43) ^ word
    }

    #[inline]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        // Fill bytes using next_u64
        let mut left = dest;
        while left.len() >= 8 {
            let bytes = self.next_u64().to_le_bytes();
            left[..8].copy_from_slice(&bytes);
            left = &mut left[8..];
        }
        if !left.is_empty() {
            let bytes = self.next_u64().to_le_bytes();
            left.copy_from_slice(&bytes[..left.len()]);
        }
    }

    #[inline]
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_fast_rng_deterministic() {
        let mut rng1 = FastRng::new(12345);
        let mut rng2 = FastRng::new(12345);

        for _ in 0..100 {
            assert_eq!(rng1.random(), rng2.random());
        }
    }

    #[test]
    fn test_fast_rng_range() {
        let mut rng = FastRng::new(42);

        for _ in 0..10000 {
            let val = rng.random();
            assert!(val >= 0.0 && val < 1.0, "Value {} out of range [0, 1)", val);
        }
    }

    #[test]
    fn test_fast_rng_as_rand_rng() {
        // Test that FastRng works with rand's Rng trait
        let mut rng = FastRng::new(12345);

        let _: f64 = rng.gen();
        let _: u32 = rng.gen();
        let _: bool = rng.gen();
    }

    #[test]
    fn test_fast_rng_reseed() {
        let mut rng = FastRng::new(12345);
        let first_val = rng.random();

        // Generate more values
        for _ in 0..100 {
            rng.random();
        }

        // Reseed and verify we get the same sequence
        rng.reseed(12345);
        assert_eq!(rng.random(), first_val);
    }
}
