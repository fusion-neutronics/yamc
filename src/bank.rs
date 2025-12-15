// Particle banking system
// Corresponds to OpenMC's bank.cpp
//
// Handles secondary particle queuing from multi-neutron reactions (n,2n), (n,3n), etc.
// Future: Will also handle fission bank for eigenvalue/criticality calculations

use crate::particle::Particle;
use std::collections::VecDeque;

/// Local particle bank/queue for storing secondary particles during transport
/// 
/// In OpenMC terminology:
/// - source_bank: Primary source particles for a batch
/// - fission_bank: Fission sites for eigenvalue calculations (not yet implemented)
/// - Secondary particles: Stored in thread-local queues during transport
///
/// This struct manages the secondary particle queue for a single history.
pub struct ParticleBank {
    /// Queue of particles to be transported (primary + secondaries)
    queue: VecDeque<Particle>,
}

impl ParticleBank {
    /// Create a new empty particle bank
    pub fn new() -> Self {
        ParticleBank {
            queue: VecDeque::new(),
        }
    }

    /// Create a particle bank with an initial capacity
    pub fn with_capacity(capacity: usize) -> Self {
        ParticleBank {
            queue: VecDeque::with_capacity(capacity),
        }
    }

    /// Add a primary particle to the bank
    pub fn add_source_particle(&mut self, particle: Particle) {
        self.queue.push_back(particle);
    }

    /// Bank a secondary particle from a reaction (e.g., from n,2n or n,3n)
    pub fn bank_secondary(&mut self, particle: Particle) {
        self.queue.push_back(particle);
    }

    /// Get the next particle from the bank for transport
    /// Returns None if the bank is empty
    pub fn pop_particle(&mut self) -> Option<Particle> {
        self.queue.pop_front()
    }

    /// Check if the bank is empty
    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    /// Get the number of particles in the bank
    pub fn len(&self) -> usize {
        self.queue.len()
    }

    /// Clear all particles from the bank
    pub fn clear(&mut self) {
        self.queue.clear();
    }

    /// Reserve capacity for additional particles
    pub fn reserve(&mut self, additional: usize) {
        self.queue.reserve(additional);
    }
}

impl Default for ParticleBank {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_particle_bank_basic() {
        let mut bank = ParticleBank::new();
        assert!(bank.is_empty());
        assert_eq!(bank.len(), 0);

        // Create a test particle
        let particle = Particle {
            energy: 1.0e6,
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            alive: true,
            id: 1,
            current_cell_index: None,
        };

        bank.add_source_particle(particle.clone());
        assert_eq!(bank.len(), 1);
        assert!(!bank.is_empty());

        let retrieved = bank.pop_particle().unwrap();
        assert_eq!(retrieved.energy, 1.0e6);
        assert!(bank.is_empty());
    }

    #[test]
    fn test_particle_bank_secondary() {
        let mut bank = ParticleBank::new();

        let primary = Particle {
            energy: 14.0e6,
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            alive: true,
            id: 1,
            current_cell_index: None,
        };

        let secondary = Particle {
            energy: 7.0e6,
            position: [0.0, 0.0, 0.0],
            direction: [1.0, 0.0, 0.0],
            alive: true,
            id: 2,
            current_cell_index: None,
        };

        bank.add_source_particle(primary);
        bank.bank_secondary(secondary);

        assert_eq!(bank.len(), 2);

        // Primary should come out first (FIFO)
        let p1 = bank.pop_particle().unwrap();
        assert_eq!(p1.energy, 14.0e6);

        // Then secondary
        let p2 = bank.pop_particle().unwrap();
        assert_eq!(p2.energy, 7.0e6);

        assert!(bank.is_empty());
    }

    #[test]
    fn test_particle_bank_with_capacity() {
        let bank = ParticleBank::with_capacity(100);
        assert!(bank.is_empty());
        assert!(bank.queue.capacity() >= 100);
    }

    #[test]
    fn test_particle_bank_clear() {
        let mut bank = ParticleBank::new();
        
        let particle = Particle {
            energy: 1.0e6,
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            alive: true,
            id: 1,
            current_cell_index: None,
        };

        bank.add_source_particle(particle.clone());
        bank.add_source_particle(particle.clone());
        bank.add_source_particle(particle);

        assert_eq!(bank.len(), 3);
        
        bank.clear();
        assert!(bank.is_empty());
        assert_eq!(bank.len(), 0);
    }
}
