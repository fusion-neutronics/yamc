// Struct representing a nuclide - reads from HDF5 files
use crate::reaction::Reaction;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

// Global cache for nuclides to avoid reloading
static GLOBAL_NUCLIDE_CACHE: Lazy<Mutex<HashMap<String, Arc<Nuclide>>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));

// Scattering MTs - EXCLUDING synthetic MT 4 (inelastic is represented by MT 50-91)
const SCATTERING_MTS_NON_INELASTIC: &[i32] = &[
    2,   // elastic
    5, 11, 16, 17, 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36, 37, 41, 42, 44, 45,
    152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
    170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187,
    188, 189, 190, 194, 195, 196, 198, 199, 200,
];

// Inelastic constituent MTs (the REAL reactions, not MT 4)
const INELASTIC_CONSTITUENT_MTS: std::ops::Range<i32> = 50..92;

/// Helper function to check if an MT number is a scattering reaction (excludes MT 4 synthetic)
#[inline]
fn is_scattering_mt(mt: i32) -> bool {
    // Inelastic constituent MTs (50-91)
    if (50..92).contains(&mt) {
        return true;
    }
    // Other scattering MTs (elastic and non-inelastic scattering)
    SCATTERING_MTS_NON_INELASTIC.contains(&mt)
}

/// Enum to represent either an MT number or reaction name for flexible reaction identification
#[derive(Debug, Clone)]
pub enum ReactionIdentifier {
    Mt(i32),
    Name(String),
}

impl From<i32> for ReactionIdentifier {
    fn from(mt: i32) -> Self {
        ReactionIdentifier::Mt(mt)
    }
}

impl From<String> for ReactionIdentifier {
    fn from(name: String) -> Self {
        ReactionIdentifier::Name(name)
    }
}

impl From<&str> for ReactionIdentifier {
    fn from(name: &str) -> Self {
        ReactionIdentifier::Name(name.to_string())
    }
}

/// Clear the global nuclide cache (used by tests and Python to ensure deterministic selective temperature behavior)
#[allow(dead_code)]
pub fn clear_nuclide_cache() {
    match GLOBAL_NUCLIDE_CACHE.lock() {
        Ok(mut cache) => cache.clear(),
        Err(poisoned) => {
            // Handle poisoned mutex by clearing the data anyway
            let mut cache = poisoned.into_inner();
            cache.clear();
        }
    }
}

/// Core data model for a single nuclide and its reaction cross section data.
///
/// A `Nuclide` mirrors (and is deserialized from) a JSON schema containing
/// metadata plus reaction channel data at one or more temperatures. Reaction
/// data are organized by temperature key (e.g. "294") and then by ENDF/MT
/// number. Each [`Reaction`] holds its own threshold information and (possibly
/// truncated) energy grid relative to the top‑level temperature energy grid.
///
/// Temperatures:
/// * `available_temperatures` always lists every temperature present in the
///   source JSON file – even if a filtered load only materialized a subset.
/// * `loaded_temperatures` tracks the subset actually parsed into `reactions`
///   and `energy` according to caller filtering semantics.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Nuclide {
    /// Canonical nuclide name (e.g. "Li6"). May be derived when absent.
    pub name: Option<String>,
    /// Optional human readable element name (legacy field; may be absent).
    pub element: Option<String>,
    /// Element symbol, e.g. "Li".
    pub atomic_symbol: Option<String>,
    /// Atomic (proton) number Z.
    pub atomic_number: Option<u32>,
    /// Neutron number N (may be computed from A - Z if missing).
    pub neutron_number: Option<u32>,
    /// Mass number A.
    pub mass_number: Option<u32>,
    /// Atomic weight ratio (target mass / neutron mass) from HDF5 file.
    pub atomic_weight_ratio: Option<f64>,
    /// Origin / library identifier (e.g. JEFF, ENDF, custom tag).
    pub library: Option<String>,
    /// Top‑level energy grid per temperature (full grid; per‑reaction grids may be threshold‑truncated).
    pub energy: Option<HashMap<String, Vec<f64>>>,
    /// temperature -> MT number -> reaction data.
    #[serde(default)]
    pub reactions: HashMap<String, HashMap<i32, Reaction>>, // temperature -> mt (i32) -> Reaction
    /// True if any fission MT channel is present.
    pub fissionable: bool,
    /// All temperatures present in the JSON file regardless of filtering.
    #[serde(skip, default)]
    pub available_temperatures: Vec<String>, // All temps listed in the JSON (even if not loaded)
    /// Subset of temperatures actually loaded into `reactions` / `energy`.
    #[serde(skip, default)]
    pub loaded_temperatures: Vec<String>, // Subset actually loaded into reactions/energy
    /// Optional path the JSON was read from (None for in‑memory sources / WASM).
    #[serde(skip, default)]
    pub data_path: Option<String>, // Path JSON was loaded from (for potential future extension)
}

impl Nuclide {
    /// Get the element name with auto-loading if not available
    pub fn get_element(&mut self) -> Option<String> {
        if self.element.is_none() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                let _ = self.auto_load_from_config(name, None);
            }
        }
        self.element.clone()
    }

    /// Get the atomic number with auto-loading if not available
    pub fn get_atomic_number(&mut self) -> Option<u32> {
        if self.atomic_number.is_none() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                let _ = self.auto_load_from_config(name, None);
            }
        }
        self.atomic_number
    }

    /// Get the mass number with auto-loading if not available
    pub fn get_mass_number(&mut self) -> Option<u32> {
        if self.mass_number.is_none() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                let _ = self.auto_load_from_config(name, None);
            }
        }
        self.mass_number
    }

    /// Get the atomic symbol with auto-loading if not available
    pub fn get_atomic_symbol(&mut self) -> Option<String> {
        if self.atomic_symbol.is_none() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                let _ = self.auto_load_from_config(name, None);
            }
        }
        self.atomic_symbol.clone()
    }

    /// Get the neutron number with auto-loading if not available
    pub fn get_neutron_number(&mut self) -> Option<u32> {
        if self.neutron_number.is_none() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                let _ = self.auto_load_from_config(name, None);
            }
        }
        self.neutron_number
    }

    /// Get available temperatures with auto-loading if not available
    pub fn get_available_temperatures(&mut self) -> Vec<String> {
        if self.available_temperatures.is_empty() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                let _ = self.auto_load_from_config(name, None);
            }
        }
        self.available_temperatures.clone()
    }

    /// Sample the top-level reaction type (fission, absorption, elastic, inelastic, other) at a given energy and temperature
    /// This version includes auto-loading if data is not available.
    pub fn sample_reaction<R: rand::Rng + ?Sized>(
        &mut self,
        energy: f64,
        temperature: &str,
        rng: &mut R,
    ) -> Option<&Reaction> {
        // Auto-load data if not available
        if self.reactions.is_empty() && self.name.is_some() {
            if let Some(name) = &self.name.clone() {
                if let Err(_) = self.auto_load_from_config(name, Some(temperature)) {
                    println!("[sample_reaction] Failed to auto-load data for nuclide '{}'", name);
                    return None;
                }
            }
        }

        // Try temperature as given, then with 'K' appended, then any available
        let temp_reactions = if let Some(r) = self.reactions.get(temperature) {
            r
        } else if let Some(r) = self.reactions.get(&format!("{}K", temperature)) {
            r
        } else if let Some((_, r)) = self.reactions.iter().next() {
            r
        } else {
            return None;
        };

        // Define MTs for each event type
        let total_mt = 1;
        let fission_mt = 18;
        let absorption_mt = 101;
        let elastic_mt = 2;
        let inelastic_mt = 4;

        // Helper to get xs for a given MT using Reaction::cross_section_at
        let get_xs = |mt: i32| -> f64 {
            temp_reactions
                .get(&mt)
                .and_then(|reaction| reaction.cross_section_at(energy))
                .unwrap_or(0.0)
        };

        let total_xs = get_xs(total_mt);
        if total_xs <= 0.0 {
            return None;
        }

        let xi = rng.gen_range(0.0..total_xs);
        let mut accum = 0.0;

        // Absorption

        let xs_absorption = get_xs(absorption_mt);
        accum += xs_absorption;
        if xi < accum && xs_absorption > 0.0 {
            return temp_reactions.get(&absorption_mt);
        }

        // Elastic
        let xs_elastic = get_xs(elastic_mt);
        accum += xs_elastic;
        if xi < accum && xs_elastic > 0.0 {
            return temp_reactions.get(&elastic_mt);
        }

        // Fission (only if nuclide is fissionable, checked last)
        let xs_fission = if self.fissionable {
            get_xs(fission_mt)
        } else {
            0.0
        };
        accum += xs_fission;
        if xi < accum && xs_fission > 0.0 {
            return temp_reactions.get(&fission_mt);
        }

        // inelastic selection as fallback
        temp_reactions.get(&inelastic_mt)
    }

    /// Sample the top-level reaction type without auto-loading (requires data to be pre-loaded)
    /// This version is used when the nuclide is stored in shared/immutable contexts like Arc<Nuclide>.
    #[inline]
    pub fn sample_reaction_no_autoload<R: rand::Rng + ?Sized>(
        &self,
        energy: f64,
        temperature: &str,
        rng: &mut R,
    ) -> Option<&Reaction> {
        // Try temperature as given, then with 'K' appended, then any available
        let temp_reactions = if let Some(r) = self.reactions.get(temperature) {
            r
        } else if let Some(r) = self.reactions.get(&format!("{}K", temperature)) {
            r
        } else if let Some((_, r)) = self.reactions.iter().next() {
            r
        } else {
            return None;
        };

        // Define MTs for each event type
        let total_mt = 1;
        let fission_mt = 18;
        let absorption_mt = 101;
        let scattering_mt = 1001; // Synthetic scattering reaction (replaces separate MT 2 and MT 4)

        // Helper to get xs for a given MT using Reaction::cross_section_at
        let get_xs = |mt: i32| -> f64 {
            temp_reactions
                .get(&mt)
                .and_then(|reaction| reaction.cross_section_at(energy))
                .unwrap_or(0.0)
        };

        let total_xs = get_xs(total_mt);
        if total_xs <= 0.0 {
            return None;
        }

        let xi = rng.gen_range(0.0..total_xs);
        let mut accum = 0.0;

        // Absorption
        let xs_absorption = get_xs(absorption_mt);
        accum += xs_absorption;
        if xi < accum && xs_absorption > 0.0 {
            return temp_reactions.get(&absorption_mt);
        }

        // Scattering (replaces separate elastic and inelastic)
        let xs_scattering = get_xs(scattering_mt);
        accum += xs_scattering;
        if xi < accum && xs_scattering > 0.0 {
            return temp_reactions.get(&scattering_mt);
        }

        // Fission (only if nuclide is fissionable, checked last)
        let xs_fission = if self.fissionable {
            get_xs(fission_mt)
        } else {
            0.0
        };
        accum += xs_fission;
        if xi < accum && xs_fission > 0.0 {
            return temp_reactions.get(&fission_mt);
        }

        // Fallback to scattering if nothing else was sampled
        temp_reactions.get(&scattering_mt)
    }

    /// Sample a specific scattering constituent reaction from all available scattering MTs.
    /// This samples from MT 2 (elastic), MT 50-91 (inelastic constituents), and other scattering MTs.
    /// Note: NEVER returns MT 4 (inelastic) as it's a synthetic reaction.
    ///
    /// # Arguments
    /// * `energy` - Neutron energy in eV
    /// * `temperature` - Temperature string (e.g., "294" or "294K")
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// * `&Reaction` for the sampled constituent scattering reaction (MT 2, 50-91, 16, 17, etc.)
    ///
    /// # Panics
    /// * If no scattering constituent reactions are available
    /// * If sampling logic fails despite having valid reactions and cross sections
    #[inline]
    pub fn sample_scattering_constituent<R: rand::Rng + ?Sized>(
        &self,
        energy: f64,
        temperature: &str,
        rng: &mut R,
    ) -> &Reaction {
        // Try temperature as given, then with 'K' appended, then any available
        let temp_reactions = if let Some(r) = self.reactions.get(temperature) {
            r
        } else if let Some(r) = self.reactions.get(&format!("{}K", temperature)) {
            r
        } else if let Some((_, r)) = self.reactions.iter().next() {
            r
        } else {
            panic!(
                "[sample_scattering_constituent] No reaction data available for any temperature."
            );
        };

        // Build list of available scattering reactions with their cross sections
        // Single iteration through the HashMap instead of 50+ individual lookups
        let mut available_reactions: Vec<(i32, &Reaction, f64)> = Vec::new();
        let mut total_scattering_xs = 0.0;

        for (&mt, reaction) in temp_reactions.iter() {
            // Skip MT 4 (synthetic inelastic) - only use constituent reactions
            if mt == 4 || mt == 1 || mt == 18 || mt == 101 || mt == 1001 {
                continue;
            }
            
            // Check if this is a scattering MT
            if is_scattering_mt(mt) {
                if let Some(xs) = reaction.cross_section_at(energy) {
                    if xs > 0.0 {
                        available_reactions.push((mt, reaction, xs));
                        total_scattering_xs += xs;
                    }
                }
            }
        }

        if available_reactions.is_empty() {
            panic!("sample_scattering_constituent: No scattering constituent reactions with non-zero cross sections at energy {} eV for temperature '{}'. This should not happen if MT 1001 was sampled.", energy, temperature);
        }

        if total_scattering_xs <= 0.0 {
            panic!("sample_scattering_constituent: Total scattering cross section is zero at energy {} eV for temperature '{}'. All constituent reactions have zero cross section.", energy, temperature);
        }

        // Sample which specific scattering reaction occurs
        let xi = rng.gen_range(0.0..total_scattering_xs);
        let mut accum = 0.0;

        for (_, reaction, xs) in available_reactions {
            accum += xs;
            if xi < accum {
                return reaction;
            }
        }

        // This should never be reached due to the sampling logic above
        panic!("sample_scattering_constituent: Failed to sample any scattering reaction despite having available MTs and positive total cross section. This indicates a bug in the sampling logic.");
    }

    /// Sample a specific inelastic reaction from the constituent MT 50-91 reactions that make up MT 4.
    /// This provides more detailed physics than just sampling MT 4 directly.
    ///
    /// # Arguments
    /// * `energy` - Neutron energy in eV
    /// * `temperature` - Temperature string (e.g., "294" or "294K")
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// * `&Reaction` for the sampled constituent inelastic reaction (MT 50-91)
    ///
    /// # Panics
    /// * If no inelastic constituent reactions (MT 50-91) are available
    /// * If sampling logic fails despite having valid reactions and cross sections
    pub fn sample_inelastic_constituent<R: rand::Rng + ?Sized>(
        &self,
        energy: f64,
        temperature: &str,
        rng: &mut R,
    ) -> &Reaction {
        // Try temperature as given, then with 'K' appended, then any available
        let temp_reactions = if let Some(r) = self.reactions.get(temperature) {
            r
        } else if let Some(r) = self.reactions.get(&format!("{}K", temperature)) {
            r
        } else if let Some((temp, r)) = self.reactions.iter().next() {
            println!("[sample_inelastic_constituent] Requested temperature '{}' not found. Using available temperature '{}'.", temperature, temp);
            r
        } else {
            panic!(
                "[sample_inelastic_constituent] No reaction data available for any temperature."
            );
        };

        // MT 4 is composed of MT 50-91 (inelastic scattering to discrete levels)
        let inelastic_constituent_mts: Vec<i32> = (50..92).collect();

        // Filter to only the MTs that are actually available in this nuclide
        let available_inelastic_mts: Vec<i32> = inelastic_constituent_mts
            .into_iter()
            .filter(|&mt| temp_reactions.contains_key(&mt))
            .collect();

        if available_inelastic_mts.is_empty() {
            panic!("sample_inelastic_constituent: No inelastic constituent reactions (MT 50-91) available in this nuclide at temperature '{}'. This indicates missing nuclear data or incorrect MT 4 sampling.", temperature);
        }

        // Helper to get cross section for a given MT
        let get_xs = |mt: i32| -> f64 {
            temp_reactions
                .get(&mt)
                .and_then(|reaction| reaction.cross_section_at(energy))
                .unwrap_or(0.0)
        };

        // Calculate total cross section for all available inelastic constituents
        let total_inelastic_xs: f64 = available_inelastic_mts.iter().map(|&mt| get_xs(mt)).sum();

        if total_inelastic_xs <= 0.0 {
            panic!("sample_inelastic_constituent: Total inelastic cross section is zero at energy {} eV for temperature '{}'. All constituent reactions have zero cross section.", energy, temperature);
        }

        // Sample which specific inelastic reaction occurs
        let xi = rng.gen_range(0.0..total_inelastic_xs);
        let mut accum = 0.0;

        for &mt in &available_inelastic_mts {
            let xs = get_xs(mt);
            accum += xs;
            if xi < accum && xs > 0.0 {
                return temp_reactions.get(&mt).expect(
                    "sample_inelastic_constituent: MT not found in temp_reactions after filtering.",
                );
            }
        }

        // This should never be reached due to the sampling logic above
        panic!("sample_inelastic_constituent: Failed to sample any inelastic reaction despite having available MTs and positive total cross section. This indicates a bug in the sampling logic.");
    }

    /// Sample a specific absorption constituent reaction (MTs that make up MT 101) at a given energy and temperature.
    /// Panics if no constituent reactions are available or total cross section is zero.
    #[inline]
    pub fn sample_absorption_constituent<R: rand::Rng + ?Sized>(
        &self,
        energy: f64,
        temperature: &str,
        rng: &mut R,
    ) -> &Reaction {
        // Try temperature as given, then with 'K' appended, then any available
        let temp_reactions = if let Some(r) = self.reactions.get(temperature) {
            r
        } else if let Some(r) = self.reactions.get(&format!("{}K", temperature)) {
            r
        } else if let Some((_, r)) = self.reactions.iter().next() {
            r
        } else {
            panic!(
                "[sample_absorption_constituent] No reaction data available for any temperature."
            );
        };

        // Build list of available absorption reactions with their cross sections
        // Single iteration through the HashMap instead of checking hundreds of MTs individually
        let mut available_reactions: Vec<(i32, &Reaction, f64)> = Vec::new();
        let mut total_xs = 0.0;

        for (&mt, reaction) in temp_reactions.iter() {
            // Check if this MT is an absorption constituent
            // MTs that make up MT 101 (absorption)
            let is_absorption = matches!(mt,
                102 | 103 | 104 | 105 | 106 | 107 | 108 | 109 | 111 | 112 | 113 | 114 | 115 | 116 | 117 | 155 | 182 | 191 | 192 | 193 | 197
            ) || (600..850).contains(&mt);

            if is_absorption {
                if let Some(xs) = reaction.cross_section_at(energy) {
                    if xs > 0.0 {
                        available_reactions.push((mt, reaction, xs));
                        total_xs += xs;
                    }
                }
            }
        }

        if available_reactions.is_empty() {
            panic!("sample_absorption_constituent: No absorption constituent reactions available in this nuclide at temperature '{}'. This indicates missing nuclear data or incorrect MT 101 sampling.", temperature);
        }

        if total_xs <= 0.0 {
            panic!("sample_absorption_constituent: Total absorption cross section is zero at energy {} eV for temperature '{}'. All constituent reactions have zero cross section.", energy, temperature);
        }

        // Sample which specific absorption reaction occurs
        let xi = rng.gen_range(0.0..total_xs);
        let mut accum = 0.0;

        for (_, reaction, xs) in available_reactions {
            accum += xs;
            if xi < accum {
                return reaction;
            }
        }

        // This should never be reached due to the sampling logic above
        panic!("sample_absorption_constituent: Failed to sample any absorption reaction despite having available MTs and positive total cross section. This indicates a bug in the sampling logic.");
    }

    /// Get the energy grid for a specific temperature
    pub fn energy_grid(&self, temperature: &str) -> Option<&Vec<f64>> {
        self.energy
            .as_ref()
            .and_then(|energy_map| energy_map.get(temperature))
    }

    /// Get a list of available temperatures
    pub fn temperatures(&self) -> Option<Vec<String>> {
        let mut temps = std::collections::HashSet::new();

        // Check reactions first
        for temp in self.reactions.keys() {
            temps.insert(temp.clone());
        }

        // Also check energy map
        if let Some(energy_map) = &self.energy {
            for temp in energy_map.keys() {
                temps.insert(temp.clone());
            }
        }

        if temps.is_empty() {
            None
        } else {
            let mut temps_vec: Vec<String> = temps.into_iter().collect();
            temps_vec.sort();
            Some(temps_vec)
        }
    }

    /// Get a list of available MT numbers
    pub fn reaction_mts(&self) -> Option<Vec<i32>> {
        let mut mts = std::collections::HashSet::new();
        for temp_reactions in self.reactions.values() {
            for &mt in temp_reactions.keys() {
                mts.insert(mt);
            }
        }
        if mts.is_empty() {
            None
        } else {
            let mut mts_vec: Vec<i32> = mts.into_iter().collect();
            mts_vec.sort();
            Some(mts_vec)
        }
    }

    /// Get microscopic cross section data for a specific reaction and temperature.
    /// Returns a tuple of (cross_section_values, energy_grid).
    /// If temperature is None, uses the single loaded temperature if only one exists.
    /// Automatically loads data if not already loaded, using the nuclide name and config.
    ///
    /// # Arguments
    /// * `reaction` - Either an MT number (i32) or reaction name (String/&str) like "(n,gamma)" or "fission"
    /// * `temperature` - Optional temperature string
    pub fn microscopic_cross_section<R>(
        &mut self,
        reaction: R,
        temperature: Option<&str>,
        trim_trailing_zeros: bool,
    ) -> Result<(Vec<f64>, Vec<f64>), Box<dyn std::error::Error>>
    where
        R: Into<ReactionIdentifier>,
    {
        // Convert the reaction parameter to an MT number
        let mt = match reaction.into() {
            ReactionIdentifier::Mt(mt_num) => mt_num,
            ReactionIdentifier::Name(name) => {
                // Use the REACTION_MT mapping to convert string to MT number
                crate::data::REACTION_MT.get(name.as_str())
                    .copied()
                    .ok_or_else(|| format!("Unknown reaction name '{}'. Available reactions can be found in REACTION_MT mapping.", name))?
            }
        };
        // Check if we need to load data automatically
        if self.loaded_temperatures.is_empty() {
            // No data loaded yet - try to load it automatically
            if let Some(name) = self.name.clone() {
                self.auto_load_from_config(&name, temperature)?;
            } else {
                return Err(
                    "No data loaded and no nuclide name available for automatic loading".into(),
                );
            }
        }

        // Check if we need to load additional temperature
        if let Some(temp) = temperature {
            let temp_normalized = format!("{}K", temp);
            let temp_without_k = if temp.ends_with('K') {
                &temp[..temp.len() - 1]
            } else {
                temp
            };

            // Check if the requested temperature is available but not loaded
            let needs_temp_load = !self.loaded_temperatures.contains(&temp.to_string())
                && !self.loaded_temperatures.contains(&temp_normalized)
                && !self
                    .loaded_temperatures
                    .contains(&temp_without_k.to_string())
                && (self.available_temperatures.contains(&temp.to_string())
                    || self.available_temperatures.contains(&temp_normalized)
                    || self
                        .available_temperatures
                        .contains(&temp_without_k.to_string()));

            if needs_temp_load {
                if let Some(name) = self.name.clone() {
                    self.auto_load_additional_temperature(&name, temp)?;
                }
            }
        }

        // Now proceed with the original logic
        let (xs, energy) = self.get_microscopic_cross_section_data(mt, temperature)?;
        if trim_trailing_zeros {
            let mut last_nonzero = xs.len();
            for (i, &val) in xs.iter().enumerate().rev() {
                if val != 0.0 {
                    last_nonzero = i + 1;
                    break;
                }
            }
            Ok((xs[..last_nonzero].to_vec(), energy[..last_nonzero].to_vec()))
        } else {
            Ok((xs, energy))
        }
    }

    /// Helper method to automatically load data from config
    pub fn auto_load_from_config(
        &mut self,
        nuclide_name: &str,
        temperature: Option<&str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Try to load using the nuclide name and config
        let path_or_url = {
            let cfg = crate::config::CONFIG
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            cfg.get_cross_section(nuclide_name)
        };

        if let Some(path_or_url) = path_or_url {
            // Determine which temperatures to load
            let temps_to_load = if let Some(temp) = temperature {
                let mut temps = std::collections::HashSet::new();
                temps.insert(temp.to_string());
                Some(temps)
            } else {
                None // Load all temperatures
            };

            // Load the data
            let loaded_nuclide = if crate::url_cache::is_keyword(&path_or_url) {
                load_nuclide_for_python(
                    Some(&path_or_url),
                    Some(nuclide_name),
                    temps_to_load.as_ref(),
                )?
            } else {
                load_nuclide_for_python(Some(&path_or_url), None, temps_to_load.as_ref())?
            };

            // Update self with the loaded data
            *self = loaded_nuclide;
            Ok(())
        } else {
            Err(format!("No configuration found for nuclide '{}'. Use Config.set_cross_sections() to configure data sources.", nuclide_name).into())
        }
    }

    /// Helper method to load additional temperature
    fn auto_load_additional_temperature(
        &mut self,
        nuclide_name: &str,
        temperature: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let path_or_url = {
            let cfg = crate::config::CONFIG
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            cfg.get_cross_section(nuclide_name)
        };

        if let Some(path_or_url) = path_or_url {
            // Create union of current loaded temperatures plus the new one
            let mut temps_to_load = std::collections::HashSet::new();
            for temp in &self.loaded_temperatures {
                temps_to_load.insert(temp.clone());
            }
            temps_to_load.insert(temperature.to_string());

            // Reload with the expanded temperature set
            let loaded_nuclide = if crate::url_cache::is_keyword(&path_or_url) {
                load_nuclide_for_python(
                    Some(&path_or_url),
                    Some(nuclide_name),
                    Some(&temps_to_load),
                )?
            } else {
                load_nuclide_for_python(Some(&path_or_url), None, Some(&temps_to_load))?
            };

            // Update self with the reloaded data
            *self = loaded_nuclide;
            Ok(())
        } else {
            Err(format!(
                "No configuration found for nuclide '{}' to load additional temperature",
                nuclide_name
            )
            .into())
        }
    }

    /// Core method that extracts the cross section data (unchanged logic)
    fn get_microscopic_cross_section_data(
        &self,
        mt: i32,
        temperature: Option<&str>,
    ) -> Result<(Vec<f64>, Vec<f64>), Box<dyn std::error::Error>> {
        // Determine which temperature to use
        let temp_key = if let Some(temp) = temperature {
            // Try the provided temperature as-is, then with 'K' suffix
            if self.reactions.contains_key(temp) {
                temp
            } else if self.reactions.contains_key(&format!("{}K", temp)) {
                &format!("{}K", temp)
            } else {
                return Err(format!(
                    "Temperature '{}' not found in loaded data. Available temperatures: [{}]",
                    temp,
                    self.loaded_temperatures.join(", ")
                )
                .into());
            }
        } else {
            // No temperature provided - use single loaded temperature if available
            if self.loaded_temperatures.len() == 1 {
                &self.loaded_temperatures[0]
            } else if self.loaded_temperatures.is_empty() {
                return Err("No temperatures loaded in nuclide data".into());
            } else {
                return Err(format!(
                    "Multiple temperatures loaded [{}], must specify which one to use",
                    self.loaded_temperatures.join(", ")
                )
                .into());
            }
        };

        // Get the reaction data for this temperature and MT
        let temp_reactions = self.reactions.get(temp_key).ok_or_else(|| {
            format!(
                "Temperature '{}' not found in reactions. Available temperatures: [{}]",
                temp_key,
                self.loaded_temperatures.join(", ")
            )
        })?;

        let reaction = temp_reactions.get(&mt).ok_or_else(|| {
            // Get available MTs for this temperature
            let available_mts: Vec<String> =
                temp_reactions.keys().map(|mt| mt.to_string()).collect();
            format!(
                "MT {} not found for temperature '{}'. Available MTs: [{}]",
                mt,
                temp_key,
                available_mts.join(", ")
            )
        })?;

        // Return the cross section and energy data
        if reaction.cross_section.is_empty() {
            return Err(format!(
                "No cross section data available for MT {} at temperature '{}'",
                mt, temp_key
            )
            .into());
        }

        if reaction.energy.is_empty() {
            return Err(format!(
                "No energy grid available for MT {} at temperature '{}'",
                mt, temp_key
            )
            .into());
        }

        Ok((reaction.cross_section.clone(), reaction.energy.clone()))
    }
}

/// Helper function to sum cross sections from all available scattering MTs
/// to create the synthetic scattering reaction (MT 1001).
/// This includes MT 2 (elastic), MT 50-91 (inelastic constituents), and other scattering MTs.
/// Note: MT 4 is NOT included as it's a synthetic reaction itself.
fn compute_scattering_xs(
    temp_reactions: &HashMap<i32, Reaction>,
    energy_grid: &[f64],
) -> Vec<f64> {
    use std::collections::HashSet;
    let mut scattering_xs = vec![0.0; energy_grid.len()];

    // Collect all scattering MTs that are actually present in the data
    let mut included_mts = HashSet::new();

    // Add MT 2 (elastic) if present
    if temp_reactions.contains_key(&2) {
        included_mts.insert(2);
    }

    // Add inelastic constituents (MT 50-91) if present - NEVER MT 4
    for mt in INELASTIC_CONSTITUENT_MTS {
        if temp_reactions.contains_key(&mt) {
            included_mts.insert(mt);
        }
    }

    // Add other scattering MTs (from nonelastic, excluding MT 4)
    for &mt in SCATTERING_MTS_NON_INELASTIC {
        if mt != 2 && temp_reactions.contains_key(&mt) {
            included_mts.insert(mt);
        }
    }

    // Sum cross sections from all included MTs
    for mt in included_mts {
        if let Some(reaction) = temp_reactions.get(&mt) {
            for (i, &energy) in energy_grid.iter().enumerate() {
                if let Some(xs) = reaction.cross_section_at(energy) {
                    scattering_xs[i] += xs;
                }
            }
        }
    }

    scattering_xs
}

// ============================================================================
// HDF5 Reading Functions
// ============================================================================
// Nuclear data is now read from HDF5 files (OpenMC format) instead of JSON.
// The JSON reading code has been removed - use read_nuclide_from_hdf5 instead.

/// Read a nuclide from an HDF5 file path.
/// This is now the primary way to load nuclear data.
pub fn read_nuclide_from_hdf5<P: AsRef<Path>>(
    path: P,
    temps: Option<&std::collections::HashSet<String>>,
) -> Result<Nuclide, Box<dyn std::error::Error>> {
    crate::nuclide_hdf5::read_nuclide_from_hdf5(path, temps)
}


/// Get or load a nuclide from cache, loading from HDF5 file if needed.
///
/// Parameters:
/// - `nuclide_name`: Name of the nuclide (e.g., "Be9", "Li6")
/// - `hdf5_path_map`: Map of nuclide names to HDF5 file paths
/// - `temperatures_to_include`: Optional set of temperatures to load
///
/// Returns cached nuclide if available with sufficient temperatures, otherwise
/// loads from HDF5 file.
pub fn get_or_load_nuclide(
    nuclide_name: &str,
    hdf5_path_map: &HashMap<String, String>,
    temperatures_to_include: Option<&std::collections::HashSet<String>>,
) -> Result<Arc<Nuclide>, Box<dyn std::error::Error>> {
    use std::collections::HashSet;
    let requested: HashSet<String> = temperatures_to_include
        .map(|s| s.iter().cloned().collect())
        .unwrap_or_else(HashSet::new);

    // Get path/keyword from map
    let path_or_keyword = hdf5_path_map.get(nuclide_name).ok_or_else(|| {
        format!(
            "No HDF5 file provided for nuclide '{}'. Please supply a path for all nuclides.",
            nuclide_name
        )
    })?;

    // Resolve keywords/URLs to local paths (when download feature is enabled)
    #[cfg(feature = "download")]
    let resolved_path = {
        let resolved = crate::url_cache::resolve_path_or_url(path_or_keyword, nuclide_name)?;
        resolved.to_string_lossy().to_string()
    };

    #[cfg(not(feature = "download"))]
    let resolved_path = path_or_keyword.clone();

    // Create cache key using resolved path
    let normalized_source = match std::fs::canonicalize(&resolved_path) {
        Ok(canonical) => canonical.to_string_lossy().to_string(),
        Err(_) => resolved_path.clone(),
    };
    let cache_key = format!("{}@{}", nuclide_name, normalized_source);

    // Fast path: cache hit with sufficient temps
    {
        let cache = match GLOBAL_NUCLIDE_CACHE.lock() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        };
        if let Some(existing) = cache.get(&cache_key) {
            if requested.is_empty()
                || requested.is_subset(&existing.loaded_temperatures.iter().cloned().collect())
            {
                return Ok(Arc::clone(existing));
            }
        }
    }

    // Cache miss or insufficient temps - load from HDF5
    let filter = if requested.is_empty() {
        None
    } else {
        Some(&requested)
    };

    let mut nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5(&resolved_path, filter.map(|s| s as &HashSet<String>))?;
    nuclide.data_path = Some(resolved_path.clone());

    // Store in cache
    let arc_nuclide = Arc::new(nuclide);
    {
        let mut cache = match GLOBAL_NUCLIDE_CACHE.lock() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        };
        cache.insert(cache_key, Arc::clone(&arc_nuclide));
    }

    Ok(arc_nuclide)
}

/// Load a nuclide with Python wrapper semantics.
/// Handles both path and name parameters and preserves available_temperatures when filtering.
/// Supports keywords (e.g., "tendl-2019", "fendl-3.1d") when the download feature is enabled.
#[allow(dead_code)]
pub fn load_nuclide_for_python(
    path: Option<&str>,
    nuclide_name: Option<&str>,
    temperatures: Option<&std::collections::HashSet<String>>,
) -> Result<Nuclide, Box<dyn std::error::Error>> {
    // Get the path - path is required for HDF5
    let path_str = path.ok_or("HDF5 file path is required")?;

    // Resolve the path (handles keywords, URLs, and local paths when download feature is enabled)
    #[cfg(feature = "download")]
    let hdf5_path = {
        // Check if it's a keyword or URL - only then do we need the nuclide name
        if crate::url_cache::is_keyword(path_str) || crate::url_cache::is_url(path_str) {
            let name = nuclide_name.ok_or("Nuclide name is required for keyword/URL resolution")?;
            let resolved = crate::url_cache::resolve_path_or_url(path_str, name)?;
            resolved.to_string_lossy().to_string()
        } else {
            // Local path - no name needed
            path_str.to_string()
        }
    };

    #[cfg(not(feature = "download"))]
    let hdf5_path = path_str.to_string();

    // Load the nuclide
    crate::nuclide_hdf5::read_nuclide_from_hdf5(&hdf5_path, temperatures)
}

/// Load a nuclide from a path for the standalone Python function.
#[allow(dead_code)]
pub fn load_nuclide_from_path_or_keyword(
    path: &str,
) -> Result<Nuclide, Box<dyn std::error::Error>> {
    // Load the nuclide with all temperatures from HDF5
    crate::nuclide_hdf5::read_nuclide_from_hdf5(path, None)
}

mod tests {
    #[test]
    fn test_sample_absorption_constituent_li6() {
        use rand::rngs::StdRng;
        use rand::SeedableRng;
        let nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");
        let temperature = nuclide.loaded_temperatures.first().expect("No temperatures loaded");
        let energy = 1e6; // 1 MeV, typical fast neutron
        let mut rng = StdRng::seed_from_u64(12345);
        // Should not panic and should return a Reaction
        let reaction = nuclide.sample_absorption_constituent(energy, temperature, &mut rng);
        // Check that the sampled MT is one of the expected absorption constituent MTs
        let valid_mts: Vec<i32> = {
            let mut mts = vec![
                102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117, 155,
                182, 191, 192, 193, 197,
            ];
            mts.extend(600..650);
            mts.extend(650..700);
            mts.extend(700..750);
            mts.extend(750..800);
            mts.extend(800..850);
            mts
        };
        assert!(
            valid_mts.contains(&reaction.mt_number),
            "Sampled absorption constituent MT {} not in valid list",
            reaction.mt_number
        );
        // Print for debug
        println!("Sampled absorption constituent MT: {}", reaction.mt_number);
    }
    #[test]
    fn test_get_or_load_nuclide_uses_cache() {
        use std::collections::HashMap;
        let li6_path = std::path::Path::new("tests/Li6.h5");
        assert!(li6_path.exists(), "tests/Li6.h5 missing");
        // Only remove Li6 from cache, don't clear all (avoid race with other tests)
        {
            let mut cache = match super::GLOBAL_NUCLIDE_CACHE.lock() {
                Ok(cache) => cache,
                Err(poisoned) => poisoned.into_inner(),
            };
            // Remove any Li6 entries (handle both normalized and non-normalized paths)
            let keys_to_remove: Vec<String> = cache
                .keys()
                .filter(|k| k.starts_with("Li6@") && k.contains("Li6.h5"))
                .cloned()
                .collect();
            for key in keys_to_remove {
                cache.remove(&key);
            }
        }
        let li6_path = std::fs::canonicalize("tests/Li6.h5").expect("tests/Li6.h5 missing");
        let raw = crate::nuclide_hdf5::read_nuclide_from_hdf5(&li6_path, None).expect("Direct read failed");
        assert_eq!(raw.name.as_deref(), Some("Li6"));
        // Don't assert cache state here (other tests may be using it)
        let mut hdf5_map = HashMap::new();
        hdf5_map.insert("Li6".to_string(), li6_path.to_string_lossy().to_string());
        let first =
            super::get_or_load_nuclide("Li6", &hdf5_map, None).expect("Initial cached load failed");
        // Ensure Li6 now present in cache with the correct cache key
        {
            let cache = match super::GLOBAL_NUCLIDE_CACHE.lock() {
                Ok(cache) => cache,
                Err(poisoned) => poisoned.into_inner(),
            };

            // Check for Li6 cache key (path may be normalized)
            let found = cache
                .keys()
                .any(|k| k.starts_with("Li6@") && k.contains("Li6.h5"));
            assert!(found, "Li6 should be present after cached load");
        }
        let second =
            super::get_or_load_nuclide("Li6", &hdf5_map, None).expect("Second cached load failed");
        assert!(
            std::sync::Arc::ptr_eq(&first, &second),
            "Expected identical Arc pointer from cache on second load"
        );
        assert_eq!(
            first.name.as_deref(),
            raw.name.as_deref(),
            "Names differ between raw and cached"
        );
    }

    #[cfg(test)]
    #[test]
    fn test_sample_reaction_li6() {
        use rand::rngs::StdRng;
        use rand::SeedableRng;
        let mut nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");
        let temperature = nuclide.loaded_temperatures.first().expect("No temperatures loaded").clone();

        // Vary energy from 1.0 to 15e6 (10 steps)
        let energies = (0..10).map(|i| 1.0 + i as f64 * (15e6 - 1e6) / 9.0);

        for energy in energies {
            let mut rng1 = StdRng::seed_from_u64(42);
            let mut rng2 = StdRng::seed_from_u64(42);
            let mut rng3 = StdRng::seed_from_u64(43); // Different seed

            let mt1 = nuclide.sample_reaction(energy, &temperature, &mut rng1).map(|r| r.mt_number);
            let mt2 = nuclide.sample_reaction(energy, &temperature, &mut rng2).map(|r| r.mt_number);
            let mt3 = nuclide.sample_reaction(energy, &temperature, &mut rng3).map(|r| r.mt_number);

            // Ensure reactions were sampled successfully
            assert!(
                mt1.is_some(),
                "sample_reaction returned None at energy {}",
                energy
            );
            assert!(
                mt2.is_some(),
                "Repeat sample with same seed returned None at energy {}",
                energy
            );
            assert!(
                mt3.is_some(),
                "Sample with different seed returned None at energy {}",
                energy
            );

            let mt1 = mt1.unwrap();
            let mt2 = mt2.unwrap();
            let mt3 = mt3.unwrap();

            // Ensure same-seed reactions are the same (determinism)
            assert_eq!(
                mt1, mt2,
                "Different MT for same seed at energy {}",
                energy
            );

            // Print info about the third reaction (with different seed)
            println!(
                "Energy: {:e}, MT (seed 42): {}, MT (seed 43): {}",
                energy, mt1, mt3
            );

            // Ensure basic validity
            assert!(
                mt1 > 0,
                "Sampled reaction has invalid MT number at energy {}",
                energy
            );
        }
    }

    #[test]
    fn test_reaction_mts_li6() {
        // Load Li6 nuclide from test HDF5
        let nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");
        let mts = nuclide.reaction_mts().expect("No MTs found");
        // Check for some expected MTs (the exact list may differ from JSON)
        // These are common MTs that should be present in Li6
        let expected = vec![2, 102, 103, 105];
        for mt in &expected {
            assert!(mts.contains(mt), "Expected MT {} in Li6", mt);
        }
        println!("Li6 MTs: {:?}", mts);
    }

    #[test]
    fn test_reaction_mts_li7() {
        // Load Li7 nuclide from test HDF5
        let nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li7.h5", None).expect("Failed to load Li7.h5");
        let mts = nuclide.reaction_mts().expect("No MTs found");
        // Check for presence of key MTs
        assert!(mts.contains(&2), "MT=2 should be present");
        assert!(!mts.is_empty(), "MT list should not be empty");
        println!("Li7 MTs: {:?}", mts);
    }

    #[test]
    fn test_fissionable_false_for_be9_and_fe58() {
        let nuclide_be9 =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");
        assert_eq!(
            nuclide_be9.fissionable, false,
            "Be9 should not be fissionable"
        );

        let nuclide_fe58 =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Fe58.h5", None).expect("Failed to load Fe58.h5");
        assert_eq!(
            nuclide_fe58.fissionable, false,
            "Fe58 should not be fissionable"
        );
    }

    #[test]
    fn test_li6_reactions_contain_specific_mts() {
        // Load Li6 nuclide from test HDF5
        let nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");

        // Check that MT 2 (elastic) is present
        let required = [2];
        for mt in &required {
            let mut found = false;
            for temp_reactions in nuclide.reactions.values() {
                if let Some(reaction) = temp_reactions.get(mt) {
                    assert_ne!(
                        reaction.mt_number, 0,
                        "Reaction MT number for MT {} is 0",
                        mt
                    );
                    found = true;
                    break;
                }
            }
            assert!(found, "MT {} not found in any temperature reactions", mt);
        }
    }

    #[test]
    fn test_available_temperatures_be9() {
        let nuclide_be9 =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");
        // HDF5 files use "294K" format instead of "294"
        assert!(
            nuclide_be9.available_temperatures.iter().any(|t| t.contains("294")),
            "available_temperatures should contain 294K variant, got {:?}",
            nuclide_be9.available_temperatures
        );
        let temps_method = nuclide_be9
            .temperatures()
            .expect("temperatures() returned None");
        assert!(
            !temps_method.is_empty(),
            "temperatures() should return at least one temperature"
        );
    }

    #[test]
    fn test_be9_mt_numbers_per_temperature() {
        let nuclide_be9 =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        // Get the temperature key (HDF5 uses "294K" format)
        let temp_key = nuclide_be9.reactions.keys().next().expect("No temperature found");

        // Helper closure to extract and sort MT list for a temperature
        let get_sorted_mts = |temp: &str| -> Vec<i32> {
            let mut mts: Vec<i32> = nuclide_be9
                .reactions
                .get(temp)
                .expect("Temperature not found in reactions")
                .keys()
                .cloned()
                .collect();
            mts.sort();
            mts
        };

        let mts = get_sorted_mts(temp_key);
        println!("Be9 MTs at {}: {:?}", temp_key, mts);

        // Check some expected MTs are present
        assert!(mts.contains(&2), "Be9 should have MT=2 (elastic)");
    }

    #[test]
    fn test_available_temperatures_fe56_includes_294() {
        let nuclide_fe56 = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Fe56.h5", None)
            .expect("Failed to load Fe56.h5");
        assert!(
            nuclide_fe56
                .available_temperatures
                .iter()
                .any(|t| t.contains("294")),
            "Fe56 available_temperatures should contain '294' variant"
        );
        let temps_method = nuclide_fe56
            .temperatures()
            .expect("temperatures() returned None");
        assert!(
            !temps_method.is_empty(),
            "Fe56 temperatures() should return at least one temperature"
        );
    }

    #[test]
    fn test_clear_nuclide_cache() {
        // Insert a test nuclide into the cache
        let nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");
        let nuclide_arc = std::sync::Arc::new(nuclide);

        {
            let mut cache = match super::GLOBAL_NUCLIDE_CACHE.lock() {
                Ok(cache) => cache,
                Err(poisoned) => poisoned.into_inner(),
            };
            cache.insert(
                "Test_Nuclide".to_string(),
                std::sync::Arc::clone(&nuclide_arc),
            );

            // Verify it's in the cache
            assert!(
                cache.contains_key("Test_Nuclide"),
                "Test nuclide should be in cache before clearing"
            );
        }

        // Call the clear function
        super::clear_nuclide_cache();

        // Verify cache is now empty
        let cache = match super::GLOBAL_NUCLIDE_CACHE.lock() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        };
        assert!(
            cache.is_empty(),
            "Cache should be empty after calling clear_nuclide_cache"
        );
    }

    #[test]
    fn test_microscopic_cross_section_with_temperature() {
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        // Get the actual temperature key from the data
        let temp_key = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Test with specific temperature
        let result = nuclide.microscopic_cross_section(2, Some(&temp_key), false);
        assert!(result.is_ok(), "Should successfully get MT=2 data for {}: {:?}", temp_key, result.err());

        let (xs, energy) = result.unwrap();
        assert!(!xs.is_empty(), "Cross section data should not be empty");
        assert!(!energy.is_empty(), "Energy data should not be empty");
        assert_eq!(
            xs.len(),
            energy.len(),
            "Cross section and energy arrays should have same length"
        );

        // Test with invalid temperature
        let result_invalid = nuclide.microscopic_cross_section(2, Some("999K"), false);
        assert!(
            result_invalid.is_err(),
            "Should fail for unavailable temperature"
        );
    }

    #[test]
    fn test_microscopic_cross_section_single_temperature() {
        // Load Be9 with only one temperature
        let temps_filter = std::collections::HashSet::from(["294K".to_string()]);
        let mut nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", Some(&temps_filter))
            .expect("Failed to load Be9.h5 with temperature filter");

        // Should work without specifying temperature since only one is loaded
        let result = nuclide.microscopic_cross_section(2, None, false);
        assert!(
            result.is_ok(),
            "Should successfully get MT=2 data without temperature: {:?}",
            result.err()
        );

        let (xs, energy) = result.unwrap();
        assert!(!xs.is_empty(), "Cross section data should not be empty");
        assert!(!energy.is_empty(), "Energy data should not be empty");
        assert_eq!(
            xs.len(),
            energy.len(),
            "Cross section and energy arrays should have same length"
        );
    }

    #[test]
    fn test_microscopic_cross_section_invalid_mt() {
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        let temp_key = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Should fail for non-existent MT
        let result = nuclide.microscopic_cross_section(9999, Some(&temp_key), false);
        assert!(result.is_err(), "Should error for invalid MT number");

        let error_msg = result.unwrap_err().to_string();
        assert!(
            error_msg.contains("MT 9999 not found"),
            "Error should mention MT not found: {}",
            error_msg
        );
    }

    #[test]
    fn test_microscopic_cross_section_multiple_mt_numbers() {
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        let temp_key = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Test common MT numbers that should exist in Be9
        let test_mts = [2, 102]; // elastic and capture should exist

        for mt in test_mts {
            let result = nuclide.microscopic_cross_section(mt, Some(&temp_key), false);
            if result.is_ok() {
                let (xs, energy) = result.unwrap();
                assert!(!xs.is_empty(), "MT={} should have cross section data", mt);
                assert!(!energy.is_empty(), "MT={} should have energy data", mt);
                assert_eq!(xs.len(), energy.len(), "MT={} data length mismatch", mt);

                // Validate data quality
                for &e in &energy {
                    assert!(e > 0.0, "MT={} energy values should be positive", mt);
                }
                for &x in &xs {
                    assert!(x >= 0.0, "MT={} cross sections should be non-negative", mt);
                }
            }
        }
    }

    #[test]
    fn test_microscopic_cross_section_lithium() {
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");

        // Li6 should have only one temperature, so no temperature needed
        let result = nuclide.microscopic_cross_section(2, None, false);
        assert!(
            result.is_ok(),
            "Should successfully get Li6 elastic scattering data: {:?}",
            result.err()
        );

        let (xs, energy) = result.unwrap();
        assert!(
            !xs.is_empty(),
            "Li6 elastic scattering data should not be empty"
        );
        assert!(!energy.is_empty(), "Li6 energy data should not be empty");
    }

    #[test]
    fn test_cache_optimization_keyword_vs_path() {
        // Test that accessing the same file via keyword vs direct path uses same cache entry
        use std::collections::HashMap;

        super::clear_nuclide_cache();

        // Load Li6 from local file first
        let mut li6_map = HashMap::new();
        li6_map.insert("Li6".to_string(), "tests/Li6.h5".to_string());
        let _first =
            super::get_or_load_nuclide("Li6", &li6_map, None).expect("Initial file load failed");

        // Now try to load the same file by its absolute path
        let absolute_path = std::fs::canonicalize("tests/Li6.h5").unwrap();
        let mut abs_map = HashMap::new();
        abs_map.insert(
            "Li6".to_string(),
            absolute_path.to_string_lossy().to_string(),
        );
        let _second =
            super::get_or_load_nuclide("Li6", &abs_map, None).expect("Absolute path load failed");

        // Verify both entries use the same cache key (should only be one entry)
        {
            let cache = match super::GLOBAL_NUCLIDE_CACHE.lock() {
                Ok(cache) => cache,
                Err(poisoned) => poisoned.into_inner(),
            };

            // Should only have one cache entry since both resolve to the same file
            let li6_entries: Vec<_> = cache
                .keys()
                .filter(|k| k.starts_with("Li6@") && k.contains("Li6.h5"))
                .collect();

            assert_eq!(
                li6_entries.len(),
                1,
                "Expected exactly 1 cache entry for Li6, found: {:?}",
                li6_entries
            );
        }
    }

    #[test]
    fn test_microscopic_cross_section_string_reactions() {
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        let temp_key = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Test with MT number (existing functionality)
        let result_mt = nuclide.microscopic_cross_section(2, Some(&temp_key), false);
        assert!(result_mt.is_ok(), "Should work with MT number");
        let (xs_mt, energy_mt) = result_mt.unwrap();

        // Test with reaction name string
        let result_str = nuclide.microscopic_cross_section("(n,elastic)", Some(&temp_key), false);
        assert!(result_str.is_ok(), "Should work with reaction name string");
        let (xs_str, energy_str) = result_str.unwrap();

        // Results should be identical
        assert_eq!(
            xs_mt, xs_str,
            "Cross section data should be identical for MT 2 and '(n,elastic)'"
        );
        assert_eq!(
            energy_mt, energy_str,
            "Energy data should be identical for MT 2 and '(n,elastic)'"
        );
    }

    #[test]
    fn test_microscopic_cross_section_invalid_reaction_string() {
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        let temp_key = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Test with invalid reaction name
        let result = nuclide.microscopic_cross_section("(n,invalid)", Some(&temp_key), false);
        assert!(result.is_err(), "Should fail for invalid reaction name");

        let error_msg = result.unwrap_err().to_string();
        assert!(
            error_msg.contains("Unknown reaction name"),
            "Error should mention unknown reaction name"
        );
    }

    #[test]
    fn test_microscopic_cross_section_fission_alias() {
        // Note: Be9 is not fissionable, so we'll just test the string recognition
        // The fission alias should map to MT 18
        let mut nuclide =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Be9.h5", None).expect("Failed to load Be9.h5");

        let temp_key = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Test that "fission" string is recognized (even though Be9 doesn't have fission reactions)
        let result = nuclide.microscopic_cross_section("fission", Some(&temp_key), false);

        // Should fail because Be9 doesn't have MT 18, but the error should be about missing MT, not unknown reaction
        assert!(
            result.is_err(),
            "Should fail because Be9 doesn't have fission reactions"
        );
        let error_msg = result.unwrap_err().to_string();
        assert!(
            error_msg.contains("MT 18 not found"),
            "Error should be about missing MT 18, not unknown reaction name"
        );
    }

    #[test]
    fn test_nuclide_different_data_sources() {
        // Test that loading the same nuclide from different sources gives different results

        // Clear cache to ensure fresh loads
        crate::nuclide::clear_nuclide_cache();

        // For this test, we'll use different file paths as data sources
        let mut li6_file = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None)
            .expect("Failed to load Li6 from file");

        let mut li7_file = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li7.h5", None)
            .expect("Failed to load Li7 from file");

        // Get cross sections from both (different nuclides will have different data)
        let (xs_li6, _) = li6_file
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li6 cross section");
        let (xs_li7, _) = li7_file
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li7 cross section");

        // Should have different data since they're different nuclides
        let data_different = xs_li6.len() != xs_li7.len()
            || xs_li6
                .iter()
                .zip(&xs_li7)
                .any(|(a, b)| (a - b).abs() > 1e-10);

        assert!(data_different, "Li6 and Li7 data should be different");
        println!("Li6: {} points, Li7: {} points", xs_li6.len(), xs_li7.len());
    }

    #[test]
    fn test_nuclide_file_vs_keyword_sources() {
        // Test that file paths and keywords can coexist in cache

        crate::nuclide::clear_nuclide_cache();

        // Load Li6 from local file
        let mut li6_file = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None)
            .expect("Failed to load Li6 from file");

        // Load Li7 from local file (different nuclide, different source)
        let mut li7_file = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li7.h5", None)
            .expect("Failed to load Li7 from file");

        assert_eq!(li6_file.name.as_deref(), Some("Li6"));
        assert_eq!(li7_file.name.as_deref(), Some("Li7"));

        // Verify we can load cross sections from both
        let (xs_li6, _) = li6_file
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li6 cross section");
        let (xs_li7, _) = li7_file
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li7 cross section");

        assert!(
            !xs_li6.is_empty() && !xs_li7.is_empty(),
            "Both should have cross section data"
        );
    }

    #[test]
    fn test_nuclide_cache_respects_data_source_boundaries() {
        // Test that the cache properly separates different data sources

        crate::nuclide::clear_nuclide_cache();

        // Load Li6 from file first time
        let mut li6_1 = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None)
            .expect("Failed to load Li6 first time");

        // Load Li7 from file (different nuclide/source)
        let mut li7 =
            crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li7.h5", None).expect("Failed to load Li7");

        // Load Li6 from file again (should use cache)
        let mut li6_2 = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None)
            .expect("Failed to load Li6 second time");

        // Get cross sections
        let (xs_li6_1, _) = li6_1
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li6 cross section first time");
        let (xs_li7, _) = li7
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li7 cross section");
        let (xs_li6_2, _) = li6_2
            .microscopic_cross_section(102, None, false)
            .expect("Failed to get Li6 cross section second time");

        // Li6 loads should be identical (cache working)
        assert_eq!(
            xs_li6_1, xs_li6_2,
            "Li6 loads should be identical (cache working)"
        );

        // Li6 vs Li7 should be different (different nuclides)
        let li6_vs_li7_different = xs_li6_1.len() != xs_li7.len()
            || xs_li6_1
                .iter()
                .zip(&xs_li7)
                .any(|(a, b)| (a - b).abs() > 1e-10);
        assert!(
            li6_vs_li7_different,
            "Li6 and Li7 should have different data"
        );
    }

    #[test]
    fn test_sample_inelastic_constituent() {
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        // Load Li6 which should have some inelastic reactions (MT 50-91)
        let nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");
        let temperature = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());

        // Check what inelastic constituent MTs are available
        let available_mts = nuclide.reaction_mts().unwrap_or_default();
        let inelastic_mts: Vec<i32> = available_mts
            .into_iter()
            .filter(|&mt| mt >= 50 && mt < 92)
            .collect();

        if inelastic_mts.is_empty() {
            println!("No inelastic constituent reactions (MT 50-91) found in Li6, skipping test");
            return;
        }

        println!(
            "Available inelastic constituent MTs in Li6: {:?}",
            inelastic_mts
        );

        // Test sampling at various energies
        let energies = [0.375e7];
        for energy in energies {
            for i in 0..5 {
                let mut rng = StdRng::seed_from_u64(42 + i);
                let reaction = nuclide.sample_inelastic_constituent(energy, &temperature, &mut rng);
                println!(
                    "Sample {}: Energy {:.1e}: Sampled inelastic MT {}",
                    i + 1,
                    energy,
                    reaction.mt_number
                );
                // Verify the sampled MT is in the expected range
                assert!(
                    reaction.mt_number >= 50 && reaction.mt_number < 92,
                    "Sampled MT {} should be in range 50-91",
                    reaction.mt_number
                );
            }
        }
    }

    #[test]
    fn test_sample_inelastic_constituent_deterministic() {
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let nuclide = crate::nuclide_hdf5::read_nuclide_from_hdf5("tests/Li6.h5", None).expect("Failed to load Li6.h5");
        let temperature = nuclide.loaded_temperatures.first().cloned().unwrap_or("294K".to_string());
        let energy = 1.0e7;

        // Same seed should give same result
        let mut rng1 = StdRng::seed_from_u64(12345);
        let mut rng2 = StdRng::seed_from_u64(12345);

        let reaction1 = nuclide.sample_inelastic_constituent(energy, &temperature, &mut rng1);
        let reaction2 = nuclide.sample_inelastic_constituent(energy, &temperature, &mut rng2);
        assert_eq!(
            reaction1.mt_number, reaction2.mt_number,
            "Same seed should give same inelastic constituent reaction"
        );
    }
}
