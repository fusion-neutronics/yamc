// HDF5 reading functions for nuclear data
// Matches OpenMC's HDF5 file format

use hdf5::{File, Group};
use hdf5::types::{VarLenUnicode, VarLenAscii, FixedAscii};
use std::collections::HashMap;
use std::path::Path;


use crate::nuclide::{Nuclide, FissionNuData};
use crate::reaction::Reaction;
use crate::reaction_product::{
    ReactionProduct, AngleEnergyDistribution, AngleDistribution,
    EnergyDistribution, ParticleType, Yield, Tabulated1D,
};
use crate::secondary_kalbach::{KalbachMann, KMTable, Interpolation as KMInterpolation};
use crate::secondary_correlated::{CorrelatedAngleEnergy, CorrTable, Tabular, Interpolation as CorrInterpolation};

/// Read a nuclide from an HDF5 file
pub fn read_nuclide_from_hdf5<P: AsRef<Path>>(
    path: P,
    temps_filter: Option<&std::collections::HashSet<String>>,
) -> Result<Nuclide, Box<dyn std::error::Error>> {
    let file = File::open(path.as_ref())?;

    // Get the nuclide group (first group in the file, e.g., "Be9")
    let nuclide_names: Vec<String> = file.member_names()?;
    if nuclide_names.is_empty() {
        return Err("No nuclide groups found in HDF5 file".into());
    }

    let nuclide_name = &nuclide_names[0];
    let nuclide_group = file.group(nuclide_name)?;

    parse_nuclide_from_hdf5_group(&nuclide_group, nuclide_name, temps_filter)
}

/// Helper to read a string attribute from HDF5
/// Handles various HDF5 string types: variable-length unicode, variable-length ASCII,
/// and fixed-size ASCII strings of various lengths
fn read_string_attr(group: &Group, name: &str) -> Result<String, hdf5::Error> {
    let attr = group.attr(name)?;

    // Try reading as VarLenUnicode first (variable-length string)
    if let Ok(s) = attr.read_scalar::<VarLenUnicode>() {
        return Ok(s.to_string());
    }

    // Try VarLenAscii
    if let Ok(s) = attr.read_scalar::<VarLenAscii>() {
        return Ok(s.to_string());
    }

    // For fixed-size strings, we need to try from LARGEST to SMALLEST
    // because FixedAscii<4> will successfully read a 12-char string but only return 4 chars
    macro_rules! try_fixed_ascii_largest_first {
        ($($n:expr),*) => {
            $(
                if let Ok(s) = attr.read_scalar::<FixedAscii<$n>>() {
                    let trimmed = s.to_string().trim_end_matches('\0').to_string();
                    // Only return if we got a non-empty string
                    if !trimmed.is_empty() {
                        return Ok(trimmed);
                    }
                }
            )*
        }
    }

    // Try from largest to smallest to ensure we get the full string
    try_fixed_ascii_largest_first!(64, 32, 20, 16, 12, 11, 10, 9, 8, 7, 6, 5, 4);

    // If all else fails, return an error
    Err(hdf5::Error::Internal(format!(
        "Could not read string attribute '{}': unknown string type",
        name
    )))
}

/// Helper to read an f64 attribute from HDF5
fn read_f64_attr(group: &Group, name: &str) -> Result<f64, hdf5::Error> {
    let attr = group.attr(name)?;
    attr.read_scalar::<f64>()
}

/// Parse element symbol and name from nuclide name (e.g., "Li6" -> ("Li", "lithium"))
fn parse_element_from_name(name: &str, z: i32) -> (Option<String>, Option<String>) {
    // Extract symbol from name (letters before numbers)
    let symbol: String = name.chars().take_while(|c| c.is_alphabetic()).collect();

    // Element names lookup by atomic number
    const ELEMENTS: &[(i32, &str, &str)] = &[
        (1, "H", "hydrogen"), (2, "He", "helium"), (3, "Li", "lithium"),
        (4, "Be", "beryllium"), (5, "B", "boron"), (6, "C", "carbon"),
        (7, "N", "nitrogen"), (8, "O", "oxygen"), (9, "F", "fluorine"),
        (10, "Ne", "neon"), (11, "Na", "sodium"), (12, "Mg", "magnesium"),
        (13, "Al", "aluminum"), (14, "Si", "silicon"), (15, "P", "phosphorus"),
        (16, "S", "sulfur"), (17, "Cl", "chlorine"), (18, "Ar", "argon"),
        (19, "K", "potassium"), (20, "Ca", "calcium"), (21, "Sc", "scandium"),
        (22, "Ti", "titanium"), (23, "V", "vanadium"), (24, "Cr", "chromium"),
        (25, "Mn", "manganese"), (26, "Fe", "iron"), (27, "Co", "cobalt"),
        (28, "Ni", "nickel"), (29, "Cu", "copper"), (30, "Zn", "zinc"),
        (31, "Ga", "gallium"), (32, "Ge", "germanium"), (33, "As", "arsenic"),
        (34, "Se", "selenium"), (35, "Br", "bromine"), (36, "Kr", "krypton"),
        (37, "Rb", "rubidium"), (38, "Sr", "strontium"), (39, "Y", "yttrium"),
        (40, "Zr", "zirconium"), (41, "Nb", "niobium"), (42, "Mo", "molybdenum"),
        (43, "Tc", "technetium"), (44, "Ru", "ruthenium"), (45, "Rh", "rhodium"),
        (46, "Pd", "palladium"), (47, "Ag", "silver"), (48, "Cd", "cadmium"),
        (49, "In", "indium"), (50, "Sn", "tin"), (51, "Sb", "antimony"),
        (52, "Te", "tellurium"), (53, "I", "iodine"), (54, "Xe", "xenon"),
        (55, "Cs", "cesium"), (56, "Ba", "barium"), (57, "La", "lanthanum"),
        (58, "Ce", "cerium"), (59, "Pr", "praseodymium"), (60, "Nd", "neodymium"),
        (61, "Pm", "promethium"), (62, "Sm", "samarium"), (63, "Eu", "europium"),
        (64, "Gd", "gadolinium"), (65, "Tb", "terbium"), (66, "Dy", "dysprosium"),
        (67, "Ho", "holmium"), (68, "Er", "erbium"), (69, "Tm", "thulium"),
        (70, "Yb", "ytterbium"), (71, "Lu", "lutetium"), (72, "Hf", "hafnium"),
        (73, "Ta", "tantalum"), (74, "W", "tungsten"), (75, "Re", "rhenium"),
        (76, "Os", "osmium"), (77, "Ir", "iridium"), (78, "Pt", "platinum"),
        (79, "Au", "gold"), (80, "Hg", "mercury"), (81, "Tl", "thallium"),
        (82, "Pb", "lead"), (83, "Bi", "bismuth"), (84, "Po", "polonium"),
        (85, "At", "astatine"), (86, "Rn", "radon"), (87, "Fr", "francium"),
        (88, "Ra", "radium"), (89, "Ac", "actinium"), (90, "Th", "thorium"),
        (91, "Pa", "protactinium"), (92, "U", "uranium"), (93, "Np", "neptunium"),
        (94, "Pu", "plutonium"), (95, "Am", "americium"), (96, "Cm", "curium"),
    ];

    // Find element by atomic number
    if let Some((_, _, name)) = ELEMENTS.iter().find(|(num, _, _)| *num == z) {
        (Some(symbol), Some(name.to_string()))
    } else {
        (Some(symbol), None)
    }
}

/// Parse a Nuclide from an HDF5 group
fn parse_nuclide_from_hdf5_group(
    group: &Group,
    name: &str,
    temps_filter: Option<&std::collections::HashSet<String>>,
) -> Result<Nuclide, Box<dyn std::error::Error>> {
    // Read attributes
    let z: i32 = group.attr("Z")?.read_scalar()?;
    let a: i32 = group.attr("A")?.read_scalar()?;
    let awr: f64 = group.attr("atomic_weight_ratio")?.read_scalar()?;
    let _metastable: i32 = group.attr("metastable")?.read_scalar()?;

    // Parse element symbol from name (e.g., "Li6" -> "Li")
    let (element_symbol, element_name) = parse_element_from_name(name, z);

    let mut nuclide = Nuclide {
        name: Some(name.to_string()),
        element: element_name,
        atomic_symbol: element_symbol,
        atomic_number: Some(z as u32),
        neutron_number: Some((a - z) as u32),
        mass_number: Some(a as u32),
        atomic_weight_ratio: Some(awr),
        library: None,
        energy: None,
        reactions: HashMap::new(),
        fissionable: false,
        available_temperatures: Vec::new(),
        loaded_temperatures: Vec::new(),
        data_path: None,
        fission_nu: None,
        scattering_mts: HashMap::new(),
    };

    // Read available temperatures from kTs group
    let kts_group = group.group("kTs")?;
    let kt_names: Vec<String> = kts_group.member_names()?;

    let mut available_temps: Vec<String> = Vec::new();
    for kt_name in &kt_names {
        // kt_name is like "294K"
        available_temps.push(kt_name.clone());
    }
    available_temps.sort();
    nuclide.available_temperatures = available_temps.iter()
        .map(|s| s.trim_end_matches('K').to_string())
        .collect();

    // Determine which temperatures to load
    let temps_to_load: Vec<String> = if let Some(filter) = temps_filter {
        if filter.is_empty() {
            kt_names.clone()
        } else {
            kt_names.iter()
                .filter(|t| {
                    let t_num = t.trim_end_matches('K');
                    filter.contains(t_num) || filter.contains(*t)
                })
                .cloned()
                .collect()
        }
    } else {
        kt_names.clone()
    };

    // Read energy grids
    let energy_group = group.group("energy")?;
    let mut energy_map: HashMap<String, Vec<f64>> = HashMap::new();

    for temp in &temps_to_load {
        if let Ok(dataset) = energy_group.dataset(temp) {
            let energy_data: Vec<f64> = dataset.read_1d()?.to_vec();
            let temp_key = temp.trim_end_matches('K').to_string();
            energy_map.insert(temp_key, energy_data);
        }
    }
    nuclide.energy = Some(energy_map);

    // Read reactions
    let reactions_group = group.group("reactions")?;
    let reaction_names: Vec<String> = reactions_group.member_names()?;

    for temp in &temps_to_load {
        let temp_key = temp.trim_end_matches('K').to_string();
        let mut temp_reactions: HashMap<i32, Reaction> = HashMap::new();

        for rx_name in &reaction_names {
            if !rx_name.starts_with("reaction_") {
                continue;
            }

            let rx_group = reactions_group.group(rx_name)?;

            // Read reaction attributes
            let mt: i32 = rx_group.attr("mt")?.read_scalar()?;
            let q_value: f64 = rx_group.attr("Q_value")?.read_scalar()?;

            // Read center_of_mass flag (whether scattering is in CM frame)
            let scatter_in_cm: bool = rx_group.attr("center_of_mass")
                .ok()
                .and_then(|attr| attr.read_scalar::<i32>().ok())
                .map(|val| val == 1)
                .unwrap_or(false);

            // Check for fission
            if is_fission(mt) {
                nuclide.fissionable = true;
            }

            // Read cross section for this temperature
            if let Ok(temp_group) = rx_group.group(temp) {
                let xs_dataset = temp_group.dataset("xs")?;
                let xs_data: Vec<f64> = xs_dataset.read_1d()?.to_vec();

                // Read threshold_idx attribute
                let threshold_idx: usize = xs_dataset.attr("threshold_idx")
                    .map(|a| a.read_scalar::<i32>().unwrap_or(0) as usize)
                    .unwrap_or(0);

                // Calculate energy grid for this reaction
                let energy_grid = if let Some(ref e_map) = nuclide.energy {
                    if let Some(full_grid) = e_map.get(&temp_key) {
                        if threshold_idx < full_grid.len() {
                            full_grid[threshold_idx..].to_vec()
                        } else {
                            Vec::new()
                        }
                    } else {
                        Vec::new()
                    }
                } else {
                    Vec::new()
                };

                // Read products
                let products = read_products(&rx_group)?;

                let reaction = Reaction {
                    cross_section: xs_data,
                    threshold_idx,
                    interpolation: vec![2], // lin-lin
                    energy: energy_grid,
                    mt_number: mt,
                    q_value,
                    products,
                    scatter_in_cm,
                };

                temp_reactions.insert(mt, reaction);
            }
        }

        if !temp_reactions.is_empty() {
            nuclide.reactions.insert(temp_key.clone(), temp_reactions);
        }
    }

    // Read total_nu (nu-bar) data if present (for fissionable nuclides)
    if let Ok(total_nu_group) = group.group("total_nu") {
        if let Ok(yield_ds) = total_nu_group.dataset("yield") {
            // yield is shape (2, N) where row 0 is energy, row 1 is nu values
            let yield_shape = yield_ds.shape();
            if yield_shape.len() == 2 && yield_shape[0] == 2 {
                let n = yield_shape[1];
                let yield_data: Vec<f64> = yield_ds.read_raw::<f64>().unwrap_or_default();
                if yield_data.len() == 2 * n {
                    let energy: Vec<f64> = yield_data[0..n].to_vec();
                    let nu: Vec<f64> = yield_data[n..2*n].to_vec();
                    if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
                        eprintln!("  Loaded total_nu: {} points, nu range [{:.2}, {:.2}]",
                            n, nu.first().unwrap_or(&0.0), nu.last().unwrap_or(&0.0));
                    }
                    nuclide.fission_nu = Some(FissionNuData { energy, nu });
                }
            }
        }
    }

    // Synthesize hierarchical MTs for each temperature
    synthesize_hierarchical_mts(&mut nuclide);

    // Cache scattering MTs for each temperature (optimization for sample_scattering_constituent)
    populate_scattering_mts_cache(&mut nuclide);

    // Set loaded temperatures
    nuclide.loaded_temperatures = nuclide.reactions.keys().cloned().collect();
    nuclide.loaded_temperatures.sort();

    Ok(nuclide)
}

/// Synthesize hierarchical MTs (1, 3, 4, 27, 101, 1001) from constituent reactions
/// This is an optimized version that pre-interpolates reactions onto the full energy grid
fn synthesize_hierarchical_mts(nuclide: &mut Nuclide) {
    // Get the full energy grid for each temperature
    let Some(ref energy_map) = nuclide.energy else { return };

    // Collect new reactions to add (to avoid borrow issues)
    let mut new_reactions: Vec<(String, i32, Reaction)> = Vec::new();

    for (temp_key, temp_reactions) in nuclide.reactions.iter() {
        let Some(full_energy) = energy_map.get(temp_key) else { continue };
        let n = full_energy.len();
        if n == 0 { continue }

        // Pre-interpolate all reactions onto the full energy grid (O(n * num_reactions))
        // This avoids repeated binary searches later
        let mut xs_cache: HashMap<i32, Vec<f64>> = HashMap::new();
        for (&mt, rx) in temp_reactions.iter() {
            let xs_vec: Vec<f64> = full_energy.iter()
                .map(|&e| rx.cross_section_at(e).unwrap_or(0.0))
                .collect();
            xs_cache.insert(mt, xs_vec);
        }

        // Helper to get cached cross section at index
        let get_xs = |mt: i32, i: usize| -> f64 {
            xs_cache.get(&mt).and_then(|v| v.get(i).copied()).unwrap_or(0.0)
        };

        // MT 4: Total inelastic (sum of MT 50-91)
        let mut mt4_xs = vec![0.0; n];
        for i in 0..n {
            for mt in 50..92 {
                mt4_xs[i] += get_xs(mt, i);
            }
        }
        if mt4_xs.iter().any(|&x| x > 0.0) {
            new_reactions.push((temp_key.clone(), 4, Reaction {
                cross_section: mt4_xs.clone(),
                threshold_idx: 0,
                interpolation: vec![2],
                energy: full_energy.clone(),
                mt_number: 4,
                q_value: 0.0,
                products: Vec::new(),
                scatter_in_cm: false,
            }));
        }

        // MT 101: Disappearance (sum of capture/charged particle reactions)
        let disappearance_mts: Vec<i32> = (102..118).chain(155..156).chain(182..200).chain(600..850).collect();
        let mut mt101_xs = vec![0.0; n];
        for i in 0..n {
            for &mt in &disappearance_mts {
                mt101_xs[i] += get_xs(mt, i);
            }
        }
        if mt101_xs.iter().any(|&x| x > 0.0) {
            new_reactions.push((temp_key.clone(), 101, Reaction {
                cross_section: mt101_xs.clone(),
                threshold_idx: 0,
                interpolation: vec![2],
                energy: full_energy.clone(),
                mt_number: 101,
                q_value: 0.0,
                products: Vec::new(),
                scatter_in_cm: false,
            }));
        }

        // MT 27: Absorption (fission + disappearance)
        let mut mt27_xs = vec![0.0; n];
        for i in 0..n {
            mt27_xs[i] = get_xs(18, i) + mt101_xs[i];
        }
        if mt27_xs.iter().any(|&x| x > 0.0) {
            new_reactions.push((temp_key.clone(), 27, Reaction {
                cross_section: mt27_xs,
                threshold_idx: 0,
                interpolation: vec![2],
                energy: full_energy.clone(),
                mt_number: 27,
                q_value: 0.0,
                products: Vec::new(),
                scatter_in_cm: false,
            }));
        }

        // MT 3: Non-elastic (sum of all non-elastic reactions, exclude MT 2)
        let mut non_elastic_xs = vec![0.0; n];
        for i in 0..n {
            for (&mt, _) in temp_reactions.iter() {
                if mt != 2 && mt < 100 && mt != 3 && mt != 1 && mt != 4 && mt != 27 {
                    non_elastic_xs[i] += get_xs(mt, i);
                }
            }
            // Add disappearance
            non_elastic_xs[i] += mt101_xs[i];
        }
        if non_elastic_xs.iter().any(|&x| x > 0.0) {
            new_reactions.push((temp_key.clone(), 3, Reaction {
                cross_section: non_elastic_xs.clone(),
                threshold_idx: 0,
                interpolation: vec![2],
                energy: full_energy.clone(),
                mt_number: 3,
                q_value: 0.0,
                products: Vec::new(),
                scatter_in_cm: false,
            }));
        }

        // MT 1: Total (elastic + non-elastic)
        let mut mt1_xs = vec![0.0; n];
        for i in 0..n {
            mt1_xs[i] = get_xs(2, i) + non_elastic_xs[i];
        }
        if mt1_xs.iter().any(|&x| x > 0.0) {
            new_reactions.push((temp_key.clone(), 1, Reaction {
                cross_section: mt1_xs,
                threshold_idx: 0,
                interpolation: vec![2],
                energy: full_energy.clone(),
                mt_number: 1,
                q_value: 0.0,
                products: Vec::new(),
                scatter_in_cm: false,
            }));
        }

        // MT 1001: Synthetic scattering (elastic + inelastic)
        let mut mt1001_xs = vec![0.0; n];
        for i in 0..n {
            mt1001_xs[i] = get_xs(2, i) + mt4_xs[i];
        }
        if mt1001_xs.iter().any(|&x| x > 0.0) {
            new_reactions.push((temp_key.clone(), 1001, Reaction {
                cross_section: mt1001_xs,
                threshold_idx: 0,
                interpolation: vec![2],
                energy: full_energy.clone(),
                mt_number: 1001,
                q_value: 0.0,
                products: Vec::new(),
                scatter_in_cm: false,
            }));
        }
    }

    // Insert all new reactions
    for (temp_key, mt, reaction) in new_reactions {
        if let Some(temp_reactions) = nuclide.reactions.get_mut(&temp_key) {
            temp_reactions.insert(mt, reaction);
        }
    }
}

/// Pre-populate the scattering_mts cache for each temperature.
/// This caches which MTs are scattering reactions to avoid iterating all MTs
/// during sample_scattering_constituent calls.
fn populate_scattering_mts_cache(nuclide: &mut Nuclide) {
    for (temp_key, temp_reactions) in nuclide.reactions.iter() {
        let mut scattering_mts: Vec<i32> = temp_reactions
            .keys()
            .filter(|&&mt| {
                // Skip synthetic MTs (1, 4, 18, 101, 1001) - same logic as sample_scattering_constituent
                if mt == 4 || mt == 1 || mt == 18 || mt == 101 || mt == 1001 {
                    return false;
                }
                crate::nuclide::is_scattering_mt(mt)
            })
            .copied()
            .collect();
        // Sort for deterministic iteration order
        scattering_mts.sort();
        nuclide.scattering_mts.insert(temp_key.clone(), scattering_mts);
    }
}

/// Read products from a reaction group
fn read_products(rx_group: &Group) -> Result<Vec<ReactionProduct>, Box<dyn std::error::Error>> {
    let mut products = Vec::new();

    let member_names: Vec<String> = rx_group.member_names()?;
    for name in member_names {
        if !name.starts_with("product_") {
            continue;
        }

        let prod_group = rx_group.group(&name)
            .map_err(|e| format!("Failed to open product group {}: {}", name, e))?;

        // Read particle type
        let particle_str = read_string_attr(&prod_group, "particle")?;
        let particle = match particle_str.as_str() {
            "neutron" => ParticleType::Neutron,
            "photon" => ParticleType::Photon,
            _ => ParticleType::Neutron, // Default to neutron
        };

        // Read emission mode
        let emission_mode = read_string_attr(&prod_group, "emission_mode").unwrap_or_else(|_| "prompt".to_string());

        // Read yield
        // Yield can be:
        // - shape (1,) for constant yield
        // - shape (2, N) for tabulated yield where row 0 is energy, row 1 is yield values
        let yield_ds = prod_group.dataset("yield")?;
        let yield_shape = yield_ds.shape();
        let product_yield = if yield_shape.len() == 1 {
            // Constant or polynomial yield
            let yield_data: Vec<f64> = yield_ds.read_1d()?.to_vec();
            if !yield_data.is_empty() {
                Some(Yield::Polynomial { coefficients: yield_data })
            } else {
                None
            }
        } else if yield_shape.len() == 2 {
            // Tabulated yield: shape (2, N) where row 0 is x (energy), row 1 is y (yield)
            let yield_data: Vec<f64> = yield_ds.read_raw::<f64>()?;
            let n = yield_shape[1];
            let x: Vec<f64> = yield_data.iter().cloned().take(n).collect();
            let y: Vec<f64> = yield_data.iter().cloned().skip(n).take(n).collect();
            Some(Yield::Tabulated { x, y })
        } else {
            None
        };

        // Read distributions
        let n_distribution: i32 = prod_group.attr("n_distribution")?.read_scalar()?;
        let mut distributions = Vec::new();

        for i in 0..n_distribution {
            let dist_name = format!("distribution_{}", i);
            if let Ok(dist_group) = prod_group.group(&dist_name) {
                match read_distribution(&dist_group) {
                    Ok(dist) => distributions.push(dist),
                    Err(e) => {
                        // Log but continue - some distributions may be unsupported
                        eprintln!("Warning: Failed to read {}: {}", dist_name, e);
                    }
                }
            }
        }

        let product = ReactionProduct {
            particle,
            emission_mode,
            decay_rate: 0.0,
            applicability: vec![],
            distribution: distributions,
            product_yield,
        };
        products.push(product);
    }

    Ok(products)
}

/// Read an angle-energy distribution from HDF5
fn read_distribution(group: &Group) -> Result<AngleEnergyDistribution, Box<dyn std::error::Error>> {
    let dist_type = read_string_attr(group, "type")?;

    // Debug: log what distribution type is being read
    if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
        eprintln!("  read_distribution: type={}", dist_type);
    }

    match dist_type.as_str() {
        "uncorrelated" => read_uncorrelated_distribution(group),
        "correlated" => {
            let result = read_correlated_distribution(group);
            if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
                eprintln!("  read_correlated_distribution: ok={}", result.is_ok());
            }
            result
        },
        "kalbach-mann" => read_kalbach_mann_distribution(group),
        _ => {
            // Default to uncorrelated with isotropic angle
            if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
                eprintln!("  read_distribution: unknown type '{}', falling back to uncorrelated", dist_type);
            }
            Ok(AngleEnergyDistribution::UncorrelatedAngleEnergy {
                angle: AngleDistribution {
                    energy: vec![],
                    mu: vec![],
                },
                energy: None,
            })
        }
    }
}

/// Read uncorrelated angle-energy distribution
fn read_uncorrelated_distribution(group: &Group) -> Result<AngleEnergyDistribution, Box<dyn std::error::Error>> {
    // Check if this is an evaporation energy distribution first
    // Evaporation is a special case that returns its own AngleEnergyDistribution variant
    if let Ok(energy_group) = group.group("energy") {
        let energy_type = read_string_attr(&energy_group, "type").unwrap_or_default();
        if energy_type == "evaporation" {
            // Read evaporation parameters
            // theta is the nuclear temperature parameter (Tabulated1D)
            // u is the restriction energy (threshold)
            let u = read_f64_attr(&energy_group, "u").unwrap_or(0.0);

            // Read theta as Tabulated1D
            let theta = if let Ok(theta_ds) = energy_group.dataset("theta") {
                // theta is a Tabulated1D with x (energy), y (temperature), breakpoints, interpolation
                // OpenMC stores it as a group, but some files may have it as a dataset
                // For now, try to read it as a simple Tabulated1D structure
                if let Ok(theta_group) = energy_group.group("theta") {
                    let x: Vec<f64> = theta_group.dataset("x")?.read_1d()?.to_vec();
                    let y: Vec<f64> = theta_group.dataset("y")?.read_1d()?.to_vec();
                    let breakpoints: Vec<i32> = theta_group.dataset("breakpoints")
                        .and_then(|ds| ds.read_1d().map(|arr| arr.to_vec()))
                        .unwrap_or_default();
                    let interpolation: Vec<i32> = theta_group.dataset("interpolation")
                        .and_then(|ds| ds.read_1d().map(|arr| arr.to_vec()))
                        .unwrap_or_default();
                    Some(Tabulated1D::Tabulated1D { x, y, breakpoints, interpolation })
                } else {
                    // theta might be a simple 2-column dataset
                    if let Ok(theta_data) = theta_ds.read_2d::<f64>() {
                        let x: Vec<f64> = theta_data.row(0).to_vec();
                        let y: Vec<f64> = theta_data.row(1).to_vec();
                        Some(Tabulated1D::Tabulated1D {
                            x, y,
                            breakpoints: vec![],
                            interpolation: vec![]
                        })
                    } else {
                        None
                    }
                }
            } else if let Ok(theta_group) = energy_group.group("theta") {
                // theta is stored as a group
                let x: Vec<f64> = theta_group.dataset("x")?.read_1d()?.to_vec();
                let y: Vec<f64> = theta_group.dataset("y")?.read_1d()?.to_vec();
                let breakpoints: Vec<i32> = theta_group.dataset("breakpoints")
                    .and_then(|ds| ds.read_1d().map(|arr| arr.to_vec()))
                    .unwrap_or_default();
                let interpolation: Vec<i32> = theta_group.dataset("interpolation")
                    .and_then(|ds| ds.read_1d().map(|arr| arr.to_vec()))
                    .unwrap_or_default();
                Some(Tabulated1D::Tabulated1D { x, y, breakpoints, interpolation })
            } else {
                None
            };

            if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
                eprintln!("  read_uncorrelated: evaporation u={:.2e}, theta={}", u, theta.is_some());
            }

            return Ok(AngleEnergyDistribution::Evaporation { theta, u });
        }
    }

    let mut angle = AngleDistribution {
        energy: vec![],
        mu: vec![],
    };

    // Read angle distribution if present
    if let Ok(angle_group) = group.group("angle") {
        if let Ok(energy_ds) = angle_group.dataset("energy") {
            angle.energy = energy_ds.read_1d()?.to_vec();
        }
        if let Ok(mu_ds) = angle_group.dataset("mu") {
            // mu is shape (3, N) where:
            // - row 0: mu values (cosine of angle)
            // - row 1: PDF values
            // - row 2: CDF values
            // The data is offset-encoded for each incoming energy
            let mu_data: Vec<f64> = mu_ds.read_raw::<f64>()?;

            // Read offsets if available
            if let Ok(offsets_attr) = mu_ds.attr("offsets") {
                let offsets: Vec<i32> = offsets_attr.read_1d()?.to_vec();
                let n_energy = angle.energy.len();
                let n_cols = mu_data.len() / 3;

                for i in 0..n_energy {
                    let j = offsets[i] as usize;
                    let n = if i < n_energy - 1 {
                        (offsets[i + 1] - offsets[i]) as usize
                    } else {
                        n_cols - j
                    };

                    let x: Vec<f64> = (0..n).map(|k| mu_data[0 * n_cols + j + k]).collect();
                    let p: Vec<f64> = (0..n).map(|k| mu_data[1 * n_cols + j + k]).collect();

                    angle.mu.push(crate::reaction_product::Tabulated { x, p });
                }
            } else {
                // No offsets - assume single distribution or isotropic
                let ncols = mu_data.len() / 3;
                let x: Vec<f64> = mu_data[0..ncols].to_vec();
                let p: Vec<f64> = mu_data[ncols..2*ncols].to_vec();
                angle.mu.push(crate::reaction_product::Tabulated { x, p });
            }
        }
    }

    // Read energy distribution if present
    let energy_dist = if let Ok(energy_group) = group.group("energy") {
        let energy_type = read_string_attr(&energy_group, "type").unwrap_or_default();
        match energy_type.as_str() {
            "level" => {
                // Read threshold and mass_ratio attributes (OpenMC LevelInelastic format)
                let threshold: f64 = energy_group.attr("threshold")
                    .and_then(|a| a.read_scalar())
                    .unwrap_or(0.0);
                let mass_ratio: f64 = energy_group.attr("mass_ratio")
                    .and_then(|a| a.read_scalar())
                    .unwrap_or(1.0);

                if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
                    eprintln!("  level: threshold={:.4e}, mass_ratio={:.4}", threshold, mass_ratio);
                }

                Some(EnergyDistribution::LevelInelastic { threshold, mass_ratio })
            }
            "continuous" => {
                // Read continuous tabular distribution
                // The energy group contains 'energy' dataset (incoming energies)
                // and 'distribution' dataset (outgoing energy PDFs with shape (3, N))
                let energy_arr: Vec<f64> = match energy_group.dataset("energy") {
                    Ok(ds) => ds.read_1d::<f64>().map(|a| a.to_vec()).unwrap_or_default(),
                    Err(_) => vec![],
                };

                if energy_arr.is_empty() {
                    None
                } else {
                    // Read distribution dataset: shape (3, N) with e_out, p, c
                    let energy_out = if let Ok(dist_ds) = energy_group.dataset("distribution") {
                        // Get offsets attribute - tells us where each incoming energy's distribution starts
                        let offsets: Vec<i32> = dist_ds.attr("offsets")
                            .and_then(|a| a.read_1d().map(|arr| arr.to_vec()))
                            .unwrap_or_default();

                        // Read the distribution data
                        let dist_data: Vec<f64> = dist_ds.read_raw::<f64>().unwrap_or_default();
                        let n_cols = dist_data.len() / 3;

                        if std::env::var("YAMC_DEBUG_LOAD").is_ok() {
                            eprintln!("  continuous: {} energies, dist shape (3, {}), {} offsets",
                                energy_arr.len(), n_cols, offsets.len());
                        }

                        // Build TabulatedProbability for each incoming energy
                        let n_energy = energy_arr.len();
                        let mut energy_out_vec = Vec::with_capacity(n_energy);

                        for i in 0..n_energy {
                            let j = offsets.get(i).copied().unwrap_or(0) as usize;
                            let n = if i < n_energy - 1 {
                                (offsets.get(i + 1).copied().unwrap_or(0) - offsets[i]) as usize
                            } else {
                                n_cols.saturating_sub(j)
                            };

                            if n > 0 && j + n <= n_cols {
                                // Row 0: e_out, Row 1: p (PDF), Row 2: c (CDF - we use the PDF and convert)
                                let e_out: Vec<f64> = (0..n).map(|k| dist_data[0 * n_cols + j + k]).collect();
                                let p: Vec<f64> = (0..n).map(|k| dist_data[1 * n_cols + j + k]).collect();

                                energy_out_vec.push(crate::reaction_product::TabulatedProbability::Tabulated {
                                    x: e_out,
                                    p,
                                });
                            }
                        }
                        energy_out_vec
                    } else {
                        vec![]
                    };

                    Some(EnergyDistribution::ContinuousTabular {
                        energy: energy_arr,
                        energy_out,
                    })
                }
            }
            _ => None,
        }
    } else {
        None
    };

    Ok(AngleEnergyDistribution::UncorrelatedAngleEnergy {
        angle,
        energy: energy_dist,
    })
}

/// Read correlated angle-energy distribution
fn read_correlated_distribution(group: &Group) -> Result<AngleEnergyDistribution, Box<dyn std::error::Error>> {
    // Read incoming energy grid
    let energy_ds = group.dataset("energy")?;
    let energy: Vec<f64> = energy_ds.read_1d()?.to_vec();
    let n_energy = energy.len();

    // Note: energy dataset has an 'interpolation' attribute that is shape (2, N)
    // where row 0 is breakpoints and row 1 is interpolation types
    // We don't need it for basic loading but could use it for more accurate interpolation

    // Read outgoing energy distribution data
    // Shape is (5, N) where:
    // Row 0: e_out (outgoing energy)
    // Row 1: p (PDF)
    // Row 2: c (CDF)
    // Row 3: interp_mu (angular interpolation scheme for each e_out)
    // Row 4: offset_mu (offset into mu array for each e_out)
    let eout_ds = group.dataset("energy_out")?;
    let eout_data: Vec<f64> = eout_ds.read_raw::<f64>()?;

    // Read attributes
    let offsets: Vec<i32> = eout_ds.attr("offsets")?.read_1d()?.to_vec();
    let interp: Vec<i32> = eout_ds.attr("interpolation")?.read_1d()?.to_vec();
    let n_discrete: Vec<i32> = eout_ds.attr("n_discrete_lines")?.read_1d()?.to_vec();

    // Read mu data for angular distributions
    // Shape is (3, M) where:
    // Row 0: mu values (cosine of angle)
    // Row 1: p (PDF)
    // Row 2: c (CDF)
    let mu_data: Vec<f64> = group.dataset("mu")?.read_raw::<f64>()?;
    let mu_ncols = mu_data.len() / 3;
    let eout_ncols = eout_data.len() / 5;

    // Build distributions for each incoming energy
    let mut distributions = Vec::with_capacity(n_energy);

    for i in 0..n_energy {
        let j = offsets[i] as usize;
        let n = if i < n_energy - 1 {
            (offsets[i + 1] - offsets[i]) as usize
        } else {
            eout_ncols - j
        };

        let interpolation = if interp[i] == 1 {
            CorrInterpolation::Histogram
        } else {
            CorrInterpolation::LinLin
        };

        let e_out: Vec<f64> = (0..n).map(|k| eout_data[0 * eout_ncols + j + k]).collect();
        let p: Vec<f64> = (0..n).map(|k| eout_data[1 * eout_ncols + j + k]).collect();
        let c: Vec<f64> = (0..n).map(|k| eout_data[2 * eout_ncols + j + k]).collect();

        // Read angular distributions for each outgoing energy
        let mut angle: Vec<Tabular> = Vec::with_capacity(n);
        for k in 0..n {
            let idx = j + k;
            let interp_mu = eout_data[3 * eout_ncols + idx] as i32;
            let offset_mu = eout_data[4 * eout_ncols + idx] as usize;

            // Determine size of angular distribution
            let m = if idx + 1 < eout_ncols {
                let next_offset = eout_data[4 * eout_ncols + idx + 1] as usize;
                next_offset - offset_mu
            } else {
                mu_ncols - offset_mu
            };

            // Read mu data for this outgoing energy
            let mu_x: Vec<f64> = (0..m).map(|l| mu_data[0 * mu_ncols + offset_mu + l]).collect();
            let mu_p: Vec<f64> = (0..m).map(|l| mu_data[1 * mu_ncols + offset_mu + l]).collect();
            let mu_c: Vec<f64> = (0..m).map(|l| mu_data[2 * mu_ncols + offset_mu + l]).collect();

            // Angular interpolation: 0 or 1 = histogram, 2 = lin-lin
            let mu_interp = if interp_mu == 0 || interp_mu == 1 {
                CorrInterpolation::Histogram
            } else {
                CorrInterpolation::LinLin
            };

            angle.push(Tabular {
                x: mu_x,
                p: mu_p,
                c: mu_c,
                interpolation: mu_interp,
                n_discrete: 0, // Angular distributions typically have no discrete lines
            });
        }

        let table = CorrTable {
            interpolation,
            n_discrete: n_discrete[i] as usize,
            e_out,
            p,
            c,
            angle,
        };
        distributions.push(table);
    }

    let correlated = CorrelatedAngleEnergy {
        energy,
        distributions,
    };

    Ok(AngleEnergyDistribution::CorrelatedAngleEnergy {
        correlated,
    })
}

/// Read Kalbach-Mann distribution
fn read_kalbach_mann_distribution(group: &Group) -> Result<AngleEnergyDistribution, Box<dyn std::error::Error>> {
    // Read incoming energy grid
    let energy_ds = group.dataset("energy")?;
    let energy: Vec<f64> = energy_ds.read_1d()?.to_vec();
    let n_energy = energy.len();

    // Read distribution data
    let dist_ds = group.dataset("distribution")?;
    let dist_data: Vec<f64> = dist_ds.read_raw::<f64>()?;

    // Read attributes
    let offsets: Vec<i32> = dist_ds.attr("offsets")?.read_1d()?.to_vec();
    let interp: Vec<i32> = dist_ds.attr("interpolation")?.read_1d()?.to_vec();
    let n_discrete: Vec<i32> = dist_ds.attr("n_discrete_lines")?.read_1d()?.to_vec();
    let dist_ncols = dist_data.len() / 5;

    // Build distributions for each incoming energy
    let mut distributions = Vec::with_capacity(n_energy);

    for i in 0..n_energy {
        let j = offsets[i] as usize;
        let n = if i < n_energy - 1 {
            (offsets[i + 1] - offsets[i]) as usize
        } else {
            dist_ncols - j
        };

        let interpolation = if interp[i] == 1 {
            KMInterpolation::Histogram
        } else {
            KMInterpolation::LinLin
        };

        let e_out: Vec<f64> = (0..n).map(|k| dist_data[0 * dist_ncols + j + k]).collect();
        let p: Vec<f64> = (0..n).map(|k| dist_data[1 * dist_ncols + j + k]).collect();
        let c: Vec<f64> = (0..n).map(|k| dist_data[2 * dist_ncols + j + k]).collect();
        let r: Vec<f64> = (0..n).map(|k| dist_data[3 * dist_ncols + j + k]).collect();
        let a: Vec<f64> = (0..n).map(|k| dist_data[4 * dist_ncols + j + k]).collect();

        let table = KMTable {
            interpolation,
            n_discrete: n_discrete[i] as usize,
            e_out,
            p,
            c,
            r,
            a,
        };
        distributions.push(table);
    }

    let kalbach = KalbachMann {
        energy,
        distributions,
    };

    Ok(AngleEnergyDistribution::KalbachMann {
        kalbach,
    })
}

/// Check if an MT number corresponds to fission
fn is_fission(mt: i32) -> bool {
    mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_be9_hdf5() {
        let path = "tests/Be9.h5";
        let nuclide = match read_nuclide_from_hdf5(path, None) {
            Ok(n) => n,
            Err(e) => panic!("Failed to read Be9.h5: {:?}", e),
        };

        assert_eq!(nuclide.name, Some("Be9".to_string()));
        assert_eq!(nuclide.atomic_number, Some(4));
        assert_eq!(nuclide.mass_number, Some(9));
        assert!(!nuclide.reactions.is_empty());

        // Check that we have reactions for 294K
        assert!(nuclide.reactions.contains_key("294"));

        let reactions_294 = &nuclide.reactions["294"];
        assert!(reactions_294.contains_key(&2)); // elastic

        println!("Loaded {} temperatures", nuclide.loaded_temperatures.len());
        println!("Available temperatures: {:?}", nuclide.available_temperatures);
        println!("Reactions at 294K: {:?}", reactions_294.keys().collect::<Vec<_>>());
    }
}
