use once_cell::sync::Lazy;
use std::collections::HashMap;

/// Map from element symbol to sorted vector of nuclide (isotope) identifiers.
///
/// Keys are element symbols (e.g. `"Li"`) and values are the sorted list of
/// nuclide names for naturally occurring isotopes of that element (e.g.
/// `["Li6", "Li7"]`). The mapping is derived automatically from
/// [`NATURAL_ABUNDANCE`] so it stays consistent with the set of isotopes for
/// which natural abundances are defined.
pub static ELEMENT_NUCLIDES: Lazy<HashMap<&'static str, Vec<&'static str>>> = Lazy::new(|| {
    let mut map: HashMap<&'static str, Vec<&'static str>> = HashMap::new();
    for &nuclide in NATURAL_ABUNDANCE.keys() {
        // Find the index where the first digit occurs
        let idx = nuclide
            .find(|c: char| c.is_ascii_digit())
            .unwrap_or(nuclide.len());
        let element = &nuclide[..idx]; // This is a &'static str slice
        map.entry(element).or_insert_with(Vec::new).push(nuclide);
    }
    // Sort nuclides for each element
    for nuclides in map.values_mut() {
        nuclides.sort();
    }
    map
});
// src/data.rs
// This module contains large static data tables for the materials library.
// The volume of data is significant; doc comments summarize the intent of
// each table while the literals provide the canonical numeric values.

/// Natural terrestrial isotopic abundances (fractional, summing to ~1.0 per
/// element) for stable isotopes.
///
/// Each key is a nuclide name (e.g. `"Fe56"`) and the value is its natural
/// abundance by atom fraction. Values are sourced from standard reference
/// compilations (rounded as needed). Elements with a single stable isotope are
/// assigned 1.0.
pub static NATURAL_ABUNDANCE: Lazy<HashMap<&'static str, f64>> = Lazy::new(|| {
    let mut m = HashMap::new();
    // ...existing NATURAL_ABUNDANCE data from element.rs...
    // Example entries (replace with full data):
    m.insert("Li6", 0.0759);
    m.insert("Li7", 0.9241);

    // Hydrogen
    m.insert("H1", 0.99984426);
    m.insert("H2", 0.00015574);

    // Helium
    m.insert("He3", 0.000002);
    m.insert("He4", 0.999998);

    // Lithium
    m.insert("Li6", 0.07589);
    m.insert("Li7", 0.92411);

    // Beryllium
    m.insert("Be9", 1.0);

    // Boron
    m.insert("B10", 0.1982);
    m.insert("B11", 0.8018);

    // Carbon
    m.insert("C12", 0.988922);
    m.insert("C13", 0.011078);

    // Nitrogen
    m.insert("N14", 0.996337);
    m.insert("N15", 0.003663);

    // Oxygen
    m.insert("O16", 0.9976206);
    m.insert("O17", 0.000379);
    m.insert("O18", 0.0020004);

    // Fluorine
    m.insert("F19", 1.0);

    // Neon
    m.insert("Ne20", 0.9048);
    m.insert("Ne21", 0.0027);
    m.insert("Ne22", 0.0925);

    // Sodium
    m.insert("Na23", 1.0);

    // Magnesium
    m.insert("Mg24", 0.78951);
    m.insert("Mg25", 0.1002);
    m.insert("Mg26", 0.11029);

    // Aluminum
    m.insert("Al27", 1.0);

    // Silicon
    m.insert("Si28", 0.9222968);
    m.insert("Si29", 0.0468316);
    m.insert("Si30", 0.0308716);

    // Phosphorus
    m.insert("P31", 1.0);

    // Sulfur
    m.insert("S32", 0.9504074);
    m.insert("S33", 0.0074869);
    m.insert("S34", 0.0419599);
    m.insert("S36", 0.0001458);

    // Chlorine
    m.insert("Cl35", 0.757647);
    m.insert("Cl37", 0.242353);

    // Argon
    m.insert("Ar36", 0.003336);
    m.insert("Ar38", 0.000629);
    m.insert("Ar40", 0.996035);

    // Potassium
    m.insert("K39", 0.932581);
    m.insert("K40", 0.000117);
    m.insert("K41", 0.067302);

    // Calcium
    m.insert("Ca40", 0.96941);
    m.insert("Ca42", 0.00647);
    m.insert("Ca43", 0.00135);
    m.insert("Ca44", 0.02086);
    m.insert("Ca46", 0.00004);
    m.insert("Ca48", 0.00187);

    // Scandium
    m.insert("Sc45", 1.0);

    // Titanium
    m.insert("Ti46", 0.0825);
    m.insert("Ti47", 0.0744);
    m.insert("Ti48", 0.7372);
    m.insert("Ti49", 0.0541);
    m.insert("Ti50", 0.0518);

    // Vanadium
    m.insert("V50", 0.0025);
    m.insert("V51", 0.9975);

    // Chromium
    m.insert("Cr50", 0.04345);
    m.insert("Cr52", 0.83789);
    m.insert("Cr53", 0.09501);
    m.insert("Cr54", 0.02365);

    // Manganese
    m.insert("Mn55", 1.0);

    // Iron
    m.insert("Fe54", 0.05845);
    m.insert("Fe56", 0.91754);
    m.insert("Fe57", 0.02119);
    m.insert("Fe58", 0.00282);

    // Cobalt
    m.insert("Co59", 1.0);

    // Nickel
    m.insert("Ni58", 0.680769);
    m.insert("Ni60", 0.262231);
    m.insert("Ni61", 0.011399);
    m.insert("Ni62", 0.036345);
    m.insert("Ni64", 0.009256);

    // Copper
    m.insert("Cu63", 0.6915);
    m.insert("Cu65", 0.3085);

    // Zinc
    m.insert("Zn64", 0.4917);
    m.insert("Zn66", 0.2773);
    m.insert("Zn67", 0.0404);
    m.insert("Zn68", 0.1845);
    m.insert("Zn70", 0.0061);

    // Gallium
    m.insert("Ga69", 0.60108);
    m.insert("Ga71", 0.39892);

    // Germanium
    m.insert("Ge70", 0.2052);
    m.insert("Ge72", 0.2745);
    m.insert("Ge73", 0.0776);
    m.insert("Ge74", 0.3652);
    m.insert("Ge76", 0.0775);

    // Arsenic
    m.insert("As75", 1.0);

    // Selenium
    m.insert("Se74", 0.0086);
    m.insert("Se76", 0.0923);
    m.insert("Se77", 0.076);
    m.insert("Se78", 0.2369);
    m.insert("Se80", 0.498);
    m.insert("Se82", 0.0882);

    // Bromine
    m.insert("Br79", 0.50686);
    m.insert("Br81", 0.49314);

    // Krypton
    m.insert("Kr78", 0.00355);
    m.insert("Kr80", 0.02286);
    m.insert("Kr82", 0.11593);
    m.insert("Kr83", 0.115);
    m.insert("Kr84", 0.56987);
    m.insert("Kr86", 0.17279);

    // Rubidium
    m.insert("Rb85", 0.7217);
    m.insert("Rb87", 0.2783);

    // Strontium
    m.insert("Sr84", 0.0056);
    m.insert("Sr86", 0.0986);
    m.insert("Sr87", 0.07);
    m.insert("Sr88", 0.8258);

    // Yttrium
    m.insert("Y89", 1.0);

    // Zirconium
    m.insert("Zr90", 0.5145);
    m.insert("Zr91", 0.1122);
    m.insert("Zr92", 0.1715);
    m.insert("Zr94", 0.1738);
    m.insert("Zr96", 0.028);

    // Niobium
    m.insert("Nb93", 1.0);

    // Molybdenum
    m.insert("Mo92", 0.14649);
    m.insert("Mo94", 0.09187);
    m.insert("Mo95", 0.15873);
    m.insert("Mo96", 0.16673);
    m.insert("Mo97", 0.09582);
    m.insert("Mo98", 0.24292);
    m.insert("Mo100", 0.09744);

    // Ruthenium
    m.insert("Ru96", 0.0554);
    m.insert("Ru98", 0.0187);
    m.insert("Ru99", 0.1276);
    m.insert("Ru100", 0.126);
    m.insert("Ru101", 0.1706);
    m.insert("Ru102", 0.3155);
    m.insert("Ru104", 0.1862);

    // Rhodium
    m.insert("Rh103", 1.0);

    // Palladium
    m.insert("Pd102", 0.0102);
    m.insert("Pd104", 0.1114);
    m.insert("Pd105", 0.2233);
    m.insert("Pd106", 0.2733);
    m.insert("Pd108", 0.2646);
    m.insert("Pd110", 0.1172);

    // Silver
    m.insert("Ag107", 0.51839);
    m.insert("Ag109", 0.48161);

    // Cadmium
    m.insert("Cd106", 0.01245);
    m.insert("Cd108", 0.00888);
    m.insert("Cd110", 0.1247);
    m.insert("Cd111", 0.12795);
    m.insert("Cd112", 0.24109);
    m.insert("Cd113", 0.12227);
    m.insert("Cd114", 0.28754);
    m.insert("Cd116", 0.07512);

    // Indium
    m.insert("In113", 0.04281);
    m.insert("In115", 0.95719);

    // Tin
    m.insert("Sn112", 0.0097);
    m.insert("Sn114", 0.0066);
    m.insert("Sn115", 0.0034);
    m.insert("Sn116", 0.1454);
    m.insert("Sn117", 0.0768);
    m.insert("Sn118", 0.2422);
    m.insert("Sn119", 0.0859);
    m.insert("Sn120", 0.3258);
    m.insert("Sn122", 0.0463);
    m.insert("Sn124", 0.0579);

    // Antimony
    m.insert("Sb121", 0.5721);
    m.insert("Sb123", 0.4279);

    // Tellurium
    m.insert("Te120", 0.0009);
    m.insert("Te122", 0.0255);
    m.insert("Te123", 0.0089);
    m.insert("Te124", 0.0474);
    m.insert("Te125", 0.0707);
    m.insert("Te126", 0.1884);
    m.insert("Te128", 0.3174);
    m.insert("Te130", 0.3408);

    // Iodine
    m.insert("I127", 1.0);

    // Xenon
    m.insert("Xe124", 0.00095);
    m.insert("Xe126", 0.00089);
    m.insert("Xe128", 0.0191);
    m.insert("Xe129", 0.26401);
    m.insert("Xe130", 0.04071);
    m.insert("Xe131", 0.21232);
    m.insert("Xe132", 0.26909);
    m.insert("Xe134", 0.10436);
    m.insert("Xe136", 0.08857);

    // Cesium
    m.insert("Cs133", 1.0);

    // Barium
    m.insert("Ba130", 0.0011);
    m.insert("Ba132", 0.001);
    m.insert("Ba134", 0.0242);
    m.insert("Ba135", 0.0659);
    m.insert("Ba136", 0.0785);
    m.insert("Ba137", 0.1123);
    m.insert("Ba138", 0.717);

    // Lanthanum
    m.insert("La138", 0.0008881);
    m.insert("La139", 0.9991119);

    // Cerium
    m.insert("Ce136", 0.00186);
    m.insert("Ce138", 0.00251);
    m.insert("Ce140", 0.88449);
    m.insert("Ce142", 0.11114);

    // Praseodymium
    m.insert("Pr141", 1.0);

    // Neodymium
    m.insert("Nd142", 0.27153);
    m.insert("Nd143", 0.12173);
    m.insert("Nd144", 0.23798);
    m.insert("Nd145", 0.08293);
    m.insert("Nd146", 0.17189);
    m.insert("Nd148", 0.05756);
    m.insert("Nd150", 0.05638);

    // Samarium
    m.insert("Sm144", 0.0308);
    m.insert("Sm147", 0.15);
    m.insert("Sm148", 0.1125);
    m.insert("Sm149", 0.1382);
    m.insert("Sm150", 0.0737);
    m.insert("Sm152", 0.2674);
    m.insert("Sm154", 0.2274);

    // Europium
    m.insert("Eu151", 0.4781);
    m.insert("Eu153", 0.5219);

    // Gadolinium
    m.insert("Gd152", 0.002);
    m.insert("Gd154", 0.0218);
    m.insert("Gd155", 0.148);
    m.insert("Gd156", 0.2047);
    m.insert("Gd157", 0.1565);
    m.insert("Gd158", 0.2484);
    m.insert("Gd160", 0.2186);

    // Terbium
    m.insert("Tb159", 1.0);

    // Dysprosium
    m.insert("Dy156", 0.00056);
    m.insert("Dy158", 0.00095);
    m.insert("Dy160", 0.02329);
    m.insert("Dy161", 0.18889);
    m.insert("Dy162", 0.25475);
    m.insert("Dy163", 0.24896);
    m.insert("Dy164", 0.2826);

    // Holmium
    m.insert("Ho165", 1.0);

    // Erbium
    m.insert("Er162", 0.00139);
    m.insert("Er164", 0.01601);
    m.insert("Er166", 0.33503);
    m.insert("Er167", 0.22869);
    m.insert("Er168", 0.26978);
    m.insert("Er170", 0.1491);

    // Thulium
    m.insert("Tm169", 1.0);

    // Ytterbium
    m.insert("Yb168", 0.00123);
    m.insert("Yb170", 0.02982);
    m.insert("Yb171", 0.14086);
    m.insert("Yb172", 0.21686);
    m.insert("Yb173", 0.16103);
    m.insert("Yb174", 0.32025);
    m.insert("Yb176", 0.12995);

    // Lutetium
    m.insert("Lu175", 0.97401);
    m.insert("Lu176", 0.02599);

    // Hafnium
    m.insert("Hf174", 0.0016);
    m.insert("Hf176", 0.0526);
    m.insert("Hf177", 0.186);
    m.insert("Hf178", 0.2728);
    m.insert("Hf179", 0.1362);
    m.insert("Hf180", 0.3508);

    // Tantalum
    m.insert("Ta180_m1", 0.0001201);
    m.insert("Ta181", 0.9998799);

    // Tungsten
    m.insert("W180", 0.0012);
    m.insert("W182", 0.265);
    m.insert("W183", 0.1431);
    m.insert("W184", 0.3064);
    m.insert("W186", 0.2843);

    // Rhenium
    m.insert("Re185", 0.374);
    m.insert("Re187", 0.626);

    // Osmium
    m.insert("Os184", 0.0002);
    m.insert("Os186", 0.0159);
    m.insert("Os187", 0.0196);
    m.insert("Os188", 0.1324);
    m.insert("Os189", 0.1615);
    m.insert("Os190", 0.2626);
    m.insert("Os192", 0.4078);

    // Iridium
    m.insert("Ir191", 0.373);
    m.insert("Ir193", 0.627);

    // Platinum
    m.insert("Pt190", 0.00012);
    m.insert("Pt192", 0.00782);
    m.insert("Pt194", 0.32864);
    m.insert("Pt195", 0.33775);
    m.insert("Pt196", 0.25211);
    m.insert("Pt198", 0.07356);

    // Gold
    m.insert("Au197", 1.0);

    // Mercury
    m.insert("Hg196", 0.0015);
    m.insert("Hg198", 0.1004);
    m.insert("Hg199", 0.1694);
    m.insert("Hg200", 0.2314);
    m.insert("Hg201", 0.1317);
    m.insert("Hg202", 0.2974);
    m.insert("Hg204", 0.0682);

    // Thallium
    m.insert("Tl203", 0.29524);
    m.insert("Tl205", 0.70476);

    // Lead
    m.insert("Pb204", 0.014);
    m.insert("Pb206", 0.241);
    m.insert("Pb207", 0.221);
    m.insert("Pb208", 0.524);

    // Bismuth
    m.insert("Bi209", 1.0);

    // Thorium
    m.insert("Th230", 0.0002);
    m.insert("Th232", 0.9998);

    // Protactinium
    m.insert("Pa231", 1.0);

    // Uranium
    m.insert("U234", 0.000054);
    m.insert("U235", 0.007204);
    m.insert("U238", 0.992742);
    m
});

/// Mapping from element symbol to its lowercase English name.
///
/// Provided for convenience when presenting user‑facing descriptions and for
/// validating element inputs (case sensitive symbol keys matching the raw
/// nuclear data tables).
pub static ELEMENT_NAMES: Lazy<HashMap<&'static str, &'static str>> = Lazy::new(|| {
    let mut names = HashMap::new();
    names.insert("H", "hydrogen");
    names.insert("He", "helium");
    names.insert("Li", "lithium");
    names.insert("Be", "beryllium");
    names.insert("B", "boron");
    names.insert("C", "carbon");
    names.insert("N", "nitrogen");
    names.insert("O", "oxygen");
    names.insert("F", "fluorine");
    names.insert("Ne", "neon");
    names.insert("Na", "sodium");
    names.insert("Mg", "magnesium");
    names.insert("Al", "aluminum");
    names.insert("Si", "silicon");
    names.insert("P", "phosphorus");
    names.insert("S", "sulfur");
    names.insert("Cl", "chlorine");
    names.insert("Ar", "argon");
    names.insert("K", "potassium");
    names.insert("Ca", "calcium");
    names.insert("Sc", "scandium");
    names.insert("Ti", "titanium");
    names.insert("V", "vanadium");
    names.insert("Cr", "chromium");
    names.insert("Mn", "manganese");
    names.insert("Fe", "iron");
    names.insert("Co", "cobalt");
    names.insert("Ni", "nickel");
    names.insert("Cu", "copper");
    names.insert("Zn", "zinc");
    names.insert("Ga", "gallium");
    names.insert("Ge", "germanium");
    names.insert("As", "arsenic");
    names.insert("Se", "selenium");
    names.insert("Br", "bromine");
    names.insert("Kr", "krypton");
    names.insert("Rb", "rubidium");
    names.insert("Sr", "strontium");
    names.insert("Y", "yttrium");
    names.insert("Zr", "zirconium");
    names.insert("Nb", "niobium");
    names.insert("Mo", "molybdenum");
    names.insert("Tc", "technetium");
    names.insert("Ru", "ruthenium");
    names.insert("Rh", "rhodium");
    names.insert("Pd", "palladium");
    names.insert("Ag", "silver");
    names.insert("Cd", "cadmium");
    names.insert("In", "indium");
    names.insert("Sn", "tin");
    names.insert("Sb", "antimony");
    names.insert("Te", "tellurium");
    names.insert("I", "iodine");
    names.insert("Xe", "xenon");
    names.insert("Cs", "cesium");
    names.insert("Ba", "barium");
    names.insert("La", "lanthanum");
    names.insert("Ce", "cerium");
    names.insert("Pr", "praseodymium");
    names.insert("Nd", "neodymium");
    names.insert("Pm", "promethium");
    names.insert("Sm", "samarium");
    names.insert("Eu", "europium");
    names.insert("Gd", "gadolinium");
    names.insert("Tb", "terbium");
    names.insert("Dy", "dysprosium");
    names.insert("Ho", "holmium");
    names.insert("Er", "erbium");
    names.insert("Tm", "thulium");
    names.insert("Yb", "ytterbium");
    names.insert("Lu", "lutetium");
    names.insert("Hf", "hafnium");
    names.insert("Ta", "tantalum");
    names.insert("W", "tungsten");
    names.insert("Re", "rhenium");
    names.insert("Os", "osmium");
    names.insert("Ir", "iridium");
    names.insert("Pt", "platinum");
    names.insert("Au", "gold");
    names.insert("Hg", "mercury");
    names.insert("Tl", "thallium");
    names.insert("Pb", "lead");
    names.insert("Bi", "bismuth");
    names.insert("Po", "polonium");
    names.insert("At", "astatine");
    names.insert("Rn", "radon");
    names.insert("Fr", "francium");
    names.insert("Ra", "radium");
    names.insert("Ac", "actinium");
    names.insert("Th", "thorium");
    names.insert("Pa", "protactinium");
    names.insert("U", "uranium");
    names.insert("Np", "neptunium");
    names.insert("Pu", "plutonium");
    names.insert("Am", "americium");
    names.insert("Cm", "curium");
    names.insert("Bk", "berkelium");
    names.insert("Cf", "californium");
    names.insert("Es", "einsteinium");
    names.insert("Fm", "fermium");
    names.insert("Md", "mendelevium");
    names.insert("No", "nobelium");
    names.insert("Lr", "lawrencium");
    names
});


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_lithium_natural_abundance() {
        let li6 = NATURAL_ABUNDANCE.get("Li6").copied().unwrap_or(0.0);
        let li7 = NATURAL_ABUNDANCE.get("Li7").copied().unwrap_or(0.0);
        let sum = li6 + li7;
        assert!(
            (li6 - 0.0759).abs() < 1e-4,
            "Li6 abundance incorrect: {}",
            li6
        );
        assert!(
            (li7 - 0.9241).abs() < 1e-4,
            "Li7 abundance incorrect: {}",
            li7
        );
        assert!(
            (sum - 1.0).abs() < 1e-3,
            "Li6 + Li7 should sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_element_nuclides_li_and_be() {
        let li_nuclides = ELEMENT_NUCLIDES.get("Li").unwrap();
        assert_eq!(li_nuclides, &vec!["Li6", "Li7"]);
        let be_nuclides = ELEMENT_NUCLIDES.get("Be").unwrap();
        assert_eq!(be_nuclides, &vec!["Be9"]);
    }
}

/// A static HashMap that maps ENDF MT reaction numbers to their descriptive names
/// following the ENDF style naming convention.
#[allow(dead_code)]
pub static REACTION_NAME: Lazy<HashMap<i32, &'static str>> = Lazy::new(|| {
    [
        (1, "(n,total)"),
        (2, "(n,elastic)"),
        (3, "(n,nonelastic)"),
        (4, "(n,level)"),
        (5, "(n,misc)"),
        (11, "(n,2nd)"),
        (16, "(n,2n)"),
        (17, "(n,3n)"),
        (18, "(n,fission)"),
        (19, "(n,f)"),
        (20, "(n,nf)"),
        (21, "(n,2nf)"),
        (22, "(n,na)"),
        (23, "(n,n3a)"),
        (24, "(n,2na)"),
        (25, "(n,3na)"),
        (27, "(n,absorption)"),
        (28, "(n,np)"),
        (29, "(n,n2a)"),
        (30, "(n,2n2a)"),
        (32, "(n,nd)"),
        (33, "(n,nt)"),
        (34, "(n,n3He)"),
        (35, "(n,nd2a)"),
        (36, "(n,nt2a)"),
        (37, "(n,4n)"),
        (38, "(n,3nf)"),
        (41, "(n,2np)"),
        (42, "(n,3np)"),
        (44, "(n,n2p)"),
        (45, "(n,npa)"),
        (51, "(n,n1)"),
        (52, "(n,n2)"),
        (53, "(n,n3)"),
        (54, "(n,n4)"),
        (55, "(n,n5)"),
        (56, "(n,n6)"),
        (57, "(n,n7)"),
        (58, "(n,n8)"),
        (59, "(n,n9)"),
        (60, "(n,n10)"),
        (61, "(n,n11)"),
        (62, "(n,n12)"),
        (63, "(n,n13)"),
        (64, "(n,n14)"),
        (65, "(n,n15)"),
        (66, "(n,n16)"),
        (67, "(n,n17)"),
        (68, "(n,n18)"),
        (69, "(n,n19)"),
        (70, "(n,n20)"),
        (71, "(n,n21)"),
        (72, "(n,n22)"),
        (73, "(n,n23)"),
        (74, "(n,n24)"),
        (75, "(n,n25)"),
        (76, "(n,n26)"),
        (77, "(n,n27)"),
        (78, "(n,n28)"),
        (79, "(n,n29)"),
        (80, "(n,n30)"),
        (81, "(n,n31)"),
        (82, "(n,n32)"),
        (83, "(n,n33)"),
        (84, "(n,n34)"),
        (85, "(n,n35)"),
        (86, "(n,n36)"),
        (87, "(n,n37)"),
        (88, "(n,n38)"),
        (89, "(n,n39)"),
        (90, "(n,n40)"),
        (91, "(n,nc)"),
        (101, "(n,disappear)"),
        (102, "(n,gamma)"),
        (103, "(n,p)"),
        (104, "(n,d)"),
        (105, "(n,t)"),
        (106, "(n,3He)"),
        (107, "(n,a)"),
        (108, "(n,2a)"),
        (109, "(n,3a)"),
        (111, "(n,2p)"),
        (112, "(n,pa)"),
        (113, "(n,t2a)"),
        (114, "(n,d2a)"),
        (115, "(n,pd)"),
        (116, "(n,pt)"),
        (117, "(n,da)"),
        (152, "(n,5n)"),
        (153, "(n,6n)"),
        (154, "(n,2nt)"),
        (155, "(n,ta)"),
        (156, "(n,4np)"),
        (157, "(n,3nd)"),
        (158, "(n,nda)"),
        (159, "(n,2npa)"),
        (160, "(n,7n)"),
        (161, "(n,8n)"),
        (162, "(n,5np)"),
        (163, "(n,6np)"),
        (164, "(n,7np)"),
        (165, "(n,4na)"),
        (166, "(n,5na)"),
        (167, "(n,6na)"),
        (168, "(n,7na)"),
        (169, "(n,4nd)"),
        (170, "(n,5nd)"),
        (171, "(n,6nd)"),
        (172, "(n,3nt)"),
        (173, "(n,4nt)"),
        (174, "(n,5nt)"),
        (175, "(n,6nt)"),
        (176, "(n,2n3He)"),
        (177, "(n,3n3He)"),
        (178, "(n,4n3He)"),
        (179, "(n,3n2p)"),
        (180, "(n,3n2a)"),
        (181, "(n,3npa)"),
        (182, "(n,dt)"),
        (183, "(n,npd)"),
        (184, "(n,npt)"),
        (185, "(n,ndt)"),
        (186, "(n,np3He)"),
        (187, "(n,nd3He)"),
        (188, "(n,nt3He)"),
        (189, "(n,nta)"),
        (190, "(n,2n2p)"),
        (191, "(n,p3He)"),
        (192, "(n,d3He)"),
        (193, "(n,3Hea)"),
        (194, "(n,4n2p)"),
        (195, "(n,4n2a)"),
        (196, "(n,4npa)"),
        (197, "(n,3p)"),
        (198, "(n,n3p)"),
        (199, "(n,3n2pa)"),
        (200, "(n,5n2p)"),
        (203, "(n,Xp)"),
        (204, "(n,Xd)"),
        (205, "(n,Xt)"),
        (206, "(n,X3He)"),
        (207, "(n,Xa)"),
        (301, "heating"),
        (444, "damage-energy"),
        (600, "(n,p0)"),
        (601, "(n,p1)"),
        (602, "(n,p2)"),
        (603, "(n,p3)"),
        (604, "(n,p4)"),
        (605, "(n,p5)"),
        (606, "(n,p6)"),
        (607, "(n,p7)"),
        (608, "(n,p8)"),
        (609, "(n,p9)"),
        (610, "(n,p10)"),
        (611, "(n,p11)"),
        (612, "(n,p12)"),
        (613, "(n,p13)"),
        (614, "(n,p14)"),
        (615, "(n,p15)"),
        (616, "(n,p16)"),
        (617, "(n,p17)"),
        (618, "(n,p18)"),
        (619, "(n,p19)"),
        (620, "(n,p20)"),
        (621, "(n,p21)"),
        (622, "(n,p22)"),
        (623, "(n,p23)"),
        (624, "(n,p24)"),
        (625, "(n,p25)"),
        (626, "(n,p26)"),
        (627, "(n,p27)"),
        (628, "(n,p28)"),
        (629, "(n,p29)"),
        (630, "(n,p30)"),
        (631, "(n,p31)"),
        (632, "(n,p32)"),
        (633, "(n,p33)"),
        (634, "(n,p34)"),
        (635, "(n,p35)"),
        (636, "(n,p36)"),
        (637, "(n,p37)"),
        (638, "(n,p38)"),
        (639, "(n,p39)"),
        (640, "(n,p40)"),
        (641, "(n,p41)"),
        (642, "(n,p42)"),
        (643, "(n,p43)"),
        (644, "(n,p44)"),
        (645, "(n,p45)"),
        (646, "(n,p46)"),
        (647, "(n,p47)"),
        (648, "(n,p48)"),
        (649, "(n,pc)"),
        (650, "(n,d0)"),
        (651, "(n,d1)"),
        (652, "(n,d2)"),
        (653, "(n,d3)"),
        (654, "(n,d4)"),
        (655, "(n,d5)"),
        (656, "(n,d6)"),
        (657, "(n,d7)"),
        (658, "(n,d8)"),
        (659, "(n,d9)"),
        (660, "(n,d10)"),
        (661, "(n,d11)"),
        (662, "(n,d12)"),
        (663, "(n,d13)"),
        (664, "(n,d14)"),
        (665, "(n,d15)"),
        (666, "(n,d16)"),
        (667, "(n,d17)"),
        (668, "(n,d18)"),
        (669, "(n,d19)"),
        (670, "(n,d20)"),
        (671, "(n,d21)"),
        (672, "(n,d22)"),
        (673, "(n,d23)"),
        (674, "(n,d24)"),
        (675, "(n,d25)"),
        (676, "(n,d26)"),
        (677, "(n,d27)"),
        (678, "(n,d28)"),
        (679, "(n,d29)"),
        (680, "(n,d30)"),
        (681, "(n,d31)"),
        (682, "(n,d32)"),
        (683, "(n,d33)"),
        (684, "(n,d34)"),
        (685, "(n,d35)"),
        (686, "(n,d36)"),
        (687, "(n,d37)"),
        (688, "(n,d38)"),
        (689, "(n,d39)"),
        (690, "(n,d40)"),
        (691, "(n,d41)"),
        (692, "(n,d42)"),
        (693, "(n,d43)"),
        (694, "(n,d44)"),
        (695, "(n,d45)"),
        (696, "(n,d46)"),
        (697, "(n,d47)"),
        (698, "(n,d48)"),
        (699, "(n,dc)"),
        (700, "(n,t0)"),
        (701, "(n,t1)"),
        (702, "(n,t2)"),
        (703, "(n,t3)"),
        (704, "(n,t4)"),
        (705, "(n,t5)"),
        (706, "(n,t6)"),
        (707, "(n,t7)"),
        (708, "(n,t8)"),
        (709, "(n,t9)"),
        (710, "(n,t10)"),
        (711, "(n,t11)"),
        (712, "(n,t12)"),
        (713, "(n,t13)"),
        (714, "(n,t14)"),
        (715, "(n,t15)"),
        (716, "(n,t16)"),
        (717, "(n,t17)"),
        (718, "(n,t18)"),
        (719, "(n,t19)"),
        (720, "(n,t20)"),
        (721, "(n,t21)"),
        (722, "(n,t22)"),
        (723, "(n,t23)"),
        (724, "(n,t24)"),
        (725, "(n,t25)"),
        (726, "(n,t26)"),
        (727, "(n,t27)"),
        (728, "(n,t28)"),
        (729, "(n,t29)"),
        (730, "(n,t30)"),
        (731, "(n,t31)"),
        (732, "(n,t32)"),
        (733, "(n,t33)"),
        (734, "(n,t34)"),
        (735, "(n,t35)"),
        (736, "(n,t36)"),
        (737, "(n,t37)"),
        (738, "(n,t38)"),
        (739, "(n,t39)"),
        (740, "(n,t40)"),
        (741, "(n,t41)"),
        (742, "(n,t42)"),
        (743, "(n,t43)"),
        (744, "(n,t44)"),
        (745, "(n,t45)"),
        (746, "(n,t46)"),
        (747, "(n,t47)"),
        (748, "(n,t48)"),
        (749, "(n,tc)"),
        (750, "(n,3He0)"),
        (751, "(n,3He1)"),
        (752, "(n,3He2)"),
        (753, "(n,3He3)"),
        (754, "(n,3He4)"),
        (755, "(n,3He5)"),
        (756, "(n,3He6)"),
        (757, "(n,3He7)"),
        (758, "(n,3He8)"),
        (759, "(n,3He9)"),
        (760, "(n,3He10)"),
        (761, "(n,3He11)"),
        (762, "(n,3He12)"),
        (763, "(n,3He13)"),
        (764, "(n,3He14)"),
        (765, "(n,3He15)"),
        (766, "(n,3He16)"),
        (767, "(n,3He17)"),
        (768, "(n,3He18)"),
        (769, "(n,3He19)"),
        (770, "(n,3He20)"),
        (771, "(n,3He21)"),
        (772, "(n,3He22)"),
        (773, "(n,3He23)"),
        (774, "(n,3He24)"),
        (775, "(n,3He25)"),
        (776, "(n,3He26)"),
        (777, "(n,3He27)"),
        (778, "(n,3He28)"),
        (779, "(n,3He29)"),
        (780, "(n,3He30)"),
        (781, "(n,3He31)"),
        (782, "(n,3He32)"),
        (783, "(n,3He33)"),
        (784, "(n,3He34)"),
        (785, "(n,3He35)"),
        (786, "(n,3He36)"),
        (787, "(n,3He37)"),
        (788, "(n,3He38)"),
        (789, "(n,3He39)"),
        (790, "(n,3He40)"),
        (791, "(n,3He41)"),
        (792, "(n,3He42)"),
        (793, "(n,3He43)"),
        (794, "(n,3He44)"),
        (795, "(n,3He45)"),
        (796, "(n,3He46)"),
        (797, "(n,3He47)"),
        (798, "(n,3He48)"),
        (799, "(n,3Hec)"),
        (800, "(n,a0)"),
        (801, "(n,a1)"),
        (802, "(n,a2)"),
        (803, "(n,a3)"),
        (804, "(n,a4)"),
        (805, "(n,a5)"),
        (806, "(n,a6)"),
        (807, "(n,a7)"),
        (808, "(n,a8)"),
        (809, "(n,a9)"),
        (810, "(n,a10)"),
        (811, "(n,a11)"),
        (812, "(n,a12)"),
        (813, "(n,a13)"),
        (814, "(n,a14)"),
        (815, "(n,a15)"),
        (816, "(n,a16)"),
        (817, "(n,a17)"),
        (818, "(n,a18)"),
        (819, "(n,a19)"),
        (820, "(n,a20)"),
        (821, "(n,a21)"),
        (822, "(n,a22)"),
        (823, "(n,a23)"),
        (824, "(n,a24)"),
        (825, "(n,a25)"),
        (826, "(n,a26)"),
        (827, "(n,a27)"),
        (828, "(n,a28)"),
        (829, "(n,a29)"),
        (830, "(n,a30)"),
        (831, "(n,a31)"),
        (832, "(n,a32)"),
        (833, "(n,a33)"),
        (834, "(n,a34)"),
        (835, "(n,a35)"),
        (836, "(n,a36)"),
        (837, "(n,a37)"),
        (838, "(n,a38)"),
        (839, "(n,a39)"),
        (840, "(n,a40)"),
        (841, "(n,a41)"),
        (842, "(n,a42)"),
        (843, "(n,a43)"),
        (844, "(n,a44)"),
        (845, "(n,a45)"),
        (846, "(n,a46)"),
        (847, "(n,a47)"),
        (848, "(n,a48)"),
        (849, "(n,ac)"),
        (875, "(n,2n0)"),
        (876, "(n,2n1)"),
        (877, "(n,2n2)"),
        (878, "(n,2n3)"),
        (879, "(n,2n4)"),
        (880, "(n,2n5)"),
        (881, "(n,2n6)"),
        (882, "(n,2n7)"),
        (883, "(n,2n8)"),
        (884, "(n,2n9)"),
        (885, "(n,2n10)"),
        (886, "(n,2n11)"),
        (887, "(n,2n12)"),
        (888, "(n,2n13)"),
        (889, "(n,2n14)"),
        (890, "(n,2n15)"),
        (891, "(n,2nc)"),
        (901, "heating-local"),
    ]
    .iter()
    .cloned()
    .collect()
});

/// A static HashMap that maps ENDF reaction descriptive names to their MT numbers
/// for reverse lookup functionality. Includes special case for 'fission' -> 18.
pub static REACTION_MT: Lazy<HashMap<&'static str, i32>> = Lazy::new(|| {
    [
        ("(n,total)", 1),
        ("(n,elastic)", 2),
        ("(n,nonelastic)", 3),
        ("(n,level)", 4),
        ("(n,misc)", 5),
        ("(n,2nd)", 11),
        ("(n,2n)", 16),
        ("(n,3n)", 17),
        ("(n,fission)", 18),
        ("(n,f)", 19),
        ("(n,nf)", 20),
        ("(n,2nf)", 21),
        ("(n,na)", 22),
        ("(n,n3a)", 23),
        ("(n,2na)", 24),
        ("(n,3na)", 25),
        ("(n,absorption)", 27),
        ("(n,np)", 28),
        ("(n,n2a)", 29),
        ("(n,2n2a)", 30),
        ("(n,nd)", 32),
        ("(n,nt)", 33),
        ("(n,n3He)", 34),
        ("(n,nd2a)", 35),
        ("(n,nt2a)", 36),
        ("(n,4n)", 37),
        ("(n,3nf)", 38),
        ("(n,2np)", 41),
        ("(n,3np)", 42),
        ("(n,n2p)", 44),
        ("(n,npa)", 45),
        ("(n,n1)", 51),
        ("(n,n2)", 52),
        ("(n,n3)", 53),
        ("(n,n4)", 54),
        ("(n,n5)", 55),
        ("(n,n6)", 56),
        ("(n,n7)", 57),
        ("(n,n8)", 58),
        ("(n,n9)", 59),
        ("(n,n10)", 60),
        ("(n,n11)", 61),
        ("(n,n12)", 62),
        ("(n,n13)", 63),
        ("(n,n14)", 64),
        ("(n,n15)", 65),
        ("(n,n16)", 66),
        ("(n,n17)", 67),
        ("(n,n18)", 68),
        ("(n,n19)", 69),
        ("(n,n20)", 70),
        ("(n,n21)", 71),
        ("(n,n22)", 72),
        ("(n,n23)", 73),
        ("(n,n24)", 74),
        ("(n,n25)", 75),
        ("(n,n26)", 76),
        ("(n,n27)", 77),
        ("(n,n28)", 78),
        ("(n,n29)", 79),
        ("(n,n30)", 80),
        ("(n,n31)", 81),
        ("(n,n32)", 82),
        ("(n,n33)", 83),
        ("(n,n34)", 84),
        ("(n,n35)", 85),
        ("(n,n36)", 86),
        ("(n,n37)", 87),
        ("(n,n38)", 88),
        ("(n,n39)", 89),
        ("(n,n40)", 90),
        ("(n,nc)", 91),
        ("(n,disappear)", 101),
        ("(n,gamma)", 102),
        ("(n,p)", 103),
        ("(n,d)", 104),
        ("(n,t)", 105),
        ("(n,3He)", 106),
        ("(n,a)", 107),
        ("(n,2a)", 108),
        ("(n,3a)", 109),
        ("(n,2p)", 111),
        ("(n,pa)", 112),
        ("(n,t2a)", 113),
        ("(n,d2a)", 114),
        ("(n,pd)", 115),
        ("(n,pt)", 116),
        ("(n,da)", 117),
        ("(n,5n)", 152),
        ("(n,6n)", 153),
        ("(n,2nt)", 154),
        ("(n,ta)", 155),
        ("(n,4np)", 156),
        ("(n,3nd)", 157),
        ("(n,nda)", 158),
        ("(n,2npa)", 159),
        ("(n,7n)", 160),
        ("(n,8n)", 161),
        ("(n,5np)", 162),
        ("(n,6np)", 163),
        ("(n,7np)", 164),
        ("(n,4na)", 165),
        ("(n,5na)", 166),
        ("(n,6na)", 167),
        ("(n,7na)", 168),
        ("(n,4nd)", 169),
        ("(n,5nd)", 170),
        ("(n,6nd)", 171),
        ("(n,3nt)", 172),
        ("(n,4nt)", 173),
        ("(n,5nt)", 174),
        ("(n,6nt)", 175),
        ("(n,2n3He)", 176),
        ("(n,3n3He)", 177),
        ("(n,4n3He)", 178),
        ("(n,3n2p)", 179),
        ("(n,3n2a)", 180),
        ("(n,3npa)", 181),
        ("(n,dt)", 182),
        ("(n,npd)", 183),
        ("(n,npt)", 184),
        ("(n,ndt)", 185),
        ("(n,np3He)", 186),
        ("(n,nd3He)", 187),
        ("(n,nt3He)", 188),
        ("(n,nta)", 189),
        ("(n,2n2p)", 190),
        ("(n,p3He)", 191),
        ("(n,d3He)", 192),
        ("(n,3Hea)", 193),
        ("(n,4n2p)", 194),
        ("(n,4n2a)", 195),
        ("(n,4npa)", 196),
        ("(n,3p)", 197),
        ("(n,n3p)", 198),
        ("(n,3n2pa)", 199),
        ("(n,5n2p)", 200),
        ("(n,Xp)", 203),
        ("(n,Xd)", 204),
        ("(n,Xt)", 205),
        ("(n,X3He)", 206),
        ("(n,Xa)", 207),
        ("heating", 301),
        ("damage-energy", 444),
        ("(n,p0)", 600),
        ("(n,p1)", 601),
        ("(n,p2)", 602),
        ("(n,p3)", 603),
        ("(n,p4)", 604),
        ("(n,p5)", 605),
        ("(n,p6)", 606),
        ("(n,p7)", 607),
        ("(n,p8)", 608),
        ("(n,p9)", 609),
        ("(n,p10)", 610),
        ("(n,p11)", 611),
        ("(n,p12)", 612),
        ("(n,p13)", 613),
        ("(n,p14)", 614),
        ("(n,p15)", 615),
        ("(n,p16)", 616),
        ("(n,p17)", 617),
        ("(n,p18)", 618),
        ("(n,p19)", 619),
        ("(n,p20)", 620),
        ("(n,p21)", 621),
        ("(n,p22)", 622),
        ("(n,p23)", 623),
        ("(n,p24)", 624),
        ("(n,p25)", 625),
        ("(n,p26)", 626),
        ("(n,p27)", 627),
        ("(n,p28)", 628),
        ("(n,p29)", 629),
        ("(n,p30)", 630),
        ("(n,p31)", 631),
        ("(n,p32)", 632),
        ("(n,p33)", 633),
        ("(n,p34)", 634),
        ("(n,p35)", 635),
        ("(n,p36)", 636),
        ("(n,p37)", 637),
        ("(n,p38)", 638),
        ("(n,p39)", 639),
        ("(n,p40)", 640),
        ("(n,p41)", 641),
        ("(n,p42)", 642),
        ("(n,p43)", 643),
        ("(n,p44)", 644),
        ("(n,p45)", 645),
        ("(n,p46)", 646),
        ("(n,p47)", 647),
        ("(n,p48)", 648),
        ("(n,pc)", 649),
        ("(n,d0)", 650),
        ("(n,d1)", 651),
        ("(n,d2)", 652),
        ("(n,d3)", 653),
        ("(n,d4)", 654),
        ("(n,d5)", 655),
        ("(n,d6)", 656),
        ("(n,d7)", 657),
        ("(n,d8)", 658),
        ("(n,d9)", 659),
        ("(n,d10)", 660),
        ("(n,d11)", 661),
        ("(n,d12)", 662),
        ("(n,d13)", 663),
        ("(n,d14)", 664),
        ("(n,d15)", 665),
        ("(n,d16)", 666),
        ("(n,d17)", 667),
        ("(n,d18)", 668),
        ("(n,d19)", 669),
        ("(n,d20)", 670),
        ("(n,d21)", 671),
        ("(n,d22)", 672),
        ("(n,d23)", 673),
        ("(n,d24)", 674),
        ("(n,d25)", 675),
        ("(n,d26)", 676),
        ("(n,d27)", 677),
        ("(n,d28)", 678),
        ("(n,d29)", 679),
        ("(n,d30)", 680),
        ("(n,d31)", 681),
        ("(n,d32)", 682),
        ("(n,d33)", 683),
        ("(n,d34)", 684),
        ("(n,d35)", 685),
        ("(n,d36)", 686),
        ("(n,d37)", 687),
        ("(n,d38)", 688),
        ("(n,d39)", 689),
        ("(n,d40)", 690),
        ("(n,d41)", 691),
        ("(n,d42)", 692),
        ("(n,d43)", 693),
        ("(n,d44)", 694),
        ("(n,d45)", 695),
        ("(n,d46)", 696),
        ("(n,d47)", 697),
        ("(n,d48)", 698),
        ("(n,dc)", 699),
        ("(n,t0)", 700),
        ("(n,t1)", 701),
        ("(n,t2)", 702),
        ("(n,t3)", 703),
        ("(n,t4)", 704),
        ("(n,t5)", 705),
        ("(n,t6)", 706),
        ("(n,t7)", 707),
        ("(n,t8)", 708),
        ("(n,t9)", 709),
        ("(n,t10)", 710),
        ("(n,t11)", 711),
        ("(n,t12)", 712),
        ("(n,t13)", 713),
        ("(n,t14)", 714),
        ("(n,t15)", 715),
        ("(n,t16)", 716),
        ("(n,t17)", 717),
        ("(n,t18)", 718),
        ("(n,t19)", 719),
        ("(n,t20)", 720),
        ("(n,t21)", 721),
        ("(n,t22)", 722),
        ("(n,t23)", 723),
        ("(n,t24)", 724),
        ("(n,t25)", 725),
        ("(n,t26)", 726),
        ("(n,t27)", 727),
        ("(n,t28)", 728),
        ("(n,t29)", 729),
        ("(n,t30)", 730),
        ("(n,t31)", 731),
        ("(n,t32)", 732),
        ("(n,t33)", 733),
        ("(n,t34)", 734),
        ("(n,t35)", 735),
        ("(n,t36)", 736),
        ("(n,t37)", 737),
        ("(n,t38)", 738),
        ("(n,t39)", 739),
        ("(n,t40)", 740),
        ("(n,t41)", 741),
        ("(n,t42)", 742),
        ("(n,t43)", 743),
        ("(n,t44)", 744),
        ("(n,t45)", 745),
        ("(n,t46)", 746),
        ("(n,t47)", 747),
        ("(n,t48)", 748),
        ("(n,tc)", 749),
        ("(n,3He0)", 750),
        ("(n,3He1)", 751),
        ("(n,3He2)", 752),
        ("(n,3He3)", 753),
        ("(n,3He4)", 754),
        ("(n,3He5)", 755),
        ("(n,3He6)", 756),
        ("(n,3He7)", 757),
        ("(n,3He8)", 758),
        ("(n,3He9)", 759),
        ("(n,3He10)", 760),
        ("(n,3He11)", 761),
        ("(n,3He12)", 762),
        ("(n,3He13)", 763),
        ("(n,3He14)", 764),
        ("(n,3He15)", 765),
        ("(n,3He16)", 766),
        ("(n,3He17)", 767),
        ("(n,3He18)", 768),
        ("(n,3He19)", 769),
        ("(n,3He20)", 770),
        ("(n,3He21)", 771),
        ("(n,3He22)", 772),
        ("(n,3He23)", 773),
        ("(n,3He24)", 774),
        ("(n,3He25)", 775),
        ("(n,3He26)", 776),
        ("(n,3He27)", 777),
        ("(n,3He28)", 778),
        ("(n,3He29)", 779),
        ("(n,3He30)", 780),
        ("(n,3He31)", 781),
        ("(n,3He32)", 782),
        ("(n,3He33)", 783),
        ("(n,3He34)", 784),
        ("(n,3He35)", 785),
        ("(n,3He36)", 786),
        ("(n,3He37)", 787),
        ("(n,3He38)", 788),
        ("(n,3He39)", 789),
        ("(n,3He40)", 790),
        ("(n,3He41)", 791),
        ("(n,3He42)", 792),
        ("(n,3He43)", 793),
        ("(n,3He44)", 794),
        ("(n,3He45)", 795),
        ("(n,3He46)", 796),
        ("(n,3He47)", 797),
        ("(n,3He48)", 798),
        ("(n,3Hec)", 799),
        ("(n,a0)", 800),
        ("(n,a1)", 801),
        ("(n,a2)", 802),
        ("(n,a3)", 803),
        ("(n,a4)", 804),
        ("(n,a5)", 805),
        ("(n,a6)", 806),
        ("(n,a7)", 807),
        ("(n,a8)", 808),
        ("(n,a9)", 809),
        ("(n,a10)", 810),
        ("(n,a11)", 811),
        ("(n,a12)", 812),
        ("(n,a13)", 813),
        ("(n,a14)", 814),
        ("(n,a15)", 815),
        ("(n,a16)", 816),
        ("(n,a17)", 817),
        ("(n,a18)", 818),
        ("(n,a19)", 819),
        ("(n,a20)", 820),
        ("(n,a21)", 821),
        ("(n,a22)", 822),
        ("(n,a23)", 823),
        ("(n,a24)", 824),
        ("(n,a25)", 825),
        ("(n,a26)", 826),
        ("(n,a27)", 827),
        ("(n,a28)", 828),
        ("(n,a29)", 829),
        ("(n,a30)", 830),
        ("(n,a31)", 831),
        ("(n,a32)", 832),
        ("(n,a33)", 833),
        ("(n,a34)", 834),
        ("(n,a35)", 835),
        ("(n,a36)", 836),
        ("(n,a37)", 837),
        ("(n,a38)", 838),
        ("(n,a39)", 839),
        ("(n,a40)", 840),
        ("(n,a41)", 841),
        ("(n,a42)", 842),
        ("(n,a43)", 843),
        ("(n,a44)", 844),
        ("(n,a45)", 845),
        ("(n,a46)", 846),
        ("(n,a47)", 847),
        ("(n,a48)", 848),
        ("(n,ac)", 849),
        ("(n,2n0)", 875),
        ("(n,2n1)", 876),
        ("(n,2n2)", 877),
        ("(n,2n3)", 878),
        ("(n,2n4)", 879),
        ("(n,2n5)", 880),
        ("(n,2n6)", 881),
        ("(n,2n7)", 882),
        ("(n,2n8)", 883),
        ("(n,2n9)", 884),
        ("(n,2n10)", 885),
        ("(n,2n11)", 886),
        ("(n,2n12)", 887),
        ("(n,2n13)", 888),
        ("(n,2n14)", 889),
        ("(n,2n15)", 890),
        ("(n,2nc)", 891),
        ("heating-local", 901),
        ("fission", 18),
    ]
    .iter()
    .cloned()
    .collect()
});
