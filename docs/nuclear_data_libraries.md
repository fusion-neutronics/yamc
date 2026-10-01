# Nuclear data libraries

Six libraries are available by keyword. Point a setting at one and the sections it
needs are downloaded on first use and cached under ``~/.cache/yamc``.

<!-- doctest: skip -->
```python
yamc.cross_section_data = "endf-b8.1"
```

| keyword | what it is |
|---------|-------------------------------------------------------------------|
| ``endf-b8.1`` | ENDF/B-VIII.1. The default for photon data and for every unset transmutation part. |
| ``jeff-4.0`` | JEFF-4.0. |
| ``jendl-5.0`` | JENDL-5. |
| ``tendl-2025`` | TENDL-2025. The widest isotope coverage of the six. |
| ``tendl-2017`` | TENDL-2017, for comparing a result against the older release. |
| ``fendl-3.2d`` | FENDL-3.2d. Cross sections only, no transmutation data. |

Every setting also takes a path instead of a keyword, and ``cross_section_data``
additionally takes a directory or a dict keyed by nuclide; the four
``transmutation_*`` parts take one value each, since a network is assembled per
subsection rather than per nuclide. Which setting does what, and what each does
when left alone, is in [Nuclear data settings](nuclear_data.md#nuclear-data-settings).

## What each library publishes

| library | cross sections | ``decay`` | ``reactions`` | ``fission_yields`` | ``branching`` |
| --- | --- | :---: | :---: | :---: | :---: |
| ``endf-b8.1`` | neutron + photon | yes | yes | yes | yes |
| ``jeff-4.0`` | neutron | yes | yes | yes | yes |
| ``jendl-5.0`` | neutron + photon | yes | yes | yes | yes |
| ``fendl-3.2d`` | neutron + photon | -- | -- | -- | -- |
| ``tendl-2025`` | neutron | -- | yes | -- | yes |
| ``tendl-2017`` | neutron | -- | yes | -- | yes |

The ``cross sections`` column is what ``yamc.cross_section_data`` can be set to.
That setting is not transport-only: ``Material.transmute()`` and
``Model.simulate_transmutation()`` fold the same cross sections against the
spectrum to get their reaction rates, so an activation run needs it as much as a
transport run does. The four remaining columns are the transmutation subsections,
each set independently through its own ``yamc.transmutation_*`` setting. A library
with no subsections still works for transport, it just cannot supply a chain.

TENDL is neutron-only, so it has no decay or fission-yield evaluations. Take
those from a library that has them and source ``reactions`` and ``branching`` from
TENDL, which is what its far wider parent coverage is worth having for. Pointing a
transmutation part at a library that does not provide it fails immediately with a
message listing the subsections that library does provide.

``endf-b8.1``, ``jeff-4.0`` and ``jendl-5.0`` each supply a complete,
single-library transmutation chain.

### Photon data is chosen separately

Photons default to ``endf-b8.1`` whatever the neutron library is, which keeps
fluorescence and Auger electrons working when the neutron data is TENDL or FENDL.
Only three libraries can serve as a photon source at all:

- ``endf-b8.1`` and ``jendl-5.0`` publish the photoatomic and atomic-relaxation
  pair a coupled run needs. JENDL-5's is EPICS2017 adopted by JAEA, so it is the
  same data ENDF/B-VIII.1 carries.
- ``fendl-3.2d`` publishes photon cross sections with **no atomic relaxation**, so
  fluorescence is silently suppressed if you point photons at it.
- ``jeff-4.0`` and both TENDL releases publish no photoatomic data at all. JEFF's
  photon sublibrary is photonuclear, which is a different quantity.

There is no global photon setting to move photons with. The only way off the
default is per material, through ``photon_data`` on
[read_nuclear_data](materials.md#loading-nuclear-data), whose values may be
library keywords as well as paths.

## How big each library is

Counted from the built trees. ``--`` means the library does not publish that at
all, rather than publishing an empty one.

| library | neutron | photon | reaction edges | branching channels | decay nuclides | covariance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ``endf-b8.1`` | 557 | 100 | 5,175 | 915 | 3,820 | 42% |
| ``jeff-4.0`` | 593 | -- | 13,513 | 8,303 | 3,851 | 98% |
| ``jendl-5.0`` | 795 | 100 | 9,902 | 5,349 | 4,070 | 13% |
| ``fendl-3.2d`` | 192 | 61 | -- | -- | -- | 35% |
| ``tendl-2025`` | 2,848 | -- | 67,797 | 8,416 | -- | pending |
| ``tendl-2017`` | 2,801 | -- | 68,827 | 8,365 | -- | pending |

``neutron`` and ``photon`` are evaluations, one per nuclide and one per element.
``reaction edges`` is how many ``(parent, reaction, product)`` triples the network
holds, which is what decides whether a product is reachable at all. ``covariance``
is the share of neutron evaluations carrying cross-section covariance.

The spread between the two ends is the thing worth noticing: TENDL carries five
times the evaluations of ENDF/B-VIII.1 and thirteen times the reaction edges, and
has none of the decay data. That is why a TENDL ``reactions`` network with
ENDF/B-VIII.1 ``decay`` and ``fission_yields`` beside it is the usual arrangement.

### Covariance is presence, not coverage

The percentage counts evaluations that carry cross-section covariance, which is
not the same as covariance landing where a given run needs it. Depth varies as
much as presence: JENDL-5 carries covariance for 13% of its evaluations, and its
Fe56 covers 56 channels against ENDF/B-VIII.1's 7.

Cross sections (MF=33) are also the only covariance family that reaches the
converted data. The ENDF parser reads resonance-parameter covariance (MF=32),
angular-distribution covariance (MF=34) and radionuclide-production covariance
(MF=40, the covariance of each product state's production cross section) as
well, but none of them is written to the converted data yet, so nothing samples
them. ENDF/B-VIII.1, FENDL-3.2d and TENDL state no correlation in MF=40 between
a product's states, or with MF=33, so even once it is stored the
isomeric-branching split's uncertainty is only partly stated. An evaluation that
keeps its resolved-resonance uncertainty in MF=32 and leaves MF=33 at zero over
that range arrives with no capture, elastic or fission uncertainty below the top
of its resolved resonances, and that range reaches down to thermal. The
uncertainty report counts rate from energies where the stated variance is zero
as uncovered, so the gap shows in ``rate_fraction_covered``. MF=33 blocks stated
with LB=0 to 2 are also used today over only the lower half of their energy
grid, and short-range LB=8 blocks are folded as relative where the evaluation
states them in barns squared, as [MF=33 blocks read in
part](materials.md#what-it-does-not-cover) describes.

Beyond covariance, the decay data carries each nuclide's half-life with the
standard deviation its evaluation states, and its decay energy with a standard
deviation built from the evaluation's sigmas on the beta, gamma and alpha
components. The evaluation states no correlation between them, so they are
treated as independent and the total's sigma is their quadrature sum. The
``half_life`` and ``decay_energy`` uncertainty sources draw from these. The
converter also stores the sigma each evaluation states on a decay mode's
branching ratio (``branching_ratio_uncertainty`` in ``decay_modes.arrow``). The
published decay data predates it, so a downloaded library carries it only from
the next data version, and no source samples it yet. Fission yields are in the
same state: the converter stores each evaluation's yields with their stated
sigmas (``yield_uncertainties`` in ``fission_yields/evaluated_yields.arrow``),
the published data predates the file, and no source samples them. The sigmas
stated on decay photon lines are not in the converted data yet. [What
it does not cover](materials.md#what-it-does-not-cover) lists what a
nuclear-data uncertainty run holds at nominal.

### Energy range

Every library reaches past 20 MeV for at least some evaluations, and none is
uniform. Counting evaluations by the ceiling they were converted to:

| library | ceiling |
| --- | --- |
| ``endf-b8.1`` | 452 at 20 MeV, 34 at 150, 34 at 200, 28 at 30, rest at 28 to 65 |
| ``jeff-4.0`` | 491 at 200 MeV, 65 at 20, 16 at 150, 10 at 30, rest at 60 to 90 |
| ``jendl-5.0`` | 579 at 200 MeV, 216 at 20 |
| ``fendl-3.2d`` | 109 at 200 MeV, 72 at 150, 11 at 60 to 65, none at 20 |
| ``tendl-2025`` | 200 MeV throughout |
| ``tendl-2017`` | 200 MeV throughout |

ENDF/B-VIII.1 is the only one that is mostly a 20 MeV library, and even its Fe56
and W186 reach 150 MeV. So for a high-energy source it is worth checking the
ceiling of the nuclides that carry the reaction rather than the library's
headline.

The two halves of the code treat that differently, which is worth knowing before
reading a result. In transport a particle above a nuclide's ceiling finds no
cross section up there, and that reads as a rate too low rather than as an
error. In a transmutation collapse, from yamc-core 0.11.0, it is refused: a group
spectrum whose flux reaches past an evaluation's last energy point is rejected
before the fold, naming the nuclide and the energy its evaluation stops at,
because there is nothing up there to fold against and any answer would be
invented. Before that release it folded against the last tabulated value,
carrying it forward, which is a different and wrong answer whenever the top
groups are populated. Groups above
the ceiling carrying no flux cost nothing, which is what makes a padded
spectrum safe.

### What a library can quietly lack

Two things a library can be missing that neither the tables above nor its cross
sections reveal. Both matter for decay heat and shutdown dose rather than for
transport.

**Mean decay energies that are placeholders.** A nuclide whose decay scheme
nobody has evaluated still gets a decay file, with a third of each beta or
electron-capture branch's Q booked to the light particles and a third to the
photons. For an electron-capture emitter that is meaningless, since the
neutrino carries almost all of Q: ENDF/B-VIII.1 books Sn111 at 1.63 MeV per
decay against the 0.69 MeV JENDL-5 evaluates from its decay scheme. The counts
are 1144 records for `endf-b8.1`, 1529 for `jeff-4.0` and 375 for
`jendl-5.0`. They are exotic nuclides, but early decay heat is carried by
exactly the short-lived products that tend to be unevaluated. Every one is
named in the decay subsection's `provenance.json`, in data converted by
yamc-core 0.10.0 or newer; libraries published before that carry the numbers
without the record that says which of them are placeholders.

**Reactions that do not list the isomer they make.** A reaction that can leave
its product in a metastable state has to say so, one entry per final state. One
that lists the ground state alone puts the isomer outside the network
altogether, so no code can make it. `Os190 (n,n')` to Os190m at 1706 keV and
`Ir191 (n,2n)` to Ir190m2 at 377 keV are both carried by `jeff-4.0` and
`tendl-2025` and both absent from `tendl-2017`, which is why an osmium foil's
heat five minutes after a 14 MeV irradiation comes out at a third of the
measured value on that library. Nothing about such a channel looks wrong from
outside: it is present, its cross section is reasonable, and only the state
list is short.

### Temperatures and size

All six are converted at the same six temperatures, 250, 293.6, 600, 900, 1200 and
2500 K, so a library swap never changes which temperature is available.

Size is per nuclide rather than per library, because only the nuclides in the
model are fetched. Fe56 is 65 MB on ``endf-b8.1`` and 25 MB on ``tendl-2025``;
U235 is 199 MB and 63 MB. Transport pulls the ``distributions`` section on top of
what an activation run needs, which is 15 MB of Fe56's 65 MB on ENDF/B-VIII.1.
