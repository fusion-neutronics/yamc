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

Every setting also takes a path, a directory or a per-nuclide dict instead of a
keyword. See [Options for setting nuclear data source](nuclear_data.md#options-for-setting-nuclear-data-source).

## What each library publishes

| library | cross sections | ``decay`` | ``reactions`` | ``fission_yields`` | ``branching`` |
| --- | --- | :---: | :---: | :---: | :---: |
| ``endf-b8.1`` | neutron + photon | yes | yes | yes | yes |
| ``jeff-4.0`` | neutron | yes | yes | yes | yes |
| ``jendl-5.0`` | neutron + photon | yes | yes | yes | yes |
| ``fendl-3.2d`` | neutron + photon | -- | -- | -- | -- |
| ``tendl-2025`` | neutron | -- | yes | -- | yes |
| ``tendl-2017`` | neutron | -- | yes | -- | yes |

The ``cross sections`` column is what ``yamc.cross_section_data`` can be set to for
transport; the four remaining columns are the transmutation subsections, each set
independently through its own ``yamc.transmutation_*`` setting. A library with no
subsections still works for transport, it just cannot supply a chain.

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

Cross sections are also the only covariance family that reaches the converted
data. The ENDF parser reads angular-distribution and radionuclide-production
covariance too (MF=34 and MF=40, against MF=33 for cross sections), and nothing
consumes them.

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
and W186 reach 150 MeV. A source with particles above a nuclide's ceiling finds
no cross section up there, which reads as a rate that is too low rather than as
an error, so for a high-energy source it is worth checking the ceiling of the
nuclides that carry the reaction rather than the library's headline.

### Temperatures and size

All six are converted at the same six temperatures, 250, 293.6, 600, 900, 1200 and
2500 K, so a library swap never changes which temperature is available.

Size is per nuclide rather than per library, because only the nuclides in the
model are fetched. Fe56 is 65 MB on ``endf-b8.1`` and 25 MB on ``tendl-2025``;
U235 is 199 MB and 63 MB. Transport pulls the ``distributions`` section on top of
what an activation run needs, which is 15 MB of Fe56's 65 MB on ENDF/B-VIII.1.
