# Calibration-panel freeze artefacts

Built by [`scripts/build_panel.py`](../scripts/build_panel.py) to discharge two items
on the freeze checklist of [`docs/prereg_calibration_panel.md`](../docs/prereg_calibration_panel.md):
**§8.7** (conserved vs non-conserved position assignment, with the alignment it came
from) and **§8.4** (the panel sequence list as a TSV with a recorded sha256).

Rebuild from the committed caches — no network, no arguments:

```
python scripts/build_panel.py
```

`--fetch` refreshes the two raw inputs instead of using the caches. Both are committed
precisely so that a later UniProt release cannot silently move the position classes
under a frozen pre-registration.

`.gitattributes` marks `panel/**` as `-text`, so git never rewrites line endings here.
Without that, the sha256 values below would only match on the platform that computed
them, which would make the freeze hashes useless.

## Files

| File | What it is | sha256 |
|---|---|---|
| `ss18_orthologs_raw.tsv` | raw UniProt stream, 1,064 entries | `49b039e55b6b300a380b0fcb923d87738ca5a5175eb4f1e8b9f7bc0ef39ac5c9` |
| `7vrb.cif` | PDB 7VRB coordinates | `44b7cd0254dce3e9e95a6231e192b572aa1633187f92727d4251f63aecb01380` |
| `snh_orthologs.tsv` | the alignment: one ungapped 70-aa SNH per species, 175 rows | `7c844012211ba56fa3bb129aae6f7db2eedb22e40295e59c06a305c02bc8e0cb` |
| `snh_positions.tsv` | per-position conservation × BRG1 contact × secondary structure, 70 rows | `97da087ed55484feba9970c9c546202faf8c2edb1a31cb10f9a6b5f0fc540a4d` |
| `calibration_panel_v1.tsv` | the panel, 52 constructs | `ce71982137d12cb467cd18b080daca9c8205dda56ad7901dad6fca4ea461760b` |

## Input 1 — the orthologue alignment

UniProt release **2026_04** (04-November-2026), query:

```
gene_exact:SS18 AND taxonomy_id:7742 AND keyword:KW-1185
```

(vertebrate SS18, reference proteomes only) → 1,064 entries over 466 species. One
representative per species, reviewed entries preferred, then longest. Each full-length
sequence is locally aligned (BLOSUM62, gap −11/−1) to the mouse 70-aa parent and kept
only if the alignment covers all 70 positions with no gap on either side.

**175 species survive.** Of the 291 dropped, 281 are dropped because the SNH is
*incomplete in the entry*, not because it is divergent — these are N-terminally
truncated gene predictions that begin inside the domain (`A0A218V136`, for instance,
begins at the residue aligning to mouse position 53). Slicing a 70-mer out of those
would have injected non-homologous sequence into every conservation estimate. 10 more
are dropped for gapped alignments.

Both parents round-trip to their exact known sequences (`Q62280` → mouse,
`Q15532` → human); the build aborts if they do not.

Clade spread: Mammalia 85, Actinopterygii 46, Aves 20, Lepidosauria 13, Amphibia 5,
Testudines 4, Crocodylia 1, Coelacanthimorpha 1. Mammalia supplies 85 of 175 (49%),
well above its share of living vertebrate species, so `frac_identical_to_mouse` is
mammal-weighted; the per-clade column in `snh_orthologs.tsv` is there if a
clade-balanced recomputation is ever wanted. Minimum identity to mouse is 0.800
(several teleost fish), so the domain is highly conserved and the variable positions
stand out against a near-invariant background.

## Input 2 — the BRG1 contact surface

**PDB 7VRB**, "Structure of the Human BRG1/SS18 complex", X-ray, 2.39 Å. It is a
*single-chain fusion*: entity positions 3–44 are SMARCA4 (`P51532`) 172–213, the HSA
helix, joined by a TEV linker to entity positions 56–143, which are SS18 (`Q15532`)
14–101. So SS18 residue *p* sits at entity position *p* + 42, and SNH positions 14–70
are resolved while 1–13 are absent from the construct.

For every resolved SNH residue the script takes the minimum heavy-atom distance to any
SMARCA4 residue, computed separately within a chain and between chains across all four
copies. **22 residues contact at ≤ 5 Å intra-chain and only 3 inter-chain** — so the
fusion folds back on itself in the biological topology rather than pairing with a
neighbour in the lattice, and the intra-chain distance is the one used. The 3
inter-chain contacts are at positions 68–70, where the construct continues past the
70-aa SNH boundary; they are crystal packing.

Secondary structure is taken from the UniProt `Q15532` HELIX features, which derive
from this same structure: **α1 = 18–40, α2 = 44–69**, with 1–17 coil and 41–43 the
inter-helix linker.

## The position classes (§8.7)

| Class | Rule | n |
|---|---|---|
| invariant | identical to mouse in **all 175** orthologues | 50 |
| variable | identical in ≤ 98% | 15 |
| unassigned | in between | 5 |

The unassigned band (positions 6, 10, 14, 19, 37 — all at 98.3–99.4%) is deliberately
left out of both panel classes, so neither is contaminated by borderline positions.

Two results fell out that were not designed for, and both are checks on the method:

- **Every one of the 22 BRG1-contacting residues is invariant across all 175 species.**
  Conservation and the contact surface were measured from completely independent
  sources — an alignment and a crystal structure — and they agree exactly.
- **All four human/mouse differences (27, 38, 44, 53) are off the interface** (6.0, 8.6,
  10.0 and 6.8 Å) and are the four least-conserved positions in the whole domain. That
  is why the two parents are near-isogenic in function, as `CLAUDE.md` states.

A third result is worth the designer's attention: **at positions 44 and 53 the mouse
residue is the minority allele**, present in 5.1% and 2.9% of vertebrates, while the
human residue is the vertebrate consensus (T at 94.9%, M at 88.6%). Mouse SS18 is the
outlier at two of the four positions, not human. So `EV_SNH_A44T` and `EV_SNH_I53M`
convert the mouse parent toward the consensus and decompose two quarters of the
human/mouse difference; they are partially redundant with `EV_SNH_WT_HUMAN`, which
carries all four at once.

## The panel (§8.4)

52 constructs, of which **49 enter the primary ρ** (the 3 floor members are excluded by
prereg §5). 49 is one above the n ≥ 48 power floor §6 derives; see the §8.4 note in the
pre-registration for why no position was dropped to hit 48 exactly.

| Class | n | In ρ | Rule |
|---|---|---|---|
| `parent_wt` | 1 | yes | mouse `ED_SS18_SNH` |
| `parent_ortholog` | 1 | yes | human `ED_SS18_N70` |
| `paralog` | 2 | yes | SS18L1 / SS18L2 library records, trimmed to 1–70 |
| `conserved_interface` | 22 | yes | every invariant BRG1-contacting residue → Ala (→Gly where the parent is Ala) |
| `nonconserved_noninterface` | 15 | yes | every variable non-contacting residue → the commonest alternative in the alignment |
| `helix_disrupting` | 5 | yes | mid-helix Pro at invariant **non-contacting** positions 23, 31, 50, 57, 63 |
| `off_mechanism` | 3 | yes | `ED_VP16_TAD`, `ED_VP64`, `ED_TYCKO_2025_ADS_DUX4_1` |
| `floor` | 3 | no | no-effector, plus 2 composition-matched shuffles |

### Why these off-mechanism controls

§4.1 asks for "a minimal VP16, or a top library activator with no BAF link", and the
library has both.

- **`ED_VP16_TAD`** is the prototypical acidic AD, and its library
  `known_interactors` field reads **MED25, TBP, TAF4, TAF9** — Mediator and TFIID, with
  no mSWI/SNF subunit. For this one, off-mechanism is *recorded in the library*, not an
  assumption drawn from memory.
- **`ED_VP64`** is the VP16 tetramer `eCR_predictor`'s own designs use as the opener
  ED. It is precisely the domain `CLAUDE.md` says SNH must not evolve into, which makes
  it the most load-bearing thing on the panel to be able to tell SNH apart from. It
  shares VP16-TAD's acidic core, so these two are **one mechanism, not two**.
- **`ED_TYCKO_2025_ADS_DUX4_1`** adds a second, independent mechanism: a strong
  HT-recruit activator (−9.73 `Max Avg Act at CD2 log2(OFF:ON)`) that is not an
  mSWI/SNF subunit.

Whether any of them is truly BAF-independent is the premise §7.2 tests, not something
this panel assumes: if they fail to separate on the differential axis, the panel is
void by that control. The design checks its own premise.

### Other design notes that matter when reading the results

- The helix-disrupting prolines sit at invariant positions that are **not** at the
  interface, so a dead readout means the fold was lost, not that a contact was deleted.
  That keeps them interpretable as prereg §7.4's dead control rather than as a 23rd
  interface mutant.
- `EV_SNH_A7P` introduces a proline but is a *near-neutral* member, not a disruptor:
  position 7 is in the N-terminal coil and proline is the commonest alternative there
  (23.4% of vertebrates).
- Positions 7, 8 and 12 are **unresolved in 7VRB**, so their non-interface status is an
  inference from their being N-terminal to both helices, not a measurement. The column
  records this as `brg1_interface = unknown`.
- The paralogs are trimmed from their 100-aa and 80-aa library records to 70 aa, which
  resolves the "trim or treat separately" question §4.1 raised in favour of trimming —
  every panel member is then the same length at the same tether. The consequence is
  that their library HT-recruit scores (−7.40, −5.67) were measured on the longer
  constructs and are **not** panel measurements.

## Two library defects this surfaced

Reported, not fixed, per `CLAUDE.md`:

1. **`name` is not a key in `eCR_mod_lib`'s ED table.** 1,622 ED names are used by more
   than one record, and **2,167 of 7,701 ED records are lost if keyed by name alone**.
   `LEUTX_1` alone matches five records, **three of them repressor-subtype**, and a
   name-keyed lookup resolves it to `ED_TYCKO_2025_RDS_LEUTX_1` — a *repressor*. An
   early draft of this builder keyed by name and so selected a repressor domain as an
   "off-mechanism activator" control; it now keys on `module_id` and asserts that
   `module_id` is unique. This is defect F8a seen from the data side. Note that
   `eCR_evolver/CLAUDE.md` currently quotes "1,439 of its 6,969 ED rows", which was
   measured against an older library.
2. **`ED_TYCKO_2025_ADS_SS18L2_1` carries `uniprot_id = P10070`**, which is
   `GLI2_HUMAN`, not SS18L2. The correct accession is `Q9UHA2`, as the sibling record
   `CR_EPIFACTORS_SS18L2` already has.
