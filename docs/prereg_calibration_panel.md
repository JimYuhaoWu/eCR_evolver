# Pre-registration — reporter ↔ BRG1 calibration panel

**Repo:** eCR_evolver · **Written:** 2026-11-17 · **Updated:** 2026-11-26 ·
**Status: DRAFT, not frozen.**
See §8 for the items that must be filled before any data is collected. The two
computational items are now closed — the panel is built and its positions are
assigned from a committed alignment and crystal structure (§8.4, §8.7). What
remains is four wet-lab or sign-off decisions, none of which this repo can make.

---

## 1. What this decides

Whether a **differential reporter** readout can serve as the in-loop fitness
function for directed evolution of the SS18 SNH domain, and therefore which of
two architectures `eCR_evolver` is:

| If the panel shows | `eCR_evolver` becomes | Library size per round |
|---|---|---|
| differential reporter tracks BRG1 occupancy | a **scoring** model over deep mutational scanning data | 10⁴–10⁵ |
| it does not | a **proposer** model evaluated directly by the recruitment readout | 10¹–10² |

This is a **go/no-go on an architecture**, not a biological claim for the
manuscript. It is pre-registered because the two outcomes have very different
costs, and the temptation to read an ambiguous correlation as permission to build
the cheap version is obvious.

## 2. Why both readouts are needed at all

The two candidate readouts differ on two axes, and the second is binding:

- **What they measure.** A transcriptional reporter integrates everything
  downstream of recruitment. A BRG1 readout measures the one step the SNH domain
  performs — docking the SMARCA4 ATPase.
- **Throughput.** rTetR-tethered, FACS-sorted reporter screens run 10⁴–10⁵
  variants. BRG1 occupancy (CUT&RUN, or a proximity/BRET assay) runs 10⁰–10².

So recruitment cannot be the in-loop fitness however mechanistically correct it
is, and raw reporter activity is the wrong objective, because a 70-aa peptide's
cheapest route to higher output is acquiring a generic acidic activation patch.

The proposal under test resolves this by making the reporter **BAF-conditional**:

```
fitness(variant) = reporter(variant) - reporter(variant | BRG1 lost)
```

A variant that gained a generic AD scores alike in both arms and nets ~zero. A
variant that improved BAF docking loses its signal when BAF is gone and nets
high. The recruitment readout is then needed only twice — here, on this panel,
and at the end on the winners.

## 3. Hypothesis

**H1.** Across a variant panel spanning the activity range, the **differential**
reporter score correlates with BRG1 occupancy at the tethered locus, and does so
**better than the raw** reporter score does.

The second clause is the real content. A differential score that correlates no
better than the raw score has not bought the mechanism specificity it costs an
extra screening arm to obtain.

**H0 / the competing outcome.** Differential reporter and BRG1 occupancy are
uncorrelated or weakly correlated — i.e. the reporter axis, raw or differential,
does not report on recruitment for this domain family.

## 4. Design

### 4.1 Variant panel

Target **n = 48** entering the primary statistic (see §6 for why not 20).
**Built: `panel/calibration_panel_v1.tsv`, 52 constructs, 49 of them in the
primary ρ.** Provenance, the selection rules and the design caveats are in
[`panel/README.md`](../panel/README.md); `scripts/build_panel.py` rebuilds every
artefact from committed inputs.

Position classes come from two independent sources, measured not assumed:
a 175-species vertebrate SS18 alignment (UniProt 2026_04) for conservation, and
PDB **7VRB** — the human BRG1/SS18 complex at 2.39 Å — for the BRG1 contact
surface. **All 22 BRG1-contacting residues are invariant across all 175 species**,
and all four human/mouse differences lie off the interface.

| Class | Purpose | n (built) |
|---|---|---|
| Mouse WT SS18-SNH (`ED_SS18_SNH`) | parent; anchors both axes | 1 (replicates per §7.3) |
| Human SS18-N70 (`ED_SS18_N70`) | natural near-neutral variant (4 of 70 positions differ, all off-interface) | 1 |
| SS18L1 / SS18L2 N-terminal domains | paralogs already in the library (`ED_TYCKO_2025_ADS_SS18L1_1`, `..._SS18L2_1`). **Resolved: trimmed to 1–70** so every panel member is the same length at the same tether. Their library scores (−7.40, −5.67) were measured on the 100-aa and 80-aa constructs and are therefore *not* panel measurements | 2 |
| Conserved-position substitutions | **invariant *and* BRG1-contacting** → Ala (Gly where the parent is Ala). Expected to span intermediate to dead | 22 |
| Non-conserved-position substitutions | **variable *and* non-contacting** → the commonest alternative residue in the alignment, i.e. a natural variant. Expected near-neutral; tests that the assay is not simply reading "any mutation is bad" | 15 |
| Helix-disrupting substitutions | mid-helix Pro at invariant **non-contacting** positions 23, 31, 50, 57, 63, so a dead readout means the fold was lost rather than a contact deleted. Expected dead on both axes | 5 |
| **Off-mechanism positive control** | `ED_VP16_TAD` (the minimal VP16 AD this section asked for; its library `known_interactors` are MED25, TBP, TAF4, TAF9 — Mediator and TFIID, no mSWI/SNF subunit, so off-mechanism is *recorded* rather than assumed), `ED_VP64` (the VP16 tetramer eCR_predictor's designs use as the opener ED — the domain CLAUDE.md says SNH must not evolve into, and so the most load-bearing thing to tell SNH apart from; same acidic core as VP16-TAD, hence one mechanism not two) and `ED_TYCKO_2025_ADS_DUX4_1` (an independent strong screen activator, −9.73, not an mSWI/SNF subunit). **Must score high on raw reporter, low on differential, ~background on occupancy.** | 3 |
| No-effector / scrambled | floor on both axes; excluded from the primary ρ per §5 | 3 |

The off-mechanism control is the most important class: it is what distinguishes
"the differential axis works" from "the differential axis is noise centred on
zero". Its "no BAF link" is a literature assumption rather than a measurement,
and §7.2 is its empirical test — if these three do not separate, the panel is
void by that control, so the design checks its own premise.

One caveat to carry into the analysis: because the off-mechanism controls are
*designed* to be high-raw / low-differential, they contribute a spread point that
widens the ρ difference §5 tests. That is the intended discrimination test, not
an independent replication of it.

### 4.2 Readouts, all three measured on every panel member

1. **Reporter, BAF-intact arm.**
2. **Reporter, BRG1-loss arm** — degron (dTAG-SMARCA4) preferred over inhibitor
   (e.g. BRM014), because a degron removes the protein rather than one of its
   activities and gives a cleaner interpretation of the difference. **Reagent
   availability unconfirmed as of 2026-11-17** (§8).
3. **BRG1 occupancy at the tethered locus** — the low-throughput readout.

Arms 1 and 2 must use the same tethering system and the same locus as each other
and as the eventual screen. Otherwise the correlation measured here does not
transfer to the screen it is meant to license.

## 5. Pre-specified analysis

- **Primary statistic.** Spearman ρ between differential reporter score and BRG1
  occupancy, across all panel members except the no-effector floor class.
  Reported with a bootstrap 95% CI.
- **Secondary statistic.** Spearman ρ between **raw** reporter score and BRG1
  occupancy, same members. Report the **difference** between the two ρ values
  with a paired bootstrap CI. H1's second clause is this difference > 0.
- **Control checks** (§7) are evaluated **before** the primary statistic is read.
- No other orderings, subsets, or transformations are introduced after the data
  exist. Dropping a panel member requires an amendment stating the reason, filed
  before that member's occupancy value is unblinded.

## 6. Decision rule, and the power that makes it decidable

**Proposed thresholds — require sign-off before data (§8):**

| Primary ρ (95% CI) | Decision |
|---|---|
| ρ ≥ 0.60 and CI lower bound > 0.40 | **scoring architecture.** Differential reporter becomes the fitness function. |
| ρ < 0.40, or CI includes 0 | **proposer architecture.** Reporter axis abandoned as fitness. |
| anything else | **undecided** — expand the panel. Do not default to the cheaper build. |

**Why n = 48 and not 20.** With a 95% Fisher-z interval, an observed r = 0.6
gives:

| n | 95% CI on r = 0.6 | width |
|---|---|---|
| 20 | [0.21, 0.82] | 0.61 |
| 30 | [0.31, 0.79] | 0.48 |
| 40 | [0.35, 0.77] | 0.41 |
| 50 | [0.39, 0.75] | 0.37 |

At n = 20 the interval around an "it works" point estimate reaches down to 0.21,
inside the fail region — the experiment cannot distinguish its own two outcomes.
The lower bound only clears 0.40 at n ≈ 48. **A 20-variant panel cannot support
this decision**, which is the main thing this pre-registration adds over
"measure both on a few variants and look".

Fisher-z is used for the sizing calculation on a Pearson r, while the primary
statistic is Spearman with a bootstrap CI, which is slightly more conservative at
this n. The sizing is a floor, not a promise.

**The built panel puts 49 members in the primary ρ, not 48.** The class rules in
§4.1 are derived from the alignment and the structure — 22 invariant contacting
residues, 15 variable non-contacting ones — and they happen to sum to 49. No
position was dropped to reach exactly 48, because the only way to choose which
one to drop would be arbitrary, and 48 is a floor rather than a target. 49 clears
it.

## 7. Controls that must pass, or the panel is void

Checked before the primary statistic is read:

1. **BRG1 loss actually happened** — SMARCA4 depletion (or inhibition efficacy)
   verified independently in the loss arm.
2. **The off-mechanism control separates** — high raw reporter, low differential.
   If it does not, the differential axis failed to discriminate and the primary ρ
   is uninterpretable whatever its value.
3. **WT replicates agree** across arms, within a tolerance stated in advance.
   **Proposed, requires sign-off (§8):** in each arm, the standard deviation of
   the WT replicate scores must be **≤ 1/3 of that arm's WT-minus-floor dynamic
   range**. The tolerance is written scale-free on purpose, because the reporter
   units are not fixed until §8's tethering decision is made; it ties the noise
   floor to the signal the panel actually has to resolve, since replicate noise
   above a third of the usable range leaves a 49-point Spearman unable to read
   the gradient it is there to read.
4. **Dead controls are dead on both axes.** A helix-disrupted variant with high
   occupancy means the occupancy readout is measuring something other than
   specific recruitment.

## 8. Freeze checklist — 2 of 7 closed as of 2026-11-26

The two items this repo could settle on its own are closed. The five that remain
are wet-lab facts or sign-offs; **nothing here is frozen until all of them are.**

- [ ] **BRG1-loss reagent confirmed.** dTAG/degron believed available; *not yet
      verified*. Falling back to an inhibitor changes §4.2's interpretation and
      must be recorded here before data.
- [ ] **Tethering system fixed** — rTetR at a synthetic reporter (comparable to
      the library's existing ED scores) vs dCas9 at an endogenous locus (more
      relevant to the eventual eCR, not comparable). Recorded before data.
      *Blocks §7.3's tolerance, which is written scale-free until this is set.*
- [ ] **BRG1 occupancy assay fixed** (CUT&RUN vs proximity/BRET), and its
      background definition stated.
- [x] **Panel sequence list frozen** (2026-11-26) —
      `panel/calibration_panel_v1.tsv`, 52 constructs, 49 in the primary ρ
      sha256 `ce71982137d12cb467cd18b080daca9c8205dda56ad7901dad6fca4ea461760b`
      Rebuildable offline by `scripts/build_panel.py`; provenance, selection
      rules and caveats in [`panel/README.md`](../panel/README.md). Follows the
      sha256 practice of `eCR_navigator/docs/prereg_nanog_readout.md` §7.
      Three design questions §4.1 had left open are resolved there and are part
      of what this line asks to be signed off: the paralogs are **trimmed to
      1–70**, the off-mechanism controls are **VP16-TAD / VP64 / DUX4** from the
      library, and the panel carries **49** rather
      than exactly 48 members in the ρ (§6).
- [ ] **Thresholds in §6 signed off** (0.60 / 0.40) or replaced.
- [ ] **Replicate counts, and the WT tolerance in §7.3, stated numerically.**
      **Proposed, awaiting sign-off:** 3 biological replicates per construct in
      each reporter arm; BRG1 occupancy at 2 biological replicates per construct,
      raised to 3 for the parent and the floor class, since those two anchor the
      axis and §7.3's tolerance is computed from them. WT tolerance as stated in
      §7.3. These are the minimum that makes the panel decidable, not a budget
      estimate — if the occupancy assay cannot afford 2 × 52, that is an argument
      about panel size and belongs in §6, not a reason to lower the replicate
      count.
- [x] **Conserved vs non-conserved position assignment fixed** (2026-11-26) —
      `panel/snh_positions.tsv`
      sha256 `97da087ed55484feba9970c9c546202faf8c2edb1a31cb10f9a6b5f0fc540a4d`
      with the alignment it came from committed beside it as
      `panel/snh_orthologs.tsv` (175 species)
      sha256 `7c844012211ba56fa3bb129aae6f7db2eedb22e40295e59c06a305c02bc8e0cb`
      and the raw UniProt stream it was derived from
      (`panel/ss18_orthologs_raw.tsv`, release 2026_04)
      sha256 `49b039e55b6b300a380b0fcb923d87738ca5a5175eb4f1e8b9f7bc0ef39ac5c9`.
      A position is **invariant** only if identical in all 175 orthologues and
      **variable** only at ≤ 98% identity; the 5 positions in between are
      assigned to neither class. The BRG1 contact surface comes from PDB 7VRB
      (`panel/7vrb.cif`
      sha256 `44b7cd0254dce3e9e95a6231e192b572aa1633187f92727d4251f63aecb01380`).

## 9. What this does not test

- It does not test that an improved SNH variant reprograms anything. That is
  steps 2–4 of the arc.
- It does not establish that BRG1 occupancy *causes* the reporter output; both
  may track a third property of the variant. The panel licenses a *fitness
  proxy*, not a mechanism claim.
- It says nothing about variants acting through BAF subunits other than SMARCA4 —
  a variant that switched to SMARCA2 would read as off-mechanism here. Accepted
  scope limit; worth a line in the methods if a winner behaves oddly.
- It is not a power calculation for the evolution campaign itself.

## 10. Amendments

Numbered, dated, and appended below, each stating what changed and whether any
data had been seen at the time. Same convention as the navigator
pre-registration.

**Amendment 1 — 2026-11-26. No data of any kind had been collected; the panel
has not been built as DNA, and neither reporter arm nor any occupancy
measurement exists.**

Closed §8 items 4 and 7 by building the panel from two measured sources rather
than by hand: a 175-species vertebrate SS18 alignment and PDB 7VRB. Changes to
the design, all in §4.1 and §6:

- The conserved class became **invariant *and* BRG1-contacting** (22 positions,
  was "~15 conserved") and the non-conserved class became **variable *and*
  non-contacting** (15, was "~10"). The structure supplied the selection rule
  that the earlier draft lacked — without it, choosing 15 of 50 invariant
  positions would have been arbitrary.
- The primary ρ therefore carries **49** members, not 48. §6 treats 48 as a
  floor, and no position was dropped to hit the round number.
- Helix-disrupting positions were constrained to invariant **non-contacting**
  sites, which the earlier draft did not specify, so that class now reads as loss
  of fold rather than as further interface mutants.
- The SS18L1/SS18L2 ambiguity §4.1 flagged ("either trim them or treat them as a
  separate class") is resolved as **trim to 1–70**.
- §7.3's WT tolerance and the replicate counts are now written down numerically,
  marked as proposals pending sign-off.

Nothing in §3 (hypothesis), §5 (analysis) or §6's decision thresholds changed.
