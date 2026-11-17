# Pre-registration — reporter ↔ BRG1 calibration panel

**Repo:** eCR_evolver · **Written:** 2026-11-17 · **Status: DRAFT, not frozen.**
See §8 for the items that must be filled before any data is collected.

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

Target **n = 48** (see §6 for why not 20). Composition is given by the role each
class plays in the analysis; the sequence list itself is a §8 freeze item.

| Class | Purpose | n (target) |
|---|---|---|
| Mouse WT SS18-SNH (`ED_SS18_SNH`) | parent; anchors both axes | ≥3 replicates |
| Human SS18-N70 (`ED_SS18_N70`) | natural near-neutral variant (4 of 70 positions differ) | 1 |
| SS18L1 / SS18L2 N-terminal domains | paralogs with HT-recruit activity already in the library (`ED_TYCKO_2025_ADS_SS18L1_1`, `..._SS18L2_1`; metric `Max Avg Act at CD2 log2(OFF:ON)` = −7.40 and −5.67). **Note these records are 100 aa and 80 aa, not 70-aa-trimmed** — either trim them to the SNH boundary or treat them as a separate class | 2 |
| Conserved-position substitutions | expected to span intermediate to dead | ~15 |
| Non-conserved-position substitutions | expected near-neutral; tests that the assay is not simply reading "any mutation is bad" | ~10 |
| Helix-disrupting substitutions | expected dead on both axes | ~5 |
| **Off-mechanism positive control** | a strong acidic activation domain (e.g. a minimal VP16, or a top library activator with no BAF link). **Must score high on raw reporter, low on differential, ~background on occupancy.** | ≥3 |
| No-effector / scrambled | floor on both axes | ≥3 |

The off-mechanism control is the most important class: it is what distinguishes
"the differential axis works" from "the differential axis is noise centred on
zero".

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

## 7. Controls that must pass, or the panel is void

Checked before the primary statistic is read:

1. **BRG1 loss actually happened** — SMARCA4 depletion (or inhibition efficacy)
   verified independently in the loss arm.
2. **The off-mechanism control separates** — high raw reporter, low differential.
   If it does not, the differential axis failed to discriminate and the primary ρ
   is uninterpretable whatever its value.
3. **WT replicates agree** across arms, within a tolerance stated in advance.
4. **Dead controls are dead on both axes.** A helix-disrupted variant with high
   occupancy means the occupancy readout is measuring something other than
   specific recruitment.

## 8. Freeze checklist — all open as of 2026-11-17

- [ ] **BRG1-loss reagent confirmed.** dTAG/degron believed available; *not yet
      verified*. Falling back to an inhibitor changes §4.2's interpretation and
      must be recorded here before data.
- [ ] **Tethering system fixed** — rTetR at a synthetic reporter (comparable to
      the library's existing ED scores) vs dCas9 at an endogenous locus (more
      relevant to the eventual eCR, not comparable). Recorded before data.
- [ ] **BRG1 occupancy assay fixed** (CUT&RUN vs proximity/BRET), and its
      background definition stated.
- [ ] **Panel sequence list frozen** — the full n = 48 list as a TSV in this
      repo, with its sha256 recorded here, following the practice in
      `eCR_navigator/docs/prereg_nanog_readout.md` §7.
- [ ] **Thresholds in §6 signed off** (0.60 / 0.40) or replaced.
- [ ] **Replicate counts, and the WT tolerance in §7.3, stated numerically.**
- [ ] **Conserved vs non-conserved position assignment** fixed, with the
      alignment it came from committed here.

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

*(none yet)*
