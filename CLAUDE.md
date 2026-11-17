# CLAUDE.md — eCR_evolver

## Coding principles

- **No features beyond what was asked.** No speculative abstractions, configurability, or error handling for impossible scenarios.
- **Surgical changes.** Touch only what the task requires. Don't improve adjacent code, comments, or formatting. Match existing style.
- **Surface tradeoffs before coding.** If multiple interpretations exist, present them — don't pick silently. If something is unclear, ask.
- **If you notice unrelated dead code or issues, mention them — don't silently fix them.**

## STATUS: scaffold (2026-11-17). No model code. Blocked on the calibration panel.

The architecture is not decided, and the experiment that decides it has not run.
See [`docs/prereg_calibration_panel.md`](docs/prereg_calibration_panel.md).
**Do not start a model here before that panel reports.** Writing the wrong one of
the two architectures below means writing it twice.

## What this project is

**ML-directed evolution of effector domains**, starting with the 70-aa SNH domain
of SS18. It is **step 1 of the four-step eCR manuscript arc** (Kuang et al., same
lab, extracted SNH from SS18 and optimized it to 70 aa; this repo is the
computational half of that optimization).

The parent sequences are already in `eCR_mod_lib`:

| Record | Organism | UniProt | Length |
|---|---|---|---|
| `ED_SS18_SNH` | *Mus musculus* | Q62280 | 70 aa |
| `ED_SS18_N70` | *Homo sapiens* | Q15532 | 70 aa |

Mouse parent (the MEF→mES-relevant one):

```
MSVAFAAPRQRGKGEITPAAIQKMLDENNHLIQCIMDYQNKGKASECSQYQQILHTNLVYLATIADSNQN
```

Boundary provenance is asymmetric, and the weaker side is the parent we evolve:
the **human** record carries `source_doi = 10.1073/pnas.1802242115` and
`validation_level = ChIP-validated`; the **mouse** record carries **no DOI**
(`source = "Curated: "`, empty) and its notes read "boundary per project spec",
while still being labelled `screen-validated`. Flagging, not fixing — if the
mouse 70-aa boundary is to anchor a directed-evolution campaign, that field
should cite whatever it was taken from.

## Thesis placement — do not let this become the headline

The manuscript's thesis is **modular reconstitution**; directed evolution is a
*method inside step 1*. An improved SNH variant is an input to steps 2–4, not the
result the paper is about.

## The fitness function (decided 2026-11-17)

**Differential reporter: activity minus activity under BRG1 loss.** The library
is screened twice — once normally, once with SMARCA4 degraded or inhibited — and
selection is on the difference.

Why not the two obvious alternatives:

- **Raw reporter activity is the wrong objective** (decided 2026-11-16). SS18-SNH
  acts by docking BRG1 and recruiting BAF. The cheapest way for a 70-aa peptide
  to raise reporter output is to acquire a generic acidic activation patch — i.e.
  to evolve into VP16 and lose the mechanism the paper is about.
- **A BRG1-recruitment readout cannot be the in-loop fitness.** Reporter screens
  (rTetR-tethered, FACS-sorted — the Bintu/Tycko HT-recruit platform behind
  several of the library's integrated ED sources: HiTEff 2020, Ludwig 2023,
  Mukund 2023, Tycko 2025) run 10⁴–10⁵ variants; BRG1 occupancy by
  CUT&RUN or a proximity assay runs 10⁰–10². Three to four orders of magnitude.
  Mechanistic virtue does not buy throughput.

The differential screen gets mechanism specificity **at reporter throughput**: a
variant that gained a generic AD scores the same in both arms and nets ~zero; a
variant that genuinely improved BAF docking loses its signal when BAF is gone.

The recruitment readout appears exactly twice, at the two ends where 10¹–10² is
enough: on the calibration panel up front, and on the winners at the end.

## The architecture fork — the calibration panel decides it

| Panel outcome | What this repo becomes |
|---|---|
| differential reporter **tracks** BRG1 occupancy | a **scoring** model: deep mutational scanning of the 70 aa, trained on screen output, propose-and-rescreen rounds, library 10⁴–10⁵ |
| it **does not** | a **proposer** model: sequence- or structure-guided design of ~10¹–10² variants per round, evaluated directly by the expensive readout |

Different data pipeline, different model, different repo layout.

## Governing constraints (inherited)

- **Library round-trip.** Evolved variants must land in `eCR_mod_lib` as ED
  records under its existing schema, so they stay comparable to the ~7,700
  entries already there. This is an argument for running the screen on the same
  assay family the library's score columns come from (rTetR/HT-recruit), not a
  novel one.
- **Variant names must be unique against the library.** `eCR_predictor`'s
  `fuse.py` `_load_eds` keys ED rows by name with `setdefault` (defect F8a in
  `eCR_mod_lib/docs/DBD_EXTRACTION_FIXES.md`): a name collision is silently
  dropped and the variant never reaches a design run. Until F8a lands, check
  every emitted name against `library/module_library.tsv` before export.
- **Mechanism specificity is the objective, not activity.** Any fitness that can
  be maximized without BAF is rejected by construction.
- **Negative rounds are first-class**, as refusal is in the navigator: a round
  that fails to beat the parent is reported as a negative result, not resolved
  into a "best variant".
- **Species is a parameter.** The mouse parent drives MEF→mES; the human ortholog
  differs at 4 of the 70 positions (domain numbering: 27 E/D, 38 Y/S, 44 A/T,
  53 I/M), so the two are near-isogenic and the human record is a natural
  near-neutral variant rather than a separate design problem. Never hard-code one.

## Sibling repos (all under github.com/JimYuhaoWu)

`eCR_mod_lib` (main) · `eCR_navigator` (main) · `eCR_predictor` (master) ·
`eCR_ranker` (closed, negative result 2026-09-19) · `eCR_perturber` (feasibility
test failed 2026-09-19, disposition pending) · **eCR_evolver** (this repo)

## Where this sits in the arc

1. **Step 1 — this repo.** Optimize the SNH effector.
2. Step 2 — nominate other EDs that cooperate with Nanog (`eCR_predictor`; gated
   on F8a, since ED nomination over the library would otherwise search a library
   missing 1,439 of its 6,969 ED rows).
3. Step 3 — Nanog + better ED vs native Nanog vs extracted Nanog DBD + ED
   (gated on F7: the extracted chassis is 72 aa with 233 aa / 76% silently
   dropped, including the tryptophan-repeat dimerization module).
4. Step 4 — generalize that DBD + ED works.

Readout across steps 2–4 is `eCR_navigator/docs/prereg_nanog_readout.md`.
