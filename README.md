# eCR_evolver

**ML-directed evolution of effector domains.** Starting case: the 70-aa SNH domain of
SS18, the opener half of the eCR designs this platform builds.

Part of the eCR platform: [`eCR_mod_lib`](https://github.com/JimYuhaoWu/eCR_mod_lib) ·
[`eCR_navigator`](https://github.com/JimYuhaoWu/eCR_navigator) ·
[`eCR_predictor`](https://github.com/JimYuhaoWu/eCR_predictor) ·
[`eCR_ranker`](https://github.com/JimYuhaoWu/eCR_ranker) *(closed, negative)* ·
[`eCR_perturber`](https://github.com/JimYuhaoWu/eCR_perturber) *(feasibility failed)* ·
**eCR_evolver**

> **Status: scaffold (2026-11-17). No model code.** The architecture is undecided and the
> experiment that decides it has not run — see
> [`docs/prereg_calibration_panel.md`](docs/prereg_calibration_panel.md).

## Where it sits

Step 1 of the four-step eCR arc: optimize the effector, then nominate cooperating EDs,
then test whether DBD + ED reconstitutes the native factor, then generalize. The thesis
is **modular reconstitution** — directed evolution is a method inside step 1, not the
headline.

## The problem this repo has to solve first

SS18-SNH works by docking BRG1 and recruiting BAF. So the fitness function has to reward
*recruitment*, but the only readout with the throughput for a variant library is a
transcriptional reporter, which rewards output by any mechanism — and the cheapest way
for a 70-aa peptide to raise reporter output is to acquire a generic acidic activation
patch, i.e. to stop being an SNH domain.

The proposal is a **differential reporter**: screen the library with BAF intact and with
BRG1 removed, and select on the difference. That buys mechanism specificity at reporter
throughput. Whether it actually tracks BRG1 recruitment is an empirical question, and it
decides whether this repo is a scoring model over 10⁴–10⁵ variants or a proposer of
10¹–10². The pre-registration in `docs/` is that experiment.

## Layout

```
CLAUDE.md                          project context, constraints, the architecture fork
docs/prereg_calibration_panel.md   the go/no-go experiment on the fitness function
```
