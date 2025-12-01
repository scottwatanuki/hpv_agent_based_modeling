# Preliminary findings — Milestone 1

This section summarizes the preliminary outputs and interpretation from the baseline implementation of the HPV agent-based model. The goal of this milestone was to create a reproducible baseline simulation, validate end-to-end workflow (population initialization → simulation → analysis), and produce initial summary outputs that can be referenced in the milestone report.

## What we implemented (short)
- Modular code: `src/initialize_population.py`, `src/simulate_hpvm_abm.py`, `src/analysis_utils.py` and a scenario scaffold `src/vaccination_scenarios.py`.
- Data wiring: the ACS-like age–sex distribution is read from `data/processed/ga_age_sex_cleaned.csv` and used to synthesize a population.
- Baseline ABM: random-mixing transmission, per-contact transmission probability (β), per-step clearance, and vaccine susceptibility reduction (efficacy). A CLI and programmatic API (`run_simulation`) were provided.
- Reproducibility: runs record `results/run_meta.json` (parameters, seed, timestamp, git SHA).

## Methods (baseline configuration referenced)
- Population size: 1000 (default for quick runs)
- Timesteps: 120 (baseline run used for the milestone)
- Transmission: β = 0.02 per contact
- Clearance: per-step clearance probability = 0.05
- Partners per infectious agent per timestep: 5 (random sampling)
- Vaccination: coverage = 0.7, efficacy = 0.9 (reduces susceptibility)
- Seed: configurable via `--seed` for deterministic runs

These parameters are intentionally simple placeholders to produce reproducible, interpretable behavior for the milestone.

## Key quantitative outputs (baseline)
The baseline run artifacts are saved under `results/` and include a prevalence figure and a one-row summary CSV. Extracted values from the last baseline run:

- Timesteps: 120
- Final prevalence: 0.0
- Mean prevalence: 0.0000333333
- Peak prevalence: 0.001 (noted at timestep 0 in the summary CSV)

Figure: `results/figures/prevalence_baseline.png` (see `reports/figure_captions.md` for the caption text).

## Interpretation and immediate takeaways
- The model produced extremely low prevalence under the baseline parameterization. This outcome is consistent with:
  - a small number of seeded infections at initialization;
  - the relatively low per-contact transmission probability (β = 0.02);
  - explicit vaccine protection among 70% of agents (efficacy = 0.9) that reduces susceptibility; and
  - a non-zero clearance probability that removes infections over time.
- The reported "peak prevalence at timestep 0" is an artefact of how the seed was introduced (initial seed(s) counted at t=0). For epidemiological interpretation we will re-run with multiple seeds and average across replicates (next steps).

## Limitations relevant to these findings
- The current baseline uses random mixing (no persistent partnership network), so results do not capture network structure effects (assortativity, concurrency).
- Vaccination assignment in the baseline is Bernoulli by coverage, not age-targeted; real-world coverage varies by age and sex.
- Single-type HPV simplification: no genotype-level dynamics or multi-strain interactions.
- Parameters are placeholder values and not calibrated to surveillance targets; numbers are illustrative rather than predictive.

## Planned follow-ups that address the above and how they will change interpretation
1. Replace random mixing with a network-based contact model (degree distribution and assortativity). Expect network structure to change outbreak size and heterogeneity in prevalence across subgroups.
2. Add age-targeted vaccination scenarios using state-level NIS-Teen priors; this will enable more realistic coverage–efficacy comparisons.
3. Run replicate simulations and compute summaries with confidence intervals rather than single-run statistics.
4. Calibrate β and clearance to published prevalence/incidence targets (or to a fitted time series) to increase external validity.

## Short conclusion for the milestone report
We implemented a transparent, reproducible ABM baseline and produced initial prevalence outputs that validate the end-to-end pipeline (data → model → analysis). Results show that, under conservative transmission parameters and substantial vaccine protection, prevalence remains very low; however, these outputs are preliminary and subject to the many simplifying assumptions documented in `reports/assumptions.md`. The next milestone will expand mixing realism, add age-targeted vaccination inputs, and run systematic scenario sweeps with replicates to produce robust comparative results.

For any edits or specific phrasing required by the course template, I can adapt this section and insert it directly into the milestone report draft located in `reports/`.