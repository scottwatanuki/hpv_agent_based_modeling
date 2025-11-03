# Milestone Draft Outline — HPV-ABM-Atlanta

- Course: CSE 8803 — Data Science for Epidemiology
- Repo: https://github.com/scottwatanuki/hpv_agent_based_modeling
- Date: 2025-11-03

## 1. Background & Significance
- HPV is the most prevalent STI; high-risk types (16/18) drive most cervical cancers.
- Early vaccination with sufficient coverage and timely uptake reduces population-level burden.
- ABMs enable explicit representation of individuals and partnerships, capturing heterogeneity beyond compartmental models.

## 2. Objectives (current milestone scope)
1) Implement a tractable baseline HPV ABM with vaccination effect and clearance.
2) Produce reproducible baseline outputs (prevalence curve, summary stats).
3) Align repository structure with plan; document assumptions and usage.

Longer-term objectives (from proposal):
- Compare vaccination coverage/efficacy scenarios.
- Explore sensitivity to behavioral and epidemiological parameters.
- Extend to network-based mixing and calibration against external data.

## 3. Data
- Public sources planned (ACS, NHANES/NSFG, NIS-Teen). For this milestone:
  - Implemented: ACS-like age–sex distribution CSV at `data/processed/ga_age_sex_cleaned.csv`
    - Columns: `Age Group, Male_Estimate, Female_Estimate`
- Use: Construct synthetic population demographics; vaccination is assigned by target coverage (baseline).

## 4. Methods (Milestone Baseline)
- Population synthesis: sample individuals by age-group/sex weights; assign integer ages from group ranges.
- Vaccination: Bernoulli(coverage) independent of age/sex for baseline; efficacy reduces susceptibility.
- Transmission: random mixing per timestep; infected agents sample `k` partners; per-contact transmission prob = β × susceptibility.
- Clearance: per-timestep probability of infection clearance (memoryless).
- Parameters (baseline):
  - population_size=1000, steps=120, partners_per_step=5, β=0.02, clearance_prob=0.05, coverage=0.7, efficacy=0.9, seed=42.
- Reproducibility:
  - Deterministic seed option; outputs saved under `results/` with figure and CSV.

## 5. Preliminary Results
- Baseline run artifacts (generated):
  - Figure: `results/figures/prevalence_baseline.png`
  - Summary: `results/summary_stats_baseline.csv`
- Interpretation (brief): prevalence rises from a single seed, constrained by vaccination and clearance; shape reflects random mixing and simple parameters.
- Note: Not calibrated; numbers are illustrative for the milestone.

## 6. Assumptions & Limitations
- Static random mixing; no explicit partnership network yet.
- Single HPV-type abstraction; no genotype structure.
- Vaccination independent of age/sex and timing; no waning.
- Parameters literature-inspired placeholders; no calibration/validation in this slice.

## 7. Next Steps (toward full project)
- Scenario scaffold: JSON-driven coverage/efficacy sweeps; batch outputs and comparison plots.
- Network mixing: degree distribution and assortativity via NetworkX; concurrency options.
- Data integration: age/sex-specific vaccination priors (NIS-Teen), behavior (NHANES/NSFG).
- Calibration and sensitivity: grid/Latin hypercube; tornado plots.
- Reproducibility metadata: `results/run_meta.json` with params, seed, timestamp, git SHA.

## 8. Deliverables (this milestone)
- Code: `src/simulate_hpvm_abm.py`, `src/initialize_population.py`, `src/analysis_utils.py`.
- Outputs: baseline figure and summary CSV under `results/`.
- Docs: updated `README.md`, this outline.
- Environment: `requirements.txt`.

## 9. Timeline (high level)
- Milestone (≈30%): baseline ABM + outputs + docs (completed).
- Midpoint (≈60%): scenario runner, basic network mixing, initial sensitivity.
- Final (100%): calibrated scenarios, full analysis figures, report integration.