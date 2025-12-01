# Assumptions and Limitations (Milestone)

This milestone implements a minimal, reproducible baseline to generate preliminary outputs. Key simplifying assumptions:

## Model structure
- Static random mixing: infectious agents sample k partners uniformly each timestep; no persistent partnerships or explicit network structure yet.
- Single HPV-type abstraction: no genotype-specific dynamics or cross-protection.
- Per-contact transmission probability (β) constant across agents and time.
- Memoryless clearance: infected agents clear with constant per-step probability.
- Discrete, synchronous timesteps; events occur within-step without ordering effects beyond our loop snapshot.

## Vaccination
- Coverage assigned via Bernoulli(coverage) independent of age/sex and time (no cohort targeting, schedules, or waning).
- Efficacy reduces susceptibility only; no effect on infectiousness or duration.
- No catch-up vaccination, dropout, or multi-dose adherence modeling.

## Population & data
- Synthetic demographics derived from ACS-like age–sex distribution CSV; ages sampled uniformly within group ranges.
- No births, deaths, or migration; fixed population size.
- Behavioral heterogeneity beyond age/sex not yet modeled (e.g., partner rate distributions, concurrency).

## Parameters & calibration
- Parameter values (β, clearance_prob, partners_per_step) are literature-inspired placeholders; no calibration performed in this slice.
- No formal sensitivity analysis in this milestone; results are illustrative.

## Outputs & reproducibility
- Deterministic RNG seed option for reproducibility; run metadata recorded in `results/run_meta*.json`.
- Baseline outputs: prevalence curve and summary statistics only.

## Implications for interpretation
- Prevalence levels and dynamics reflect stylized assumptions; not validated against Atlanta-specific epidemiology.
- Vaccine impact here captures only susceptibility reduction under random mixing; network effects may alter outcomes.

## Near-term next steps to relax assumptions
- Introduce sexual contact network (degree distribution, assortativity, concurrency).
- Age/sex-specific vaccination coverage from NIS-Teen; explore waning and timing.
- Calibrate β and clearance to external targets; add sensitivity analysis.
- Extend outputs (incidence, R_eff proxies, scenario comparisons).