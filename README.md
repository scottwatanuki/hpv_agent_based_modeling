# HPV-ABM-Atlanta  
### Agent-Based Modeling of HPV Transmission and Vaccination Dynamics in the Atlanta Metropolitan Area  

---

## 📘 Overview  
This project develops an **agent-based model (ABM)** to simulate the spread of **Human Papillomavirus (HPV)** and the impact of **vaccination coverage and efficacy scenarios** within the **Atlanta Metropolitan Area**.  

The model captures individual-level heterogeneity in sexual behavior, demographic structure, and vaccine uptake to assess population-level prevalence and incidence outcomes over time.  
While focused on Atlanta for specificity, the framework is extensible to other urban regions through parameter reconfiguration.

This repository is part of the **CSE 8803: Data Science for Epidemiology** course at Georgia Tech.

---

## 🎯 Research Objectives
1. **Simulate HPV transmission** using an agent-based framework that captures heterogeneity in demographics and sexual networks.  
2. **Evaluate vaccination coverage and efficacy scenarios** to understand population-level reduction in prevalence and incident infections.  
3. **Quantify sensitivity and robustness** of model outputs to behavioral and epidemiological parameters.  
4. **Demonstrate a tractable, interpretable ABM** that balances epidemiological realism with computational feasibility.

---

## 🧠 Conceptual Background
- HPV is the most common sexually transmitted infection, with high-risk types (16/18) responsible for approximately 70% of cervical cancers globally.  
- Vaccination of adolescents before sexual debut is critical for preventing infection, but impact depends on coverage, timing, and network structure.  
- **Agent-based models (ABMs)** enable explicit representation of individuals, partnerships, and vaccination states—offering flexibility beyond traditional compartmental (ODE-based) models.

---

## 🗂️ Repository Structure
```
hpv_agent_based_modeling/
├── data/
│   ├── raw/                  # Public datasets (ACS, NSFG, NIS-Teen)
│   ├── processed/            # Cleaned CSVs used for model initialization
│   └── README.md             # Data dictionary and links to sources
│
├── src/
│   ├── initialize_population.py     # Create synthetic population from ACS data
│   ├── simulate_hpvm_abm.py         # Core ABM simulation engine
│   ├── vaccination_scenarios.py     # Functions for coverage/efficacy variations
│   └── analysis_utils.py            # Aggregation, plots, and sensitivity analysis
│
├── results/
│   ├── figures/                     # Generated plots (prevalence, incidence)
│   └── summary_stats.csv            # Key metrics across scenarios
│
├── doc/
│   ├── project_milestone.pdf        # Milestone report
│   └── final_report.pdf             # Final report
│
├── docs/
│   ├── assets/                      # Assets for the website
│   ├── index.html                   # HTML for the website
│   └── styles.css                   # CSS for the website
|
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 🧾 Data Sources

| Dataset | Description | Source |
|----------|--------------|--------|
| **ACS (2022)** | Demographic distribution by age, sex, and race for Atlanta MSA. | U.S. Census Bureau |
| **NSFG** | Sexual behavior parameters (partner count, concurrency). | CDC / NCHS |
| **CDC NIS-Teen (2023)** | HPV vaccination coverage by age and sex (state-level). | CDC |
| **Literature-derived parameters** | HPV transmission probabilities, clearance rates, vaccine efficacy. |  Stuart et al. (2024); Walker et al. 2012 |

All datasets used are **publicly available** and preprocessed under `data/processed/` for simulation initialization.

---

## ⚙️ Installation & Setup

### Requirements
- Python ≥ 3.9  
- Recommended libraries:  
  ```
  numpy
  pandas
  networkx
  matplotlib
  mesa
  seaborn
  tqdm
  ```

### Clone and setup
```bash
git clone https://github.com/scottwatanuki/hpv_agent_based_modeling.git
cd hpv_agent_based_modeling
pip install -r requirements.txt
```

### Quick start (baseline)
1) Verify data exists at `data/processed/ga_age_sex_cleaned.csv`
  - If missing, place the provided CSV there. Expected columns: `Age Group, Male_Estimate, Female_Estimate`.
2) Run baseline simulation and produce outputs under `results/`:
```bash
python src/simulate_hpvm_abm.py --coverage 0.7 --efficacy 0.9 --timesteps 120 --seed 42
```
3) Outputs:
  - Figure: `results/figures/prevalence_baseline.png`
  - Summary CSV: `results/summary_stats_baseline.csv`
  - JSON: `results/run_meta.json`

---

## ▶️ Usage

### Run the baseline simulation
```bash
python src/simulate_hpvm_abm.py --coverage 0.7 --efficacy 0.9 --timesteps 120 --seed 42
```

### Example (inside Python)
```python
from simulate_hpvm_abm import run_simulation
from analysis_utils import write_prevalence_outputs

results = run_simulation(coverage=0.7, efficacy=0.9, steps=120, seed=42)
results.plot_prevalence_curve()
write_prevalence_outputs(results.prevalence, prefix="baseline")
```

### Generate scenario comparisons
```bash
python src/vaccination_scenarios.py -scenarios scenarios.json
```

Runs a simulation from simulate_hpvm_abm for each of the three different vaccination scenarios:
- Baseline Model
- Race-Stratified Vaccination Model
- Race-Stratified Vaccination Model with Screening
---

## 📊 Output
- **Prevalence curves** over time by vaccination scenario  
- **Incidence rate** of total cancer incidence
- **Coverage–efficacy sensitivity plots**  
- **Summary statistics CSV** for downstream analysis  
- **Summary statistics JSON** for run parameters

All output figures and tables are saved under `results/`.

For the vaccination scenarios, you will see:
- Baseline Figure: `results/figures/prevalence_baseline_uniform.png`
- Race-Stratified Vaccination Figure: `results/figures/prevalence_race_stratified_vaccination.png`
- Race-Stratified Vaccination Model with Screening Figure: `results/figures/prevalence_race_stratified_with_screening.png`
- Baseline JSON: `results/run_meta_baseline_uniform.json`
- Race-Stratified Vaccination JSON: `results/run_meta_race_stratified_vaccination.json`
- Race-Stratified Vaccination Model with Screening JSON: `results/run_meta_race_stratified_with_screening.json`
- Baseline CSV: `results/summary_stats_baseline_uniform.csv`
- Race-Stratified Vaccination CSV: `results/summary_stats_race_stratified_vaccination.csv`
- Race-Stratified Vaccination Model with Screening CSV: `results/summary_stats_race_stratified_with_screening.csv`

---

## 🧩 Limitations
- Regional sexual behavior data for Atlanta are approximated using national survey priors (NSFG).  
- Vaccination behavior and network dynamics are modeled as static during simulation.  
- Large numbers of agents and complicated scenarios run slowly 
