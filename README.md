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
HPV-ABM-Atlanta/
├── data/
│   ├── raw/                  # Public datasets (ACS, NHANES, NIS-Teen)
│   ├── processed/            # Cleaned CSVs used for model initialization
│   └── README.md             # Data dictionary and links to sources
│
├── src/
│   ├── initialize_population.py     # Create synthetic population from ACS data
│   ├── generate_network.py          # Construct sexual contact network
│   ├── simulate_hpvm_abm.py         # Core ABM simulation engine
│   ├── vaccination_scenarios.py     # Functions for coverage/efficacy variations
│   └── analysis_utils.py            # Aggregation, plots, and sensitivity analysis
│
├── notebooks/
│   ├── exploratory_analysis.ipynb   # Data exploration and preprocessing
│   ├── parameter_sensitivity.ipynb  # Example scenario evaluation
│
├── results/
│   ├── figures/                     # Generated plots (prevalence, incidence)
│   └── summary_stats.csv            # Key metrics across scenarios
│
├── reports/
│   ├── project_milestone.pdf
│   └── final_report.pdf
│
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 🧾 Data Sources

| Dataset | Description | Source |
|----------|--------------|--------|
| **ACS (2022)** | Demographic distribution by age, sex, and race for Atlanta MSA. | U.S. Census Bureau |
| **NHANES / NSFG** | Sexual behavior parameters (partner count, concurrency). | CDC / NCHS |
| **CDC NIS-Teen (2023)** | HPV vaccination coverage by age and sex (state-level). | CDC |
| **BRFSS (2022)** | Behavioral health indicators (planned integration). | CDC |
| **Literature-derived parameters** | HPV transmission probabilities, clearance rates, vaccine efficacy. | Walboomers et al. (1999); Clifford et al. (2003); Stuart et al. (2024) |

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

---

## ▶️ Usage

### Run the baseline simulation
```bash
python src/simulate_hpvm_abm.py --coverage 0.7 --efficacy 0.9 --timesteps 120
```

### Example (inside Python)
```python
from src.simulate_hpvm_abm import run_simulation

results = run_simulation(coverage=0.7, efficacy=0.9, steps=120)
results.plot_prevalence_curve()
```

### Generate scenario comparisons
```bash
python src/vaccination_scenarios.py --scenarios scenario_config.json
```

---

## 📊 Output
- **Prevalence curves** over time by vaccination scenario  
- **Incidence rate** of new infections per month  
- **Coverage–efficacy sensitivity plots**  
- **Summary statistics CSV** for downstream analysis  

All output figures and tables are saved under `results/`.

---

## 🧩 Limitations
- Regional sexual behavior data for Atlanta are approximated using national survey priors (NHANES, NSFG).  
- Vaccination behavior and network dynamics are modeled as static during simulation.  
- Calibration and validation are ongoing as additional CDC and BRFSS data are incorporated.

---

## 👥 Authors

| Name | Role |
|------|------|
| **Scott Watanuki** | Data processing, visualization, and analysis |
| **Hailey Toeppner** | Model implementation, simulation, and evaluation |
| **Anish Arora** | Mathematical modeling and literature review |

Course: *CSE 8803 – Data Science for Epidemiology (Georgia Tech)*  
Instructor: *Dr. B. Aditya Prakash*

---

## 📚 References

- Stuart, R. M., Kerr, C. C., et al. *HPVsim: An agent-based model of HPV transmission and cervical disease.* PLOS Comput Biol, 2024.  
- Anderson, E., Goodreau, S. M., et al. *HIV and STI Epidemic Potential of Networks of MSM in Two Cities.* J Infect Dis, 2021.  
- Walboomers, J. M. M., et al. *Human papillomavirus is a necessary cause of invasive cervical cancer worldwide.* J Pathol, 1999.  
- Clifford, G. M., et al. *HPV types in invasive cervical cancer worldwide: a meta-analysis.* Br J Cancer, 2003.  
- Campos, N. G., et al. *An Updated Natural History Model of Cervical Cancer: Derivation of Model Parameters.* Value in Health, 2014.  
- Bénard, É., Laprise, J.-F., Brisson, M. *Potential Population-Level Effectiveness of One-Dose HPV Vaccination.* Lancet Infect Dis, 2023.  

---

## 🧾 License
This repository is licensed under the **MIT License**.  
See [`LICENSE`](LICENSE) for details.

---

## 🧠 Citation
If you use this codebase or model design in academic work, please cite:

```
@misc{HPVABMAtlanta2025,
  author = {Watanuki, Scott and Hailey, Toeppner and Anish, Arora},
  title = {HPV-ABM-Atlanta: Agent-Based Modeling of HPV Transmission and Vaccination Dynamics in the Atlanta Metropolitan Area},
  year = {2025},
  note = {Course Project, Georgia Institute of Technology}
}
```
