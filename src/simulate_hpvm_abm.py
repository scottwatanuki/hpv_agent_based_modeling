import argparse
from dataclasses import dataclass
from typing import Optional, Dict
from pathlib import Path

from src.analysis_utils import write_run_metadata
from src.initialize_population import synthesize_population, Population

import numpy as np
import matplotlib.pyplot as plt

# --- 1. POPULATION STRUCTURE DEFINITIONS ---

# Age bands
AGE_BANDS = [
    (0, 14),   # 0
    (15, 24),  # 1 (Sexually active start)
    (25, 34),  # 2
    (35, 44),  # 3
    (45, 120), # 4
]
NUM_AGE_BANDS = len(AGE_BANDS)

# Sexual Activity Classes (Based on partner acquisition rates, simplified from NATSAL/POLYMOD) - Walker et Al. source will talk ab in paper
# 0: Low, 1: Moderate, 2: Mod-High, 3: High
SEX_ACTIVITY_CLASSES = ["Low", "Mod", "Mod-High", "High"]
NUM_SEX_CLASSES = len(SEX_ACTIVITY_CLASSES)

# Distribution of population across the 4 sexual activity classes (e.g., 60% Low, 30% Mod, 7% Mod-High, 3% High)
# This is a critical input parameter based on survey data (like NATSAL/NSSHB)
# NOTE: This distribution is simplified and should be calibrated to real data.
SEX_ACTIVITY_DISTRIBUTION = np.array([0.60, 0.30, 0.07, 0.03])


def age_to_band_index(age: int) -> int:
    """Maps agent's age to the correct age band index."""
    for i, (lo, hi) in enumerate(AGE_BANDS):
        if lo <= age <= hi:
            return i
    return NUM_AGE_BANDS - 1

# --- 2. AGE-SPECIFIC & NATURAL HISTORY PARAMETERS ---

# 1. TRANSMISSION RATE STRUCTURE (ALPHA_AGE): Same as original
RELATIVE_CONTACT_BY_AGE_BAND = np.array([0.00, 1.00, 0.70, 0.40, 0.20])

# 2. CLEARANCE RATE (GAMMA): Same as original
CLEARANCE_BY_AGE_BAND = np.array([0.00, 0.060, 0.043, 0.035, 0.025])

# 3. CANCER PROGRESSION PARAMETER (Tuning knob for calibration)
# We start with the value derived from the 20-year latency period (0.004167) as a baseline,
# but this must be recalibrated after mixing is introduced.
PROB_PROGRESS_TO_CANCER = 0.0002 # Set lower to account for the log factor and prevent explosion

# --- 3. MIXING MATRIX PARAMETERS (Walker et al., 2012) ---

# Fraction choosing assortatively by age (epsilon_A in the paper)
EPSILON_A = 0.7
# Fraction choosing assortatively by sexual activity class (epsilon_S in the paper)
EPSILON_S = 0.3

# Age-Assortative Matrix (A_mix): Simplified from original
# This matrix will ONLY be used for the non-proportional component of the mixing (Semi-Assortative).
# We use a semi-assortative model where partners are preferentially chosen from adjacent age bands.
# A simplified semi-assortative probability distribution centered around age difference 0
# 5 rows (i, index of choosing agent) x 5 columns (j, index of partner age)
AGE_MIXING_PREFERENCE = np.array([
    [1.0, 0.0, 0.0, 0.0, 0.0], # 0-14 (only mix with self, though not active)
    [0.3, 1.0, 0.7, 0.0, 0.0], # 15-24
    [0.0, 0.7, 1.0, 0.6, 0.0], # 25-34
    [0.0, 0.0, 0.6, 1.0, 0.4], # 35-44
    [0.0, 0.0, 0.0, 0.4, 1.0], # 45+
])

# Normalize the AGE_MIXING_PREFERENCE matrix by row (so each row sums to 1) for easier use later
for i in range(NUM_AGE_BANDS):
    row_sum = AGE_MIXING_PREFERENCE[i].sum()
    if row_sum > 0:
        AGE_MIXING_PREFERENCE[i] /= row_sum


def calculate_mixing_matrix(N: int,
                            age_band_idx: np.ndarray,
                            sex_activity_idx: np.ndarray,
                            partnerships_offered_per_agent: np.ndarray) -> np.ndarray:
    """
    Calculates the revised mixing matrix (rho_kihjm, Equation 16 from Walker et al. 2012).

    The resulting matrix M[i, h, j, m] represents P(partner in age j, class m | agent in age i, class h).

    N_i,h,j,m = NUM_AGE_BANDS, NUM_SEX_CLASSES, NUM_AGE_BANDS, NUM_SEX_CLASSES (5x4x5x4)
    """

    # 1. Calculate P_k'jm (Total Partnerships Offered by subgroup j, m)
    # P_k'jm = N_k'jm * c_k'jm (Population in group * Average partners in group)
    # First, aggregate the total partnership offerings (P_k'jm) for each of the 20 subgroups (age x activity)
    P_subgroup = np.zeros((NUM_AGE_BANDS, NUM_SEX_CLASSES))

    # Iterate through all agents to sum up total partnerships offered in each subgroup
    for i in range(N):
        age_i = age_band_idx[i]
        class_h = sex_activity_idx[i]
        P_subgroup[age_i, class_h] += partnerships_offered_per_agent[i]

    # 2. Calculate Denominators (Sums from Walker's Equation 16)

    # Sum over all age groups (Sum_alpha P_k'alpha, m) for Term 3 (A_P and S_A)
    # Result is a vector 1xNUM_SEX_CLASSES
    sum_over_age = P_subgroup.sum(axis=0)

    # Sum over all sex activity classes (Sum_beta P_k'j, beta) for Term 2 (A_A and S_P)
    # Result is a vector NUM_AGE_BANDS x 1
    sum_over_class = P_subgroup.sum(axis=1)

    # Sum over all subgroups (Sum_alpha Sum_beta P_k'alpha, beta) for Term 4 (A_P and S_P)
    # Result is a scalar
    sum_total = P_subgroup.sum()

    # Initialize the 4D Mixing Matrix: rho_kihjm (Agent i, class h -> Partner j, class m)
    # In this implementation, the matrix will be flattened to 2D for easier indexing in the main loop:
    # M[Source Subgroup Index] -> M[Partner Subgroup Index]
    # Total subgroups = 5 * 4 = 20
    M = np.zeros((NUM_AGE_BANDS * NUM_SEX_CLASSES, NUM_AGE_BANDS * NUM_SEX_CLASSES))

    # Iterate through all 20 * 20 possible subgroup transitions (Source: i, h -> Target: j, m)
    for i in range(NUM_AGE_BANDS): # Source Age Index
        for h in range(NUM_SEX_CLASSES): # Source Class Index
            source_idx = i * NUM_SEX_CLASSES + h # Flattened source index

            for j in range(NUM_AGE_BANDS): # Target Age Index
                for m in range(NUM_SEX_CLASSES): # Target Class Index
                    target_idx = j * NUM_SEX_CLASSES + m # Flattened target index

                    rho = 0.0
                    P_jm = P_subgroup[j, m] # P_k'jm

                    # Term 1: Assortative by Age (A_A) and Assortative by Class (S_A)
                    # P(A_A & S_A) * delta_ij * delta_hm
                    term1 = EPSILON_A * EPSILON_S * (1 if i == j and h == m else 0.0)
                    rho += term1

                    # Term 2: Assortative by Age (A_A) and Proportional by Class (S_P)
                    # P(A_A & S_P) * [P_k'jm / Sum_beta P_k'j, beta] * delta_ij
                    term2 = 0.0
                    if i == j and sum_over_class[j] > 0:
                        term2 = EPSILON_A * (1.0 - EPSILON_S) * (P_jm / sum_over_class[j])
                    rho += term2

                    # Term 3: Proportional by Age (A_P) and Assortative by Class (S_A)
                    # P(A_P & S_A) * [P_k'jm / Sum_alpha P_k'alpha, m] * delta_hm
                    term3 = 0.0
                    if h == m and sum_over_age[m] > 0:
                        term3 = (1.0 - EPSILON_A) * EPSILON_S * (P_jm / sum_over_age[m])
                    rho += term3

                    # Term 4: Proportional by Age (A_P) and Proportional by Class (S_P)
                    # P(A_P & S_P) * [P_k'jm / Sum_total P_k'alpha, beta]
                    term4 = 0.0
                    if sum_total > 0:
                        term4 = (1.0 - EPSILON_A) * (1.0 - EPSILON_S) * (P_jm / sum_total)
                    rho += term4

                    M[source_idx, target_idx] = rho

    # Normalize matrix rows: Each row M[source_idx, :] must sum to 1.
    row_sums = M.sum(axis=1)
    # Only normalize rows that have positive sum (non-sexually active groups will be 0)
    M[row_sums > 0] /= row_sums[row_sums > 0, np.newaxis]

    return M


@dataclass
class SimulationResults:
    """Container for simulation outputs."""

    prevalence: list[float]
    cancer_incidence: list[float]

    def plot_prevalence_curve(self, show: bool = True, save_path: Optional[str] = None):
        plt.figure(figsize=(10, 6))
        # Updated plot function to display both metrics
        plt.plot(self.prevalence, color="tab:blue", linewidth=2, label="% Currently Infected (Prevalence)")
        plt.plot(self.cancer_incidence, color="tab:red", linestyle='--', linewidth=2, label="% Total Cancer Cases (Cumulative Incidence)")

        plt.xlabel("Time step (Months)")
        plt.ylabel("Percentage of Population")
        plt.title("HPV-ABM: Prevalence and Cancer Incidence (Revised Mixing)")
        plt.grid(True, alpha=0.3)
        plt.legend()

        if save_path:
            plt.savefig(save_path, dpi=200, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close()


def run_simulation(
    coverage: float = 0.5,
    efficacy: float = 0.9,
    steps: int = 300,
    *,
    population_size: int = 1000,
    partners_per_step: int = 5,
    BETA_SCALING_FACTOR_C: float = 0.001, # Increased 10x to sustain endemic transmission
    seed: Optional[int] = None,
    # Optional per-race overrides: mapping race -> coverage (0-1) and race -> annual screening uptake (0-1)
    coverage_by_race: Optional[Dict[str, float]] = None,
    screening_by_race: Optional[Dict[str, float]] = None,
    # Probability that a screening event detects & clears infection when it occurs
    screening_detection_effectiveness: float = 0.9,
) -> SimulationResults:
    """
    Run an HPV transmission simulation with the revised, heterogeneous, age-activity mixing matrix.
    """
    if seed is not None:
        np.random.seed(seed)

    # --- POPULATION INITIALIZATION (Modified to include sexual activity class) ---
    # synthesize_population is assumed to be able to handle these groups
    pop: Population = synthesize_population(
        n=population_size,
        coverage=coverage,
        seed=seed,
        coverage_by_race=coverage_by_race,
        screening_by_race=screening_by_race,
    )
    pop_df = pop.df

    ages = pop_df["age"].to_numpy(dtype=int)
    vaccinated = pop_df["vaccinated"].to_numpy(dtype=bool)

    N = len(pop_df)
    age_band_idx = np.array([age_to_band_index(a) for a in ages], dtype=int)

    # screening probabilities (monthly) per agent
    if "screening_prob_monthly" in pop_df.columns:
        screening_prob_monthly = pop_df["screening_prob_monthly"].to_numpy(dtype=float)
    else:
        screening_prob_monthly = np.zeros(N, dtype=float)

    # NEW: Assign sexual activity class based on the defined distribution
    sex_activity_idx = np.random.choice(
        NUM_SEX_CLASSES,
        size=N,
        p=SEX_ACTIVITY_DISTRIBUTION
    )

    # The average number of partners offered per step is now dependent on the activity class.
    # **** Needs to be calibrated ****
    ACTIVITY_PARTNER_MULTIPLIER = np.array([1, 4, 10, 30])

    # Partnerships offered by the agent in this timestep, dependent on their activity class.
    partnerships_offered_per_agent = ACTIVITY_PARTNER_MULTIPLIER[sex_activity_idx]

    # --- CALCULATE REVISED MIXING MATRIX ---
    # This calculation is run once before the simulation loop.
    REVISED_MIXING_MATRIX = calculate_mixing_matrix(
        N,
        age_band_idx,
        sex_activity_idx,
        partnerships_offered_per_agent
    )

    # --- PER-AGENT PARAMETERS ---
    # The base transmission rate is now scaled by the ACTIVITY_PARTNER_MULTIPLIER (heterogeneity)
    # The overall beta (BETA_PER_AGENT) is the calibration factor * age effect * activity effect
    BETA_PER_AGENT = BETA_SCALING_FACTOR_C * RELATIVE_CONTACT_BY_AGE_BAND[age_band_idx]
    clearance_per_agent = CLEARANCE_BY_AGE_BAND[age_band_idx]

    # --- INITIALIZE STATE VARIABLES ---
    infected = np.zeros(N, dtype=bool)
    time_infected = np.zeros(N, dtype=int)
    cancerous = np.zeros(N, dtype=bool)

    num_initial = 10
    initial_indices = np.random.choice(N, size=num_initial, replace=False)
    infected[initial_indices] = True
    time_infected[initial_indices] = 1

    prevalence = []
    cancer_incidence = []

    for t in range(steps):

        # --- PARTNER SELECTION & TRANSMISSION ---
        currently_infected = np.where(infected)[0]

        for i in currently_infected:
            source_age_idx = age_band_idx[i]
            source_class_idx = sex_activity_idx[i]
            source_flat_idx = source_age_idx * NUM_SEX_CLASSES + source_class_idx

            # Get the probability vector M[source_flat_idx, :] for this infected agent
            # This vector gives the probability of partnering with ANY other subgroup (age x activity)
            partner_probabilities_by_subgroup = REVISED_MIXING_MATRIX[source_flat_idx, :]

            if partner_probabilities_by_subgroup.sum() == 0:
                continue

            # The simplified approach: the agent finds (partners_per_step) partners,
            # with probability proportional to the total partnership count (P_subgroup)
            # weighted by the mixing matrix.

            # Identify all potential partners' flat indices
            all_partner_flat_indices = np.arange(NUM_AGE_BANDS * NUM_SEX_CLASSES)

            # --- SIMPLIFIED PARTNER SELECTION WITH MIXING MATRIX ---

            # For simplicity and to integrate the matrix, we'll revert to finding agents
            # whose characteristics match the mix.

            # Get the expected number of partners this agent will form this step
            expected_partners = partners_per_step * (
                ACTIVITY_PARTNER_MULTIPLIER[source_class_idx] / ACTIVITY_PARTNER_MULTIPLIER.mean()
            )
            num_partners_to_find = int(np.round(expected_partners))

            if num_partners_to_find == 0:
                 continue

            # We create a probability vector over the entire population (N) based on the mixing matrix.
            # 1. Map each agent in the population (j) to their target subgroup (j_age, j_class)
            # 2. Find the probability M[source_flat_idx, j_flat_idx]

            partner_agent_probabilities = np.zeros(N)
            for j in range(N):
                j_age = age_band_idx[j]
                j_class = sex_activity_idx[j]
                j_flat_idx = j_age * NUM_SEX_CLASSES + j_class

                # Probability of choosing agent j is proportional to M[source, target]
                partner_agent_probabilities[j] = REVISED_MIXING_MATRIX[source_flat_idx, j_flat_idx]

            partner_agent_probabilities[i] = 0.0 # Exclude self

            if partner_agent_probabilities.sum() == 0:
                continue

            # Normalize and sample from the entire population
            partner_agent_probabilities /= partner_agent_probabilities.sum()


            available_partners = np.where(partner_agent_probabilities > 0)[0]
            num_partners_to_find = min(num_partners_to_find, available_partners.size)

            partners = np.random.choice(
                available_partners,
                size=num_partners_to_find,
                p=partner_agent_probabilities[available_partners],
                replace=False
            )
            # --- END SIMPLIFIED PARTNER SELECTION WITH MIXING MATRIX ---

            for p in partners:
                if infected[p] or cancerous[p]:
                    continue

                base_beta = BETA_PER_AGENT[p]
                sus_multiplier = 1.0 - (efficacy if vaccinated[p] else 0.0)
                p_eff = base_beta * sus_multiplier
                if np.random.rand() < p_eff:
                    infected[p] = True
                    time_infected[p] = 0

        # --- CLEARANCE, TIME, AND CANCER PROGRESSION ---

        # 0. Screening detection: screened agents may be tested this timestep
        # and cleared with some effectiveness. Screening_prob_monthly is the
        # per-agent monthly probability of being screened/detected.
        currently_infected = np.where(infected)[0]
        if currently_infected.size > 0:
            detect_rand = np.random.rand(currently_infected.size)
            detect_mask = detect_rand < screening_prob_monthly[currently_infected]
            if detect_mask.any():
                detected_idxs = currently_infected[detect_mask]
                # treatment success
                treat_rand = np.random.rand(detected_idxs.size)
                cleared_by_screen = detected_idxs[treat_rand < screening_detection_effectiveness]
                if cleared_by_screen.size > 0:
                    infected[cleared_by_screen] = False
                    time_infected[cleared_by_screen] = 0

        # 1. Progression to Cancer *** ALSO could be fixed and calibrated*** rn using just case/total pop for probability but could be improved
        progression_pool = np.where(infected & ~cancerous)[0]

        if progression_pool.size > 0:
            duration_factor = np.log(time_infected[progression_pool] + 1)
            # The base rate is multiplied by the log factor
            progression_prob = PROB_PROGRESS_TO_CANCER * duration_factor

            rand_vals = np.random.rand(progression_pool.size)
            new_cancer_cases = progression_pool[rand_vals < progression_prob]

            if new_cancer_cases.size > 0:
                cancerous[new_cancer_cases] = True
                infected[new_cancer_cases] = False
                time_infected[new_cancer_cases] = 0

        # 2. Viral Clearance
        currently_infected = np.where(infected)[0]
        if currently_infected.size > 0:
            clear_probs = clearance_per_agent[currently_infected]
            rand_vals = np.random.rand(currently_infected.size)
            cleared = rand_vals < clear_probs

            infected[currently_infected[cleared]] = False
            time_infected[currently_infected[cleared]] = 0

        # 3. Update Time Infected
        time_infected[infected] += 1

        # --- METRIC CALCULATION ---
        prevalence.append(float(infected.mean()))
        cancer_incidence.append(float(cancerous.mean()))

    return SimulationResults(
        prevalence=prevalence,
        cancer_incidence=cancer_incidence
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run HPV ABM with Walker et al. Revised Mixing Matrix")
    parser.add_argument(
        "--coverage", type=float, default=0.4, help="Vaccination coverage (0-1)"
    )
    parser.add_argument(
        "--efficacy", type=float, default=0.9, help="Vaccine efficacy (0-1)"
    )
    parser.add_argument(
        "--timesteps", type=int, default=300, help="Number of timesteps (300 months = 25 years)"
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--beta-scale",
        dest="BETA_SCALING_FACTOR_C",
        type=float,
        default=0.0001, # New default reflecting the much higher transmission capacity of the heterogeneous network
        help="Calibration scaling factor (C) for transmission beta.",
    )
    parser.add_argument(
        "--N",
        "--population-size",
        dest="population_size",
        type=int,
        default=2000,
        help="Population size (increased to 2000 to better reflect mixing classes)",
    )
    return parser.parse_args()


def main():
    args = _parse_args()

    # NOTE: You MUST increase the timesteps and population size for this model
    # to yield stable results reflecting cancer latency.
    if args.timesteps < 240:
        print(f"Warning: Timesteps set to {args.timesteps}. Recommended minimum is 240 (20 years) for cancer modeling.")
    if args.population_size < 1500:
        print(f"Warning: Population size set to {args.population_size}. Recommended minimum is 1500 to fill all 20 subgroups.")

    results = run_simulation(
        coverage=args.coverage,
        efficacy=args.efficacy,
        steps=args.timesteps,
        population_size=args.population_size,
        seed=args.seed,
        BETA_SCALING_FACTOR_C=args.BETA_SCALING_FACTOR_C,
    )

    results.plot_prevalence_curve(show=True)

    # Write basic run metadata for reproducibility
    meta = {
        "entrypoint": "simulate_hpvm_revised.py",
        "coverage": args.coverage,
        "efficacy": args.efficacy,
        "timesteps": args.timesteps,
        "population_size": args.population_size,
        "seed": args.seed,
        # Model Parameters
        "BETA_SCALING_FACTOR_C": args.BETA_SCALING_FACTOR_C,
        "PROB_PROGRESS_TO_CANCER": PROB_PROGRESS_TO_CANCER,
        "EPSILON_A": EPSILON_A,
        "EPSILON_S": EPSILON_S,
        "SEX_ACTIVITY_DISTRIBUTION": SEX_ACTIVITY_DISTRIBUTION.tolist(),
        "partners_per_step": 2, # Baseline partner count
    }
    write_run_metadata(meta, out_path=None)


if __name__ == "__main__":
    main()