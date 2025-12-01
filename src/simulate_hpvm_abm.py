import argparse
from dataclasses import dataclass
from typing import Optional, Dict
from pathlib import Path

from src.analysis_utils import write_run_metadata
from src.initialize_population import synthesize_population, Population, load_race_distribution_from_csv, build_vaccination_coverage_map, build_screening_rates_by_race

import numpy as np
import matplotlib.pyplot as plt

AGE_BANDS = [
    (0, 14),   
    (15, 24),  
    (25, 34),  
    (35, 44),  
    (45, 100), 
]
NUM_AGE_BANDS = len(AGE_BANDS)

SEX_ACTIVITY_CLASSES = ["Low", "Mod", "Mod-High", "High"]
NUM_SEX_CLASSES = len(SEX_ACTIVITY_CLASSES)

SEX_ACTIVITY_DISTRIBUTION = np.array([0.54, 0.23, 0.17, 0.06])
SEX_ACTIVITY_DISTRIBUTION = SEX_ACTIVITY_DISTRIBUTION / SEX_ACTIVITY_DISTRIBUTION.sum()


def age_to_band_index(age: int) -> int:
    """Maps agent's age to the correct age band index."""
    for i, (lo, hi) in enumerate(AGE_BANDS):
        if lo <= age <= hi:
            return i
    return NUM_AGE_BANDS - 1

RELATIVE_CONTACT_BY_AGE_BAND = np.array([0.00, 1.00, 0.70, 0.40, 0.20])

CLEARANCE_BY_AGE_BAND = np.array([0.00, 0.008, 0.006, 0.005, 0.003])

PROB_PROGRESS_TO_CANCER = 0.004167


EPSILON_A = 0.7
EPSILON_S = 0.3


AGE_MIXING_PREFERENCE = np.array([
    [1.0, 0.0, 0.0, 0.0, 0.0], # 0-14 (only mix with self, though not active)
    [1, 307, 56, 3, 8], # 15-24
    [2, 449, 580, 118, 46], # 25-34
    [1, 256, 493, 332, 114], # 35-44
    [0, 82, 143, 148, 114] # 45+
])

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
    P_subgroup = np.zeros((NUM_AGE_BANDS, NUM_SEX_CLASSES))
    for i in range(N):
        age_i = age_band_idx[i]
        class_h = sex_activity_idx[i]
        P_subgroup[age_i, class_h] += partnerships_offered_per_agent[i]

    sum_over_age = P_subgroup.sum(axis=0)
    sum_over_class = P_subgroup.sum(axis=1)
    sum_total = P_subgroup.sum()
    M = np.zeros((NUM_AGE_BANDS * NUM_SEX_CLASSES, NUM_AGE_BANDS * NUM_SEX_CLASSES))
    for i in range(NUM_AGE_BANDS): 
        for h in range(NUM_SEX_CLASSES): 
            source_idx = i * NUM_SEX_CLASSES + h 

            for j in range(NUM_AGE_BANDS):
                for m in range(NUM_SEX_CLASSES): 
                    target_idx = j * NUM_SEX_CLASSES + m 

                    rho = 0.0
                    P_jm = P_subgroup[j, m] # P_k'jm

                    term1 = EPSILON_A * EPSILON_S * (1 if i == j and h == m else 0.0)
                    rho += term1

                    term2 = 0.0
                    if i == j and sum_over_class[j] > 0:
                        term2 = EPSILON_A * (1.0 - EPSILON_S) * (P_jm / sum_over_class[j])
                    rho += term2

                    term3 = 0.0
                    if h == m and sum_over_age[m] > 0:
                        term3 = (1.0 - EPSILON_A) * EPSILON_S * (P_jm / sum_over_age[m])
                    rho += term3

                    term4 = 0.0
                    if sum_total > 0:
                        term4 = (1.0 - EPSILON_A) * (1.0 - EPSILON_S) * (P_jm / sum_total)
                    rho += term4

                    M[source_idx, target_idx] = rho
    row_sums = M.sum(axis=1)
    M[row_sums > 0] /= row_sums[row_sums > 0, np.newaxis]

    return M


@dataclass
class SimulationResults:
    """Container for simulation outputs."""

    prevalence: list[float]
    cancer_incidence: list[float]
    final_prevalence_by_age: Optional[Dict[str, float]] = None  # Final prevalence per age band
    final_prevalence_by_race: Optional[Dict[str, float]] = None  # Final prevalence per race

    def plot_prevalence_curve(self, show: bool = True, save_path: Optional[str] = None):
        plt.figure(figsize=(10, 6))
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
    partners_per_step: int = 1,
    BETA_SCALING_FACTOR_C: float = 0.05,
    seed: Optional[int] = None,

    #needs to be calibrated and looked at later
    # Optional per-race overrides: mapping race -> coverage (0-1) and race -> annual screening uptake (0-1)
    coverage_by_race: Optional[Dict[str, float]] = None,
    screening_by_race: Optional[Dict[str, float]] = None,
    # Probability that a screening event detects & clears infection when it occurs
    screening_detection_effectiveness: float = 0.9,
    # Initial prevalence (fraction of population infected at t=0). If None, uses default small seeding (10 individuals).
    initial_prevalence: Optional[float] = None,
    # Optional list of length NUM_AGE_BANDS giving initial prevalence per age band (fractions 0-1).
    initial_prevalence_by_age: Optional[list[float]] = None,
) -> SimulationResults:
    """
    Run an HPV transmission simulation with the revised, heterogeneous, age-activity mixing matrix.
    """
    if seed is not None:
        np.random.seed(seed)

    # synthesize_population is assumed to be able to handle these groups
    # load race distribution from csv and pass to synthesize_population
    race_dist = load_race_distribution_from_csv("data/processed/atl_ga_demographics_cleaned.csv")
    
    # load comprehensive vaccination coverage map by race, sex, and age band
    coverage_map = build_vaccination_coverage_map()
    
    # load screening rates by race (for women 25+)
    screening_rates = build_screening_rates_by_race()
    
    pop: Population = synthesize_population(
        n=population_size,
        coverage=coverage,
        seed=seed,
        race_dist=race_dist,
        coverage_by_race_sex_ageband=coverage_map,
        coverage_by_race=coverage_by_race,
        screening_by_race=screening_rates if screening_by_race is None else screening_by_race,
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

    # new: assign sexual activity class based on the defined distribution
    sex_activity_idx = np.random.choice(
        NUM_SEX_CLASSES,
        size=N,
        p=SEX_ACTIVITY_DISTRIBUTION
    )

    # avg number of partners per step is now heterogeneous based on activity class
    # calibrated to NSFG data: mean ~1.5 partners/year (0.125/month) across all sexually active adults
    ACTIVITY_PARTNER_MULTIPLIER = np.array([0.9, 1.2, 1.8, 3.0])

    # partnerships offered by the agent in this timestep, dependent on their activity class.
    partnerships_offered_per_agent = ACTIVITY_PARTNER_MULTIPLIER[sex_activity_idx]

    # this calculation is run once before the simulation loop.
    REVISED_MIXING_MATRIX = calculate_mixing_matrix(
        N,
        age_band_idx,
        sex_activity_idx,
        partnerships_offered_per_agent
    )

    # The base transmission rate is now scaled by the ACTIVITY_PARTNER_MULTIPLIER (heterogeneity)
    # The overall beta (BETA_PER_AGENT) is the calibration factor * age effect * activity effect
    BETA_PER_AGENT = BETA_SCALING_FACTOR_C * RELATIVE_CONTACT_BY_AGE_BAND[age_band_idx]
    clearance_per_agent = CLEARANCE_BY_AGE_BAND[age_band_idx]

    # state variables
    infected = np.zeros(N, dtype=bool)
    time_infected = np.zeros(N, dtype=int)
    cancerous = np.zeros(N, dtype=bool)

    if initial_prevalence_by_age is not None:
        probs = np.array(initial_prevalence_by_age, dtype=float)
        if probs.shape[0] != NUM_AGE_BANDS:
            raise ValueError(f"initial_prevalence_by_age must have length {NUM_AGE_BANDS}")
        # assign infection per-agent with Bernoulli draw based on their age-band probability
        draw = np.random.rand(N)
        agent_probs = probs[age_band_idx]
        infected = draw < agent_probs
        time_infected[infected] = 1
    else:
        if initial_prevalence is not None:
            # translate prevalence fraction to integer initial infected count
            num_initial = max(1, int(np.round(initial_prevalence * N)))
        else:
            num_initial = 10

        initial_indices = np.random.choice(N, size=num_initial, replace=False)
        infected[initial_indices] = True
        time_infected[initial_indices] = 1

    prevalence = []
    cancer_incidence = []
    
    # track actual partner counts per agent per timestep
    partner_counts_per_agent = np.zeros(N, dtype=int)

    # detect sex/gender column if available so we can restrict contacts to opposite sex
    sex_col = None
    if "sex" in pop_df.columns:
        sex_col = pop_df["sex"].to_numpy()
    elif "gender" in pop_df.columns:
        sex_col = pop_df["gender"].to_numpy()

    # normalize sex values to lowercase strings when present
    if sex_col is not None:
        # convert bytes/objects to lowercase strings for robust comparison
        sex_col = np.array([str(x).strip().lower() for x in sex_col])
        # precompute boolean male indicator for fast opposite-sex masks
        is_male = np.array([s.startswith("m") for s in sex_col], dtype=bool)
    else:
        is_male = None

    # precompute each agent's flattened subgroup index (age x activity)
    flat_idx_per_agent = age_band_idx * NUM_SEX_CLASSES + sex_activity_idx

    for t in range(steps):

        # partner selection and transmission
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

            # begin partner selection with mixing matrix

            # For simplicity and to integrate the matrix, we'll revert to finding agents
            # whose characteristics match the mix.

            # Get the expected number of partners this agent will form this step
            expected_partners = partners_per_step * (
                ACTIVITY_PARTNER_MULTIPLIER[source_class_idx] / ACTIVITY_PARTNER_MULTIPLIER.mean()
            )
            num_partners_to_find = int(np.round(expected_partners))

            if num_partners_to_find == 0:
                 continue

            # we create a probability vector over the entire population (n) based on the mixing matrix
            # vectorized: prob(j) = m[source_flat_idx, flat_idx_per_agent[j]]
            partner_agent_probabilities = REVISED_MIXING_MATRIX[source_flat_idx, flat_idx_per_agent].copy()

            # apply opposite-sex constraint and exclude self
            mask = np.ones(N, dtype=bool)
            mask[i] = False
            if is_male is not None:
                mask &= (is_male != is_male[i])
            partner_agent_probabilities[~mask] = 0.0

            total_p = partner_agent_probabilities.sum()
            if total_p == 0.0:
                continue

            # normalize and sample from the entire population
            partner_agent_probabilities /= total_p


            available_partners = np.where(partner_agent_probabilities > 0)[0]
            num_partners_to_find = min(num_partners_to_find, available_partners.size)

            partners = np.random.choice(
                available_partners,
                size=num_partners_to_find,
                p=partner_agent_probabilities[available_partners],
                replace=False
            )

            # track actual partnerships formed
            partner_counts_per_agent[i] += len(partners)
            for p in partners:
                partner_counts_per_agent[p] += 1
                
                if infected[p] or cancerous[p]:
                    continue

                base_beta = BETA_PER_AGENT[p]
                sus_multiplier = 1.0 - (efficacy if vaccinated[p] else 0.0)
                p_eff = base_beta * sus_multiplier
                if np.random.rand() < p_eff:
                    infected[p] = True
                    time_infected[p] = 0

        # screening: some infected are screened and cleared this month (women 25+ have nonzero rates)
        if screening_prob_monthly is not None and screening_prob_monthly.size == N:
            # draw who gets screened this month
            got_screened = np.random.rand(N) < screening_prob_monthly
            if np.any(got_screened):
                # among screened, infections are detected/cleared with this effectiveness
                screened_infected = got_screened & infected
                if np.any(screened_infected):
                    detected = np.random.rand(screened_infected.sum()) < screening_detection_effectiveness
                    if detected.any():
                        idxs = np.where(screened_infected)[0][detected]
                        infected[idxs] = False
                        time_infected[idxs] = 0

        # clearance, progression, and cancer progression

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

        prevalence.append(float(infected.mean()))
        cancer_incidence.append(float(cancerous.mean()))

    avg_partners_per_month = partner_counts_per_agent / steps
    print(f"\n=== Partner Statistics ===")
    print(f"Mean partners/month: {avg_partners_per_month.mean():.2f}")
    print(f"Median partners/month: {np.median(avg_partners_per_month):.2f}")
    print(f"Min: {avg_partners_per_month.min():.2f}, Max: {avg_partners_per_month.max():.2f}")
    print(f"25th percentile: {np.percentile(avg_partners_per_month, 25):.2f}")
    print(f"75th percentile: {np.percentile(avg_partners_per_month, 75):.2f}")
    print(f"90th percentile: {np.percentile(avg_partners_per_month, 90):.2f}")
    print(f"95th percentile: {np.percentile(avg_partners_per_month, 95):.2f}")

    # Compute final prevalence by age group and by race
    final_prev_by_age = {}
    for age_idx in range(NUM_AGE_BANDS):
        age_lo, age_hi = AGE_BANDS[age_idx]
        mask = (age_band_idx == age_idx)
        if mask.sum() > 0:
            final_prev_by_age[f"{age_lo}-{age_hi}"] = float(infected[mask].mean())

    final_prev_by_race = {}
    if "race" in pop_df.columns:
        for race in pop_df["race"].unique():
            mask = (pop_df["race"] == race).to_numpy()
            if mask.sum() > 0:
                final_prev_by_race[race] = float(infected[mask].mean())

    return SimulationResults(
        prevalence=prevalence,
        cancer_incidence=cancer_incidence,
        final_prevalence_by_age=final_prev_by_age if final_prev_by_age else None,
        final_prevalence_by_race=final_prev_by_race if final_prev_by_race else None
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
        "--initial-prevalence",
        dest="initial_prevalence",
        type=float,
        default=None,
        help="Initial infected fraction of the population (0-1). If set, seeds that fraction at t=0.",
    )
    parser.add_argument(
        "--initial-prevalence-by-age",
        dest="initial_prevalence_by_age",
        type=str,
        default=None,
        help=("Comma-separated list of initial prevalence fractions for each age band "
              f"(length {NUM_AGE_BANDS}). Example: '0.0,0.129,0.185,0.06,0.03'"),
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
        initial_prevalence=args.initial_prevalence,
        initial_prevalence_by_age=(None if args.initial_prevalence_by_age is None else [float(x) for x in args.initial_prevalence_by_age.split(",")]),
    )

    results.plot_prevalence_curve(show=True)

    # Print final prevalence by age group
    if results.final_prevalence_by_age:
        print("\n=== Final HPV Prevalence by Age Group ===")
        for age_group, prev in results.final_prevalence_by_age.items():
            print(f"  Age {age_group}: {prev:.4f} ({prev*100:.2f}%)")

    # Print final prevalence by race
    if results.final_prevalence_by_race:
        print("\n=== Final HPV Prevalence by Race ===")
        for race, prev in results.final_prevalence_by_race.items():
            print(f"  {race}: {prev:.4f} ({prev*100:.2f}%)")

    print(f"\n=== Overall Final Prevalence ===")
    print(f"  {results.prevalence[-1]:.4f} ({results.prevalence[-1]*100:.2f}%)")

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
        "partners_per_step": 1,
        "ACTIVITY_PARTNER_MULTIPLIER": [0.9, 1.2, 1.8, 3.0],
    }
    write_run_metadata(meta, out_path=None)


if __name__ == "__main__":
    main()