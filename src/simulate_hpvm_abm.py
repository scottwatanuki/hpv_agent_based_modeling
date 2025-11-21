import argparse
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
from src.analysis_utils import write_run_metadata
#adding in population
from src.initialize_population import synthesize_population, Population

import numpy as np
import matplotlib.pyplot as plt

# --- Age bands ---

import numpy as np

AGE_BANDS = [
    (0, 14),
    (15, 24),
    (25, 34),
    (35, 44),
    (45, 120),
]

def age_to_band_index(age: int) -> int:
    for i, (lo, hi) in enumerate(AGE_BANDS):
        if lo <= age <= hi:
            return i
    return len(AGE_BANDS) - 1

# --- Age-specific parameters (placeholders for now) ---
#thinking we could add in the neural network here to estimate what to put for all of these numbers?
BETA_BY_AGE_BAND = np.array([
    0.00,   # 0–14
    0.035,   # 15–24
    0.02,  # 25–34
    0.015,  # 35–44
    0.01,  # 45+
])

CLEARANCE_BY_AGE_BAND = np.array([
    0.00,   # 0–14 (placeholder)
    0.03,   # 15–24
    0.025,   # 25–34
    0.02,   # 35–44
    0.015,   # 45+
])


@dataclass
class SimulationResults:
    """Container for simulation outputs."""

    prevalence: list[float]

    def plot_prevalence_curve(self, show: bool = True, save_path: Optional[str] = None):
        plt.figure()
        plt.plot(self.prevalence, color="tab:blue", linewidth=2)
        plt.xlabel("Time step")
        plt.ylabel("Prevalence")
        plt.title("HPV-ABM Prevalence Over Time")
        plt.grid(True, alpha=0.3)
        if save_path:
            plt.savefig(save_path, dpi=200, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close()


def run_simulation(
    coverage: float = 0.5,
    efficacy: float = 0.9,
    steps: int = 120,
    *,
    population_size: int = 1000,
    partners_per_step: int = 5,
    seed: Optional[int] = None,
) -> SimulationResults:
    """
    Run an HPV transmission simulation with static random mixing per step.

    Parameters
    - coverage: vaccination coverage (0-1)
    - efficacy: vaccine efficacy (0-1), reduces susceptibility by (1 - efficacy)
    - steps: number of timesteps to simulate
    - population_size: number of agents
    - partners_per_step: number of random partners sampled per infected per step
    - seed: optional random seed for reproducibility

    Returns - simulation with infection over time steps
    """
    if seed is not None:
        np.random.seed(seed)
    #new - to reflect age groups
    pop: Population = synthesize_population(
        n=population_size,
        coverage=coverage,
        seed=seed,
    )
    pop_df = pop.df

    ages = pop_df["age"].to_numpy(dtype=int)
    vaccinated = pop_df["vaccinated"].to_numpy(dtype=bool)

    N = len(pop_df)
    # map each agent to an age band index
    age_band_idx = np.array([age_to_band_index(a) for a in ages], dtype=int)

    # per-agent parameters derived from age
    beta_per_agent = BETA_BY_AGE_BAND[age_band_idx]          # shape (N,)
    clearance_per_agent = CLEARANCE_BY_AGE_BAND[age_band_idx]  # shape (N,)

    # Initialize population vaccination and infection
    # vaccinated = np.random.rand(N) < coverage - changed to take out so its not random anymore
    infected = np.zeros(N, dtype=bool)

    # seed 10 initial infections changed this
    num_initial = 10
    initial_indices = np.random.choice(N, size=num_initial, replace=False)
    infected[initial_indices] = True
    prevalence = []

    for _ in range(steps):
        # iterate over a snapshot of infected indices to avoid within-step feedback
        currently_infected = np.where(infected)[0]
        for i in currently_infected:
            # sample partners without replacement, excluding self when possible
            # if population is small vs partners_per_step, allow replacement
            if N > partners_per_step:
                # ensure we don't select self
                candidates = np.delete(np.arange(N), i)
                partners = np.random.choice(
                    candidates, size=partners_per_step, replace=False
                )
            else:
                partners = np.random.choice(N, size=partners_per_step, replace=True)
            #changed to reflect age speficcs
            for p in partners:
                if infected[p]:
                    continue  # already infected

                # Susceptible partner's age-specific transmission probability
                base_beta = beta_per_agent[p]

                # Vaccination reduces susceptibility
                sus_multiplier = 1.0 - (efficacy if vaccinated[p] else 0.0)
                p_eff = base_beta * sus_multiplier
                if np.random.rand() < p_eff:
                    infected[p] = True

        currently_infected = np.where(infected)[0]
        if currently_infected.size > 0:
            clear_probs = clearance_per_agent[currently_infected]
            rand_vals = np.random.rand(currently_infected.size)
            cleared = rand_vals < clear_probs
            infected[currently_infected[cleared]] = False

        prevalence.append(float(infected.mean()))

    return SimulationResults(prevalence=prevalence)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run HPV ABM baseline simulation")
    parser.add_argument(
        "--coverage", type=float, default=0.7, help="Vaccination coverage (0-1)"
    )
    parser.add_argument(
        "--efficacy", type=float, default=0.9, help="Vaccine efficacy (0-1)"
    )
    parser.add_argument(
        "--timesteps", type=int, default=120, help="Number of timesteps"
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--N",
        "--population-size",
        dest="population_size",
        type=int,
        default=1000,
        help="Population size (default 1000)",
    )
    return parser.parse_args()


def main():
    args = _parse_args()
    results = run_simulation(
        coverage=args.coverage,
        efficacy=args.efficacy,
        steps=args.timesteps,
        population_size=args.population_size,
        seed=args.seed,
    )

    # Plot to screen for now; saving will be handled in analysis utils in later tasks
    results.plot_prevalence_curve(show=True)

    # Write basic run metadata for reproducibility
    meta = {
        "entrypoint": "simulate_hpvm_abm.py",
        "coverage": args.coverage,
        "efficacy": args.efficacy,
        "timesteps": args.timesteps,
        "population_size": args.population_size,
        "seed": args.seed,
        # defaults used inside run_simulation
        "beta_by_age_band": BETA_BY_AGE_BAND.tolist(),
        "clearance_by_age_band": CLEARANCE_BY_AGE_BAND.tolist(),
        "partners_per_step": 5,
    }
    write_run_metadata(meta, out_path=None)


if __name__ == "__main__":
    main()
