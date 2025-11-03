import argparse
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
from src.analysis_utils import write_run_metadata

import numpy as np
import matplotlib.pyplot as plt


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
    coverage: float = 0.7,
    efficacy: float = 0.9,
    steps: int = 120,
    *,
    population_size: int = 1000,
    beta: float = 0.02,
    clearance_prob: float = 0.05,
    partners_per_step: int = 5,
    seed: Optional[int] = None,
) -> SimulationResults:
    """
    Run a simple HPV transmission simulation with static random mixing per step.

    Parameters
    - coverage: vaccination coverage (0-1)
    - efficacy: vaccine efficacy (0-1), reduces susceptibility by (1 - efficacy)
    - steps: number of timesteps to simulate
    - population_size: number of agents
    - beta: per-contact transmission probability
    - clearance_prob: per-step clearance probability for infected individuals
    - partners_per_step: number of random partners sampled per infected per step
    - seed: optional random seed for reproducibility

    Returns
    - SimulationResults with prevalence over time
    """
    if seed is not None:
        np.random.seed(seed)

    N = population_size

    # Initialize population vaccination and infection
    vaccinated = np.random.rand(N) < coverage
    infected = np.zeros(N, dtype=bool)
    # seed one initial infection
    infected[np.random.randint(0, N)] = True

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

            for p in partners:
                if not infected[p]:
                    susceptibility = 1.0 - (efficacy if vaccinated[p] else 0.0)
                    p_eff = beta * susceptibility
                    if np.random.rand() < p_eff:
                        infected[p] = True

            # Clearance for the infectious individual at the end of contacts
            if np.random.rand() < clearance_prob:
                infected[i] = False

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
        "beta": 0.02,
        "clearance_prob": 0.05,
        "partners_per_step": 5,
    }
    write_run_metadata(meta, out_path=None)


if __name__ == "__main__":
    main()
