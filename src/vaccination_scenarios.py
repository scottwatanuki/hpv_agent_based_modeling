from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from simulate_hpvm_abm import run_simulation
from analysis_utils import write_prevalence_outputs, write_run_metadata


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"


@dataclass
class Scenario:
    name: str
    coverage: float
    efficacy: float
    timesteps: int
    seed: Optional[int] = None
    population_size: int = 1000
    # optional per-race mappings
    coverage_by_race: Optional[Dict[str, float]] = None
    screening_by_race: Optional[Dict[str, float]] = None

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "Scenario":
        return Scenario(
            name=str(d.get("name") or d.get("id") or "scenario"),
            coverage=float(d.get("coverage", 0.7)),
            efficacy=float(d.get("efficacy", 0.9)),
            timesteps=int(d.get("timesteps", 120)),
            seed=(None if d.get("seed") is None else int(d.get("seed"))),
            population_size=int(d.get("population_size", 1000)),
            coverage_by_race=d.get("coverage_by_race"),
            screening_by_race=d.get("screening_by_race"),
        )


def load_scenarios(path: Path) -> List[Scenario]:
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict) and "scenarios" in data:
        items = data["scenarios"]
    elif isinstance(data, list):
        items = data
    else:
        raise ValueError("Scenario JSON must be a list or a dict with key 'scenarios'.")
    return [Scenario.from_dict(x) for x in items]


def run_scenarios(scenarios: List[Scenario]) -> pd.DataFrame:
    rows = []
    for sc in scenarios:
        res = run_simulation(
            coverage=sc.coverage,
            efficacy=sc.efficacy,
            steps=sc.timesteps,
            population_size=sc.population_size,
            coverage_by_race=sc.coverage_by_race,
            screening_by_race=sc.screening_by_race,
            seed=sc.seed,
        )
        paths = write_prevalence_outputs(res.prevalence, prefix=sc.name)
        # write per-scenario metadata file
        write_run_metadata(
            {
                "entrypoint": "vaccination_scenarios.py",
                "name": sc.name,
                "coverage": sc.coverage,
                "efficacy": sc.efficacy,
                "timesteps": sc.timesteps,
                "population_size": sc.population_size,
                "seed": sc.seed,
                # defaults used by run_simulation not overridden here
                "beta": 0.02,
                "clearance_prob": 0.05,
                "partners_per_step": 5,
            },
            prefix=sc.name,
        )
        rows.append(
            {
                "name": sc.name,
                "coverage": sc.coverage,
                "efficacy": sc.efficacy,
                "timesteps": sc.timesteps,
                "population_size": sc.population_size,
                "seed": sc.seed,
                "figure": paths["figure"],
                "summary_csv": paths["summary_csv"],
            }
        )
    df = pd.DataFrame(rows)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = RESULTS_DIR / "scenario_runs.csv"
    df.to_csv(out_csv, index=False)
    return df


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run vaccination scenarios from a JSON config"
    )
    p.add_argument("--scenarios", required=True, help="Path to scenarios JSON file")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    sc_path = Path(args.scenarios)
    scenarios = load_scenarios(sc_path)
    df = run_scenarios(scenarios)
    print(df)


if __name__ == "__main__":
    main()
