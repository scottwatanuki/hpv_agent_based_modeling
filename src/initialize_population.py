from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict

import numpy as np
import pandas as pd


DATA_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "processed"
    / "atl_ga_age_sex_cleaned.csv"
)


@dataclass
class Population:
    """Container for a synthetic population."""

    df: pd.DataFrame  # columns: id, age, sex, age_group, vaccinated


def _parse_age_group_to_range(age_group: str) -> Tuple[int, int]:
    """Map ACS-like age group labels to inclusive integer age ranges.

    Falls back to broad guesses if unexpected labels are encountered.
    """
    label = age_group.strip().lower()
    mapping = {
        "under 5 years": (0, 4),
        "5 to 9 years": (5, 9),
        "10 to 14 years": (10, 14),
        "15 to 19 years": (15, 19),
        "20 to 24 years": (20, 24),
        "25 to 29 years": (25, 29),
        "30 to 34 years": (30, 34),
        "35 to 39 years": (35, 39),
        "40 to 44 years": (40, 44),
        "45 to 49 years": (45, 49),
        "50 to 54 years": (50, 54),
        "55 to 59 years": (55, 59),
        "60 to 64 years": (60, 64),
        "65 to 69 years": (65, 69),
        "70 to 74 years": (70, 74),
        "75 to 79 years": (75, 79),
        "80 to 84 years": (80, 84),
        "85 years and over": (85, 95),  # cap upper bound for sampling
    }
    return mapping.get(label, (0, 95))


def _load_age_sex_distribution(csv_path: Optional[Path] = None) -> pd.DataFrame:
    """Load the age-sex distribution CSV.

    Expected columns: 'Age Group', 'Male_Estimate', 'Female_Estimate'.
    Returns a DataFrame with columns: age_group, sex, count.
    """
    path = csv_path or DATA_PATH
    df = pd.read_csv(path)
    # melt into long format
    long_df = df.melt(
        id_vars=["Age Group"],
        value_vars=["Male_Estimate", "Female_Estimate"],
        var_name="sex",
        value_name="count",
    )
    long_df["sex"] = (
        long_df["sex"]
        .str.replace("_Estimate", "", regex=False)
        .str.lower()
        .str.replace("male", "Male")
        .str.replace("female", "Female")
    )
    long_df = long_df.rename(columns={"Age Group": "age_group"})
    long_df = long_df[["age_group", "sex", "count"]]
    return long_df


def synthesize_population(
    n: int,
    coverage: float,
    *,
    csv_path: Optional[Path] = None,
    seed: Optional[int] = None,
    # Optional per-race inputs. If provided, `race_dist` is a mapping
    # race -> proportion (must sum to 1 or be normalized). `coverage_by_race`
    # maps race -> vaccination coverage (0-1). `screening_by_race` maps
    # race -> annual screening uptake (0-1). All are optional and
    # backward-compatible with existing calls.
    race_dist: Optional[Dict[str, float]] = None,
    coverage_by_race: Optional[Dict[str, float]] = None,
    screening_by_race: Optional[Dict[str, float]] = None,
) -> Population:
    """Create a synthetic population of size n using the ACS age-sex distribution.

    Vaccination is assigned via Bernoulli(coverage) independent of age/sex for this baseline.
    """
    if seed is not None:
        np.random.seed(seed)

    dist = _load_age_sex_distribution(csv_path)
    dist = dist[dist["count"] > 0].copy()
    dist["weight"] = dist["count"] / dist["count"].sum()

    # sample strata indices according to weights
    strata = dist[["age_group", "sex", "weight"]].values
    choices = np.random.choice(len(strata), size=n, p=dist["weight"].values)

    # Prepare race sampling if provided
    race_names = None
    race_probs = None
    if race_dist is not None:
        race_items = list(race_dist.items())
        race_names = [r for r, _ in race_items]
        probs = np.array([p for _, p in race_items], dtype=float)
        if probs.sum() > 0:
            probs = probs / probs.sum()
        else:
            probs = np.ones_like(probs) / len(probs)
        race_probs = probs

    rows = []
    for idx, stratum_idx in enumerate(choices):
        age_group, sex, _w = strata[stratum_idx]
        lo, hi = _parse_age_group_to_range(str(age_group))
        age = int(np.random.randint(lo, hi + 1))
        # assign race (if race_dist provided)
        if race_names is not None:
            race = str(np.random.choice(race_names, p=race_probs))
        else:
            race = "Unknown"

        # determine vaccination probability for this agent (race-specific override)
        if coverage_by_race is not None and race in coverage_by_race:
            vac_prob = float(coverage_by_race[race])
        else:
            vac_prob = coverage

        vaccinated = bool(np.random.rand() < vac_prob)

        # screening uptake given as annual probability by race -> convert to monthly
        screening_monthly = 0.0
        if screening_by_race is not None and race in screening_by_race:
            annual = float(screening_by_race[race])
            # convert to monthly (approx)
            screening_monthly = annual / 12.0

        rows.append((idx, age, str(sex), str(age_group), vaccinated, race, screening_monthly))

    pop_df = pd.DataFrame(
        rows,
        columns=["id", "age", "sex", "age_group", "vaccinated", "race", "screening_prob_monthly"],
    )
    return Population(df=pop_df)

#add into function above
AGE_GROUP_COVERAGE = {
    "13-17": 0.615,
    "18-26": 0.516,
}

def _coverage_for_age(age: int) -> float:
    if 13 <= age <= 17:
        return AGE_GROUP_COVERAGE["13-17"]
    elif 18 <= age <= 26:
        return AGE_GROUP_COVERAGE["18-26"]
    else:
        return 0.0  # or some base coverage


__all__ = ["Population", "synthesize_population"]
