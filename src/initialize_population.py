from __future__ import annotations
import pandas as pd


def load_race_distribution_from_csv(csv_path: str) -> dict:
    """load race distribution from a csv with columns 'Race/Ethnicity Category' and 'Population Estimate'.
    returns a dict: race -> proportion (excluding 'Total Population').
    """
    df = pd.read_csv(csv_path)
    df = df[df['Race/Ethnicity Category'] != 'Total Population']
    total = df['Population Estimate'].sum()
    race_dist = {row['Race/Ethnicity Category']: row['Population Estimate'] / total for _, row in df.iterrows()}
    return race_dist


def build_vaccination_coverage_map() -> dict:
    """build comprehensive vaccination coverage map by race, sex, and age band.
    
    based on:
    - ages 13-17: 45.6% (all races/sexes)
    - ages 19-21: females 58.9%, males 39.8%
    - ages 22-26: females 57.0%, males 32.1%
    - ages 27-44: race-specific rates from pmc.ncbi.nlm.nih.gov/articles/PMC10984122/
    - all other ages: 0%
    
    returns dict with keys (race, sex, age_group_label) -> coverage (0-1)
    """
    coverage_map = {}

    races = [
        "White Alone",
        "Black or African American Alone",
        "American Indian and Alaska Native Alone",
        "Asian Alone",
        "Native Hawaiian and Other Pacific Islander Alone",
        "Some Other Race Alone",
        "Two or More Races",
        "Hispanic or Latino (Any Race)"
    ]

    sexes = ["Male", "Female"]

    for race in races:
        for sex in sexes:
            coverage_map[(race, sex, "10 to 14 years")] = 0.456
            coverage_map[(race, sex, "15 to 19 years")] = 0.456

    for race in races:
        coverage_map[(race, "Female", "20 to 24 years")] = 0.589
        coverage_map[(race, "Male", "20 to 24 years")] = 0.398

    for race in races:
        coverage_map[(race, "Female", "25 to 29 years")] = 0.570
        coverage_map[(race, "Male", "25 to 29 years")] = 0.321

    race_specific_27_44 = {
        "White Alone": 0.159,
        "Black or African American Alone": 0.194,
        "Hispanic or Latino (Any Race)": 0.119,
        "American Indian and Alaska Native Alone": 0.154,
        "Asian Alone": 0.154,
        "Native Hawaiian and Other Pacific Islander Alone": 0.154,
        "Some Other Race Alone": 0.154,
        "Two or More Races": 0.154
    }

    age_bands_27_44 = ["25 to 29 years", "30 to 34 years", "35 to 39 years", "40 to 44 years"]
    for race in races:
        rate = race_specific_27_44.get(race, 0.154)
        for sex in sexes:
            for age_band in age_bands_27_44:
                coverage_map[(race, sex, age_band)] = rate

    other_bands = [
        "Under 5 years", "5 to 9 years",
        "45 to 49 years", "50 to 54 years", "55 to 59 years",
        "60 to 64 years", "65 to 69 years", "70 to 74 years",
        "75 to 79 years", "80 to 84 years", "85 years and over"
    ]
    for race in races:
        for sex in sexes:
            for age_band in other_bands:
                coverage_map[(race, sex, age_band)] = 0.0

    return coverage_map


def build_screening_rates_by_race() -> dict:
    """build cervical cancer screening rates by race for women 25+.
    
    annual screening rates (women 25+):
    - white: 63.5%
    - black: 53.2%
    - hispanic: 65.4%
    
    returns dict mapping race -> annual screening rate (0-1)
    """
    screening_map = {
        "White Alone": 0.635,
        "Black or African American Alone": 0.532,
        "Hispanic or Latino (Any Race)": 0.654,
        "American Indian and Alaska Native Alone": 0.607,
        "Asian Alone": 0.607,
        "Native Hawaiian and Other Pacific Islander Alone": 0.607,
        "Some Other Race Alone": 0.607,
        "Two or More Races": 0.607
    }
    return screening_map
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

    df: pd.DataFrame


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
        "85 years and over": (85, 95),
    }
    return mapping.get(label, (0, 95))

def calibrate():
    df = pd.read_csv("data/raw/NSFG_2022_2023_FemRespPUFData.csv")
    df1 = df[["AGE_R","P1YHSAGE"]].dropna()
    df1['age_group'] = pd.cut(df1['AGE_R'], bins=[0,15,25,35,45,99], labels=['0-14', '15-24', '25-34', '35-44', '45+'], right=False)
    df1['partners'] = pd.cut(df1['P1YHSAGE'], bins=[0,15,25,35,45,99], labels=['0-14', '15-24', '25-34', '35-44', '45+'], right=False)
    print(pd.crosstab(df1['age_group'], df1['partners']))
    df1 = pd.read_csv("data/raw/NSFG_2022_2023_MaleRespPUFData.csv")
    df2 = df1[["OPPYEARNUM", "SAMYEARNUM"]].dropna(subset=['OPPYEARNUM', 'SAMYEARNUM'], how='all')
    df2 = df2.fillna(0)
    df2['sum'] = df2['OPPYEARNUM'] + df2['SAMYEARNUM']
    df2 = df2.drop(columns=['OPPYEARNUM', 'SAMYEARNUM'])
    df = df[["OPPYEARNUM", "SAMYEARNUM"]].dropna(subset=['OPPYEARNUM', 'SAMYEARNUM'], how='all')
    df = df.fillna(0)
    df['sum'] = df['OPPYEARNUM'] + df['SAMYEARNUM']
    df = df.drop(columns=['OPPYEARNUM', 'SAMYEARNUM'])
    df = df['sum'].value_counts() + df2['sum'].value_counts()
    df = df.dropna()
    df = df.to_frame().reset_index()
    df = df[df['sum'] < 11]
    df['count_group'] = pd.cut(df['sum'], bins=[0,1,2,5,11], labels=['0','1','2-4', '5-10'], right=False)
    df = df.groupby('count_group')['count'].sum().reset_index()
    df['count'] = df['count']/df['count'].sum()
    print(df)

def _load_age_sex_distribution(csv_path: Optional[Path] = None) -> pd.DataFrame:
    """Load the age-sex distribution CSV.

    Expected columns: 'Age Group', 'Male_Estimate', 'Female_Estimate'.
    Returns a DataFrame with columns: age_group, sex, count.
    """
    path = csv_path or DATA_PATH
    df = pd.read_csv(path)
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
    coverage_by_race_and_ageband: Optional[Dict[tuple, float]] = None,
    coverage_by_race_sex_ageband: Optional[Dict[tuple, float]] = None,
) -> Population:
    """Create a synthetic population of size n using the ACS age-sex distribution.

    Vaccination is assigned via Bernoulli(coverage) independent of age/sex for this baseline.
    """
    if seed is not None:
        np.random.seed(seed)

    dist = _load_age_sex_distribution(csv_path)
    dist = dist[dist["count"] > 0].copy()
    dist["weight"] = dist["count"] / dist["count"].sum()

    strata = dist[["age_group", "sex", "weight"]].values
    choices = np.random.choice(len(strata), size=n, p=dist["weight"].values)

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
        if race_names is not None:
            race = str(np.random.choice(race_names, p=race_probs))
        else:
            race = "Unknown"
        vac_prob = coverage
        if coverage_by_race_sex_ageband is not None:
            key = (race, str(sex), str(age_group))
            if key in coverage_by_race_sex_ageband:
                vac_prob = float(coverage_by_race_sex_ageband[key])
        elif coverage_by_race_and_ageband is not None:
            key = (race, str(age_group))
            if key in coverage_by_race_and_ageband:
                vac_prob = float(coverage_by_race_and_ageband[key])
        elif coverage_by_race is not None and race in coverage_by_race:
            vac_prob = float(coverage_by_race[race])

        vaccinated = bool(np.random.rand() < vac_prob)
        screening_monthly = 0.0
        if str(sex).lower().startswith("f") and age >= 25:
            if screening_by_race is not None and race in screening_by_race:
                annual = float(screening_by_race[race])
                screening_monthly = annual / 12.0

        rows.append((idx, age, str(sex), str(age_group), vaccinated, race, screening_monthly))

    pop_df = pd.DataFrame(
        rows,
        columns=["id", "age", "sex", "age_group", "vaccinated", "race", "screening_prob_monthly"],
    )
    return Population(df=pop_df)

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
        return 0.0


__all__ = ["Population", "synthesize_population"]
