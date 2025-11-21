
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Dict

import numpy as np
import pandas as pd

from .initialize_population import Population


@dataclass #need to be calibrated
class SexualNetworkParams:
    avg_degree: int = 5
    p_within_band: float = 0.8
    seed: Optional[int] = None


def _age_band(age: int) -> int:
    if age < 15:
        return 0
    elif age < 25:
        return 1
    elif age < 35:
        return 2
    elif age < 45:
        return 3
    else:
        return 4


def build_age_mixing_network(
    population: Population,
    params: SexualNetworkParams,
) -> List[np.ndarray]:
    rng = np.random.default_rng(params.seed)

    df: pd.DataFrame = population.df
    ages = df["age"].to_numpy(dtype=int)
    N = ages.size

    band_ids = np.array([_age_band(a) for a in ages], dtype=int)

    band_to_indices: Dict[int, np.ndarray] = {}
    for idx, b in enumerate(band_ids):
        band_to_indices.setdefault(b, []).append(idx)
    band_to_indices = {b: np.array(idxs, dtype=int) for b, idxs in band_to_indices.items()}

    neighbors = [set() for _ in range(N)]
    target_edges = N * params.avg_degree // 2

    edges = 0
    attempts = 0
    max_attempts = target_edges * 20

    while edges < target_edges and attempts < max_attempts:
        attempts += 1

        i = int(rng.integers(0, N))
        bi = band_ids[i]

        if rng.random() < params.p_within_band:
            candidates = band_to_indices[bi]
        else:
            other = [idxs for b, idxs in band_to_indices.items() if b != bi]
            if not other:
                continue
            candidates = np.concatenate(other)

        if candidates.size == 0:
            continue

        j = int(candidates[rng.integers(0, candidates.size)])
        if i == j:
            continue
        if j in neighbors[i]:
            continue

        neighbors[i].add(j)
        neighbors[j].add(i)
        edges += 1

    return [np.fromiter(neigh, dtype=int) for neigh in neighbors]
