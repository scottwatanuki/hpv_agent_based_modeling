from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence, Optional, Dict, Any
from datetime import datetime, timezone
import json
import platform
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"


def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def prevalence_from_states(states: Sequence[np.ndarray]) -> List[float]:
    """Compute prevalence over time from a sequence of boolean infection states.

    states: sequence of arrays of shape (N,) with dtype=bool for each timestep.
    returns: list of prevalence values in [0,1].
    """
    if not states:
        return []
    return [float(s.mean()) for s in states]


def incidence_from_states(states: Sequence[np.ndarray]) -> List[float]:
    """Compute incidence (new infections per timestep per capita) from states.

    Assumes states[t] and states[t-1] are boolean arrays of same length N.
    Returns list where index 0 is 0.0 by definition (no prior state).
    """
    if not states:
        return []
    N = states[0].size
    inc = [0.0]
    for t in range(1, len(states)):
        prev = states[t - 1]
        curr = states[t]
        new_inf = np.logical_and(curr, np.logical_not(prev)).sum()
        inc.append(float(new_inf) / float(N))
    return inc


def plot_and_save_prevalence(
    prevalence: Sequence[float], save_path: Path, *, show: bool = False, dpi: int = 300
) -> None:
    _ensure_parent_dir(save_path)
    plt.figure()
    plt.plot(prevalence, color="tab:blue", linewidth=2)
    plt.xlabel("Time step")
    plt.ylabel("Prevalence")
    plt.title("HPV-ABM Prevalence Over Time")
    plt.grid(True, alpha=0.3)
    plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def save_summary_stats(prevalence: Sequence[float], save_path: Path) -> None:
    """Save simple summary statistics of the prevalence curve to CSV (one row)."""
    _ensure_parent_dir(save_path)
    if not prevalence:
        df = pd.DataFrame(
            [
                {
                    "timesteps": 0,
                    "final_prevalence": np.nan,
                    "mean_prevalence": np.nan,
                    "peak_prevalence": np.nan,
                    "time_of_peak": np.nan,
                }
            ]
        )
        df.to_csv(save_path, index=False)
        return

    arr = np.asarray(prevalence, dtype=float)
    peak_idx = int(np.nanargmax(arr))
    row = {
        "timesteps": int(len(arr)),
        "final_prevalence": float(arr[-1]),
        "mean_prevalence": float(np.nanmean(arr)),
        "peak_prevalence": float(arr[peak_idx]),
        "time_of_peak": peak_idx,
    }
    pd.DataFrame([row]).to_csv(save_path, index=False)


def write_prevalence_outputs(
    prevalence: Sequence[float], *, prefix: str = "baseline", dpi: int = 300
) -> dict:
    """Convenience: write prevalence plot and summary under results/.

    Returns a dict of output paths.
    """
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    fig_path = FIGURES_DIR / f"prevalence_{prefix}.png"
    csv_path = RESULTS_DIR / f"summary_stats_{prefix}.csv"

    plot_and_save_prevalence(prevalence, fig_path, show=False, dpi=dpi)
    save_summary_stats(prevalence, csv_path)

    return {"figure": str(fig_path), "summary_csv": str(csv_path)}


__all__ = [
    "prevalence_from_states",
    "incidence_from_states",
    "plot_and_save_prevalence",
    "save_summary_stats",
    "write_prevalence_outputs",
    "write_run_metadata",
]


def _git_sha() -> Optional[str]:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(ROOT))
        return sha.decode().strip()
    except Exception:
        return None


def write_run_metadata(
    meta: Dict[str, Any],
    out_path: Optional[Path] = None,
    *,
    prefix: Optional[str] = None,
) -> Path:
    """Write a small JSON file with run metadata under results/.

    If out_path is None, writes to results/run_meta.json (or run_meta_{prefix}.json if prefix provided).
    Automatically adds timestamp, python version, platform, and git SHA (if available).
    Returns the path written.
    """
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if out_path is None:
        name = f"run_meta_{prefix}.json" if prefix else "run_meta.json"
        out_path = RESULTS_DIR / name

    enriched = dict(meta)
    enriched.update(
        {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
        }
    )
    sha = _git_sha()
    if sha:
        enriched["git_sha"] = sha

    _ensure_parent_dir(out_path)
    out_path.write_text(json.dumps(enriched, indent=2))
    return out_path
