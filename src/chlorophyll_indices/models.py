"""Linear index -> chlorophyll models with honest out-of-sample error."""

from __future__ import annotations

import numpy as np
import pandas as pd


def fit(x: np.ndarray, y: np.ndarray) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        raise ValueError("Need at least three valid samples")
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    r2 = 1 - np.sum(residual**2) / np.sum((y - y.mean()) ** 2)
    return {"n": int(x.size), "slope": float(slope), "intercept": float(intercept), "r2": float(r2),
            "rmse": float(np.sqrt(np.mean(residual**2)))}


def loo_rmse(x: np.ndarray, y: np.ndarray) -> float:
    """Leave-one-out RMSE, computed exactly from the hat matrix of a simple linear fit."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    leverage = 1 / x.size + (x - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2)
    return float(np.sqrt(np.mean((residual / (1 - leverage)) ** 2)))


def leave_group_out_rmse(x: np.ndarray, y: np.ndarray, groups: np.ndarray) -> float:
    """RMSE when each group (species) is predicted by a model fitted on all the others."""
    x, y, groups = np.asarray(x, float), np.asarray(y, float), np.asarray(groups)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, groups = x[ok], y[ok], groups[ok]
    errors = []
    for g in np.unique(groups):
        test = groups == g
        slope, intercept = np.polyfit(x[~test], y[~test], 1)
        errors.append(y[test] - (slope * x[test] + intercept))
    return float(np.sqrt(np.mean(np.concatenate(errors) ** 2)))


def evaluate(df: pd.DataFrame, indices: list[str], target: str, group: str | None = None) -> pd.DataFrame:
    rows = []
    for name in indices:
        stats = fit(df[name], df[target])
        row = {"index": name, "target": target, **stats, "loo_rmse": loo_rmse(df[name], df[target])}
        if group:
            row["leave_species_out_rmse"] = leave_group_out_rmse(df[name], df[target], df[group])
        rows.append(row)
    return pd.DataFrame(rows).sort_values("r2", ascending=False).reset_index(drop=True)


def per_group_r2(df: pd.DataFrame, indices: list[str], target: str, group: str) -> pd.DataFrame:
    table = {name: {g: fit(sub[name], sub[target])["r2"] for g, sub in df.groupby(group)} for name in indices}
    return pd.DataFrame(table)
