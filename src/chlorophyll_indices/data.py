"""Readers for the course workbook layouts, returning tidy tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .indices import compute_indices, simulate_bands


def read_lab(spectra_xlsx: Path, chlorophyll_xlsx: Path, method: str = "point") -> pd.DataFrame:
    """One row per leaf: species, chlorophyll (µg/cm²) and the five indices from its spectrum.

    Spectra sheets: column B = wavelength, then one column per leaf ("Leaf_1"...), header in row 1.
    Chlorophyll sheets: "Leaf No." and "Chlorophyll Concentration (microgram/cm2)".
    Sheets are matched by species number, which tolerates names such as "Group__Species_5".
    """
    spectra = pd.read_excel(spectra_xlsx, sheet_name=None, header=None)
    chlorophyll = pd.read_excel(chlorophyll_xlsx, sheet_name=None)

    def number(name: str) -> int:
        return int("".join(ch for ch in name.split("Species")[-1] if ch.isdigit()))

    chl_by_species = {number(k): v for k, v in chlorophyll.items()}
    frames = []
    for sheet, raw in spectra.items():
        species = number(sheet)
        header = raw.iloc[0].tolist()
        leaf_cols = [i for i, h in enumerate(header) if isinstance(h, str) and h.startswith("Leaf_")]
        wavelength = pd.to_numeric(raw.iloc[1:, 1], errors="coerce")
        rows = wavelength.notna().to_numpy()
        values = raw.iloc[1:, leaf_cols].to_numpy()[rows].astype(float)
        idx = compute_indices(simulate_bands(wavelength[rows].to_numpy(), values, method))
        conc = chl_by_species[species]
        chl_col = next(c for c in conc.columns if str(c).startswith("Chlorophyll"))
        n = min(len(leaf_cols), len(conc))
        frame = pd.DataFrame({name: v[:n] for name, v in idx.items()})
        frame.insert(0, "lcc_ug_cm2", conc[chl_col].to_numpy(float)[:n])
        frame.insert(0, "leaf", np.arange(1, n + 1))
        frame.insert(0, "species", f"Species {species}")
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def read_satellite(path: Path) -> pd.DataFrame:
    """Sentinel-2 L2A reflectance at field plots with measured LCC, LAI and CCC (g/m²).

    Indices are recomputed from the band columns rather than trusted from the sheet.
    """
    raw = pd.read_excel(path, header=1)
    bands = {b: pd.to_numeric(raw[b], errors="coerce").to_numpy() for b in ("B3", "B4", "B5", "B6", "B7")}
    out = pd.DataFrame({
        "plot": raw["Name"],
        "latitude": pd.to_numeric(raw["Latitude"], errors="coerce"),
        "longitude": pd.to_numeric(raw["Longitude"], errors="coerce"),
        "date": pd.to_datetime(raw["Date(yyyy-MM-dd)"], errors="coerce").dt.date,
        "lcc_g_m2": pd.to_numeric(raw["LCC (g/m^2)"], errors="coerce"),
        "lai": pd.to_numeric(raw["LAI"], errors="coerce"),
        "ccc_g_m2": pd.to_numeric(raw["CCC (g/m^2)"], errors="coerce"),
    })
    for name, values in compute_indices(bands).items():
        out[name] = values
    return out.dropna(subset=["lcc_g_m2", "ccc_g_m2"]).reset_index(drop=True)
