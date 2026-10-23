"""Sentinel-2 band simulation from leaf spectra, and the five vegetation indices.

Index definitions follow the course/Frampton et al. (2013) convention, with B7 (783 nm) as NIR.
"""

from __future__ import annotations

import numpy as np

# Sentinel-2A centre wavelength and nominal bandwidth (nm)
S2_BANDS = {"B3": (560, 36), "B4": (665, 31), "B5": (705, 15), "B6": (740, 15), "B7": (783, 20)}

INDICES = {
    "NDVI": "(B7 − B4) / (B7 + B4)",
    "GNDVI": "(B7 − B3) / (B7 + B3)",
    "MTCI": "(B6 − B5) / (B5 − B4)",
    "IRECI": "(B7 − B4) / (B5 / B6)",
    "PSSRa": "B7 / B4",
}


def simulate_bands(wavelength_nm: np.ndarray, reflectance: np.ndarray, method: str = "point") -> dict[str, np.ndarray]:
    """Band reflectance from spectra (rows = wavelengths, columns = leaves).

    ``point`` samples the centre wavelength (reproduces the original analysis);
    ``window`` averages over the nominal bandwidth (closer to a real sensor response).
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    reflectance = np.asarray(reflectance, dtype=float)
    out = {}
    for band, (centre, width) in S2_BANDS.items():
        if method == "point":
            rows = np.flatnonzero(np.isclose(wavelength_nm, centre))
            if rows.size == 0:
                raise ValueError(f"Spectra do not include {centre} nm for {band}")
            out[band] = reflectance[rows[0]]
        elif method == "window":
            mask = (wavelength_nm >= centre - width / 2) & (wavelength_nm <= centre + width / 2)
            if not mask.any():
                raise ValueError(f"Spectra do not cover {band} ({centre} ± {width / 2} nm)")
            out[band] = reflectance[mask].mean(axis=0)
        else:
            raise ValueError("method must be 'point' or 'window'")
    return out


def compute_indices(b: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore"):
        return {
            "NDVI": (b["B7"] - b["B4"]) / (b["B7"] + b["B4"]),
            "GNDVI": (b["B7"] - b["B3"]) / (b["B7"] + b["B3"]),
            "MTCI": (b["B6"] - b["B5"]) / (b["B5"] - b["B4"]),
            "IRECI": (b["B7"] - b["B4"]) / (b["B5"] / b["B6"]),
            "PSSRa": b["B7"] / b["B4"],
        }
