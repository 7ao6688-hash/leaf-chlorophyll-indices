from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .data import read_lab, read_satellite
from .indices import INDICES
from .models import evaluate, per_group_r2


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="chlorophyll-indices",
                                description="Evaluate vegetation indices against measured leaf and canopy chlorophyll")
    p.add_argument("--spectra", type=Path, required=True, help="Workbook of leaf spectra, one sheet per species")
    p.add_argument("--chlorophyll", type=Path, required=True, help="Workbook of measured leaf chlorophyll")
    p.add_argument("--satellite", type=Path, help="Optional Sentinel-2 field-plot table with LCC/LAI/CCC")
    p.add_argument("--band-method", choices=("point", "window"), default="point")
    p.add_argument("--output-dir", type=Path, default=Path("results"))
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args(argv)

    names = list(INDICES)
    try:
        lab = read_lab(args.spectra, args.chlorophyll, args.band_method)
        sat = read_satellite(args.satellite) if args.satellite else None
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        print(f"Error reading inputs: {exc}", file=sys.stderr)
        return 1

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    pooled = evaluate(lab, names, "lcc_ug_cm2", group="species")
    species = per_group_r2(lab, names, "lcc_ug_cm2", "species")
    pooled.round(4).to_csv(out / "lab_pooled_models.csv", index=False)
    species.round(4).to_csv(out / "lab_r2_by_species.csv", index_label="species")
    summary = {
        "band_method": args.band_method,
        "lab": {"leaves": int(len(lab)), "species": int(lab.species.nunique()),
                "best_pooled": pooled.iloc[0]["index"],
                "best_within_species": species.mean().idxmax(),
                "mean_within_species_r2": species.mean().round(3).to_dict(),
                "pooled_r2": pooled.set_index("index").r2.round(3).to_dict()},
    }
    if sat is not None:
        sat_lcc = evaluate(sat, names, "lcc_g_m2")
        sat_ccc = evaluate(sat, names, "ccc_g_m2")
        sat_lcc.round(4).to_csv(out / "satellite_lcc_models.csv", index=False)
        sat_ccc.round(4).to_csv(out / "satellite_ccc_models.csv", index=False)
        summary["satellite"] = {"plots": int(len(sat)), "dates": sorted({str(d) for d in sat.date}),
                                "lcc_r2": sat_lcc.set_index("index").r2.round(3).to_dict(),
                                "ccc_r2": sat_ccc.set_index("index").r2.round(3).to_dict()}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    if not args.no_figures:
        from . import report

        report.lab_scatter(out, lab, pooled)
        report.species_heatmap(out, species)
        if sat is not None:
            report.r2_comparison(out, species.mean(), pooled, sat_lcc, sat_ccc)
            report.satellite_scatter(out, sat, sat_ccc)

    print(f"Lab: {len(lab)} leaves, {lab.species.nunique()} species")
    print(pooled[["index", "r2", "rmse", "loo_rmse", "leave_species_out_rmse"]].round(3).to_string(index=False))
    if sat is not None:
        print(f"\nSentinel-2 plots: {len(sat)}  (canopy chlorophyll)")
        print(sat_ccc[["index", "r2", "rmse", "loo_rmse"]].round(3).to_string(index=False))
    print(f"\nResults: {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
