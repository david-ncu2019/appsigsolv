"""Visualization and Report Generation."""
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec

def save_plot(series: pd.Series, components: dict, best_model: dict, comp: str, unit: str, out_root: Path, stem: str) -> Path:
    scale = 1000.0 if unit == 'mm' else 1.0

    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 20,
        'axes.labelsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'figure.titlesize': 24
    })

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1], wspace=0.2)

    idx = series.index
    periodic_cols = [c for c in components if c.endswith("yr")]
    has_jump = f"{comp}_jump" in components
    has_exp = f"{comp}_exp" in components
    has_log = f"{comp}_log" in components
    has_exp_trend = f"{comp}_exp_trend" in components

    gs_left = gs[0, 0].subgridspec(3, 1, hspace=0.3)
    ax_main = fig.add_subplot(gs_left[0:2, 0])
    ax_trend = fig.add_subplot(gs_left[2, 0], sharex=ax_main)

    ax_main.plot(idx, series.values * scale, "o", ms=5, alpha=0.75, color="#95a5a6", label="Observed")
    ax_main.plot(idx, components[f"{comp}_model"].values * scale, "-", lw=2, color="#e74c3c", label="Model")
    stats = best_model.get("_omt_stats", {"sigma_mm": 1.0, "p_value": 1.0, "n_param": 1})
    ax_main.set_title(f"{stem} {comp}: Observed vs. Model", fontweight='bold')
    ax_main.set_ylabel(f"Displacement ({unit})")
    ax_main.legend(loc="lower left", frameon=True, framealpha=0.9)
    ax_main.grid(True, ls=":", alpha=0.6)
    
    ax_trend.plot(idx, components[f"{comp}_trend"].values * scale, "-", lw=2.5, color="#2c3e50", label="Poly trend")
    if has_exp_trend:
        ax_trend.plot(idx, components[f"{comp}_exp_trend"].values * scale, "--", lw=2, color="#8e44ad", label="Exp trend")
        ax_trend.legend(loc="lower left", fontsize=12)
    ax_trend.set_title("Long-term Trend", fontweight='bold')
    ax_trend.set_ylabel(f"Trend ({unit})")
    ax_trend.grid(True, ls=":", alpha=0.6)

    right_types = []
    if periodic_cols: right_types.append("seasonal")
    if has_jump or has_exp or has_log: right_types.append("jump")
    right_types.append("residuals")
    
    n_right = len(right_types)
    gs_right = gs[0, 1].subgridspec(n_right, 1, hspace=0.3)
    
    right_axes = []
    for i, ptype in enumerate(right_types):
        ax = fig.add_subplot(gs_right[i], sharex=ax_main)
        ax.grid(True, ls=":", alpha=0.6)
        right_axes.append(ax)
        
        if ptype == "seasonal":
            for col in periodic_cols:
                label = col.replace(f"{comp}_", "")
                ax.plot(idx, components[col].values * scale, lw=2, label=label)
            ax.axhline(0, color="black", lw=1, ls="--")
            ax.set_title("Seasonal Components", fontweight='bold')
            ax.set_ylabel(f"Seasonal ({unit})")
            ax.legend(loc="upper right", fontsize=12, ncol=2)
            
        elif ptype == "jump":
            if has_jump:
                ax.plot(idx, components[f"{comp}_jump"].values * scale, "-", lw=2, color="#9b59b6", label="Jump")
            if has_exp:
                ax.plot(idx, components[f"{comp}_exp"].values * scale, "--", lw=2, color="#e67e22", label="Exp Relax")
            if has_log:
                ax.plot(idx, components[f"{comp}_log"].values * scale, ":", lw=2, color="#d35400", label="Log Relax")
            ax.axhline(0, color="black", lw=1, ls="--")
            ax.set_title("Jumps & Relaxation", fontweight='bold')
            ax.set_ylabel(f"Jump ({unit})")
            ax.legend(loc="upper left", fontsize=12)
            
        elif ptype == "residuals":
            ax.plot(idx, components[f"{comp}_noise"].reindex(series.index).values * scale, "o", ms=2, alpha=0.4, color="#3498db")
            ax.axhline(0, color="black", lw=1, ls="--")
            ax.set_title("Residual Noise", fontweight='bold')
            ax.set_ylabel(f"Residual ({unit})")

    for ax in [ax_main, ax_trend] + right_axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(2))
        for label in ax.get_xticklabels():
            label.set_rotation(0)
            label.set_horizontalalignment('center')

    plt.suptitle(f"Timeseries Decomposition: {stem} ({comp})\n"
                 f"σ={stats['sigma_mm']:.1f}{unit} | p={stats['p_value']:.4f} | n_param={stats['n_param']} | Poly-deg={best_model.get('polynomial')}",
                 y=0.98, fontweight='bold')

    png_path = out_root / f"{stem}_decomposed_{comp}.png"
    plt.savefig(png_path, dpi=100, bbox_inches="tight")
    plt.close()
    print(f"  [output] Publication-quality dynamic plot saved: {png_path}")
    return png_path

def save_report(series: pd.Series, components: dict, best_model: dict, scan_table: list, comp: str, unit: str, out_root: Path, stem: str) -> Path:
    stats = best_model.get("_omt_stats", {"sigma_mm": 1.0, "p_value": 1.0, "n_param": 1, "omt": 1.0, "K_norm": 1.0, "iterations": 1})
    model = {k: v for k, v in best_model.items() if not k.startswith("_")}
    scale = 1000.0 if unit == 'mm' else 1.0

    lines = [
        f"# Timeseries Decomposition Report: {stem} — {comp}",
        "",
        "## Accepted Model",
        "",
        f"| Parameter | Value |",
        f"|---|---|",
        f"| Component | {comp} |",
        f"| Sigma ({unit}) | {stats['sigma_mm']:.1f} |",
        f"| Polynomial degree | {model.get('polynomial', 0)} |",
        f"| Seasonal periods (yr) | {model.get('periodic', [])} |",
        f"| Jump dates | {model.get('stepDate', [])} |",
        f"| Polyline breaks | {model.get('polyline', [])} |",
        f"| Exp trend (b/day) | {model.get('exp_trend', None)} |",
        f"| Exp relaxation | {model.get('exp', {})} |",
        f"| Log relaxation | {model.get('log', {})} |",
        f"| n_params | {stats['n_param']} |",
        f"| Normalized OMT | {stats['omt']:.4f} |",
        f"| K/r threshold | {stats['K_norm']:.4f} |",
        f"| p-value | {stats['p_value']:.6f} |",
        f"| DIA iterations | {stats['iterations']} |",
        "",
        "## Variance Explained per Component (%)",
        "",
        f"| Component | Std ({unit}) | Variance Explained (%) |",
        "|---|---|---|",
    ]

    series_clean = series.dropna()
    orig_var = np.var(series_clean.values)
    for col, s in components.items():
        if col in (f"{comp}_model", f"{comp}_noise", f"{comp}_wtest"):
            continue
        s_clean = s.reindex(series_clean.index)
        pct = (np.var(s_clean.values) / orig_var * 100) if orig_var > 0 else 0.0
        lines.append(f"| {col} | {np.std(s_clean.values)*scale:.2f} | {pct:.2f} |")

    noise_clean = components[f"{comp}_noise"].reindex(series_clean.index)
    noise_std = np.std(noise_clean.values) * scale
    noise_pct = np.var(noise_clean.values) / orig_var * 100 if orig_var > 0 else 0.0
    lines.append(f"| {comp}_noise | {noise_std:.2f} | {noise_pct:.2f} |")

    lines += [
        "",
        "## Sigma Scan Summary",
        "",
        f"| Sigma ({unit}) | Accepted | p-value | n_params | n_periods | n_polylines |",
        "|---|---|---|---|---|---|",
    ]
    for row in scan_table:
        p = f"{row['p_value']:.4f}" if row["p_value"] is not None else "—"
        n = str(row["n_param"]) if row["n_param"] is not None else "—"
        np_ = str(row["n_periods"]) if row["n_periods"] is not None else "—"
        nl = str(row["n_polylines"]) if row["n_polylines"] is not None else "—"
        acc = "✓" if row["accepted"] else "✗"
        lines.append(f"| {row['sigma_mm']:.1f} | {acc} | {p} | {n} | {np_} | {nl} |")

    flagged = components[f"{comp}_wtest"][np.abs(components[f"{comp}_wtest"]) > 3.29]
    if not flagged.empty:
        lines += [
            "",
            "## Anomalous Observations (w-test > 3.29)",
            "",
            "| Date | w-stat | Residual (mm) |",
            "|---|---|---|",
        ]
        for dt, w in flagged.items():
            res_mm = components[f"{comp}_noise"].loc[dt] * 1000.0
            lines.append(f"| {dt.strftime('%Y-%m-%d')} | {w:.2f} | {res_mm:.2f} |")

    lines.append("")
    report_path = out_root / f"{stem}_report_{comp}.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  [output] Report saved: {report_path}")
    return report_path
