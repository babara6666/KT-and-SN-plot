"""
plot_sn.py — S-N Curve Plotter
AUT Master's Dissertation: Fatigue Analysis of LPBF Co-29Cr-6Mo Alloy

Reads sn_input/fatigue_data.csv and produces a single log-log S-N plot
(sn_output/SN_curve.png) with Basquin regression per group.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RUNOUT_LIMIT = 1_000_000          # cycles — specimens at this value are run-outs
DPI = 300
FIG_W_MM, FIG_H_MM = 180, 130    # figure dimensions in mm

ROOT = Path(__file__).parent
SN_INPUT  = ROOT / "sn_input"  / "fatigue_data.csv"
SN_OUTPUT = ROOT / "sn_output"
BASQUIN_CSV = SN_OUTPUT / "basquin_coefficients.csv"

REQUIRED_COLS = {"specimen_id", "group", "stress_max_MPa", "cycles", "runout"}

# Group visual styles — order determines legend order
GROUP_STYLES = {
    "180W BD//CD": {"color": "blue",    "marker": "D"},
    "180W BD⊥CD": {"color": "red",     "marker": "^"},
    "320W BD//CD": {"color": "green",   "marker": "s"},
    "320W BD⊥CD": {"color": "magenta", "marker": "o"},
}

MARKER_SIZE   = 6
LINE_WIDTH    = 1.2
ARROW_OFFSET  = 0.15   # log10 units rightward for run-out arrow


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_data(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")
    df["runout"] = df["runout"].astype(int)
    return df


def basquin_fit(log_N: np.ndarray, log_s: np.ndarray):
    """Return (b, log10_sigma_f_prime) from polyfit(log_N, log_s, 1)."""
    coeffs = np.polyfit(log_N, log_s, 1)   # [slope, intercept] = [b, log10_sf]
    return coeffs[0], coeffs[1]


def sigma_from_basquin(b: float, log_sf: float, N: float) -> float:
    """Predict σ_max at given N using Basquin relation."""
    return 10 ** (log_sf + b * np.log10(N))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    SN_OUTPUT.mkdir(parents=True, exist_ok=True)

    df = load_data(SN_INPUT)

    # Detect groups in the order they appear (fall back to GROUP_STYLES order
    # if a known group is present, preserving the visual assignment).
    groups_in_data = list(df["group"].unique())
    ordered_groups = [g for g in GROUP_STYLES if g in groups_in_data]
    unknown_groups = [g for g in groups_in_data if g not in GROUP_STYLES]
    if unknown_groups:
        print(f"Warning: unknown groups will use default style: {unknown_groups}")
    ordered_groups += unknown_groups

    # Default style for unknown groups
    _default_colors  = ["tab:orange", "tab:brown", "tab:pink", "tab:cyan"]
    _default_markers = ["P", "X", "h", "v"]

    fig, ax = plt.subplots(
        figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
        dpi=DPI,
        facecolor="white",
    )

    basquin_records = []

    for idx, grp in enumerate(ordered_groups):
        gdf = df[df["group"] == grp].copy()

        # Resolve style
        if grp in GROUP_STYLES:
            color  = GROUP_STYLES[grp]["color"]
            marker = GROUP_STYLES[grp]["marker"]
        else:
            color  = _default_colors[idx % len(_default_colors)]
            marker = _default_markers[idx % len(_default_markers)]

        failed  = gdf[gdf["runout"] == 0]
        runouts = gdf[gdf["runout"] == 1]

        # ---- Scatter: failed specimens ----
        ax.scatter(
            failed["cycles"],
            failed["stress_max_MPa"],
            color=color,
            marker=marker,
            s=MARKER_SIZE ** 2,
            zorder=3,
            label=grp,
        )

        # ---- Run-out specimens ----
        for _, row in runouts.iterrows():
            ax.scatter(
                row["cycles"],
                row["stress_max_MPa"],
                color=color,
                marker=marker,
                s=MARKER_SIZE ** 2,
                zorder=3,
            )
            # Rightward arrow annotation
            ax.annotate(
                "",
                xy=(row["cycles"] * 10 ** ARROW_OFFSET, row["stress_max_MPa"]),
                xytext=(row["cycles"], row["stress_max_MPa"]),
                arrowprops=dict(arrowstyle="->", color=color, lw=LINE_WIDTH),
            )

        # ---- Basquin regression ----
        if len(failed) < 2:
            print(f"Warning: group '{grp}' has fewer than 2 failed specimens — skipping regression.")
            basquin_records.append({"group": grp, "sigma_f_prime": None, "b": None})
            continue

        log_N = np.log10(failed["cycles"].values.astype(float))
        log_s = np.log10(failed["stress_max_MPa"].values.astype(float))
        b, log_sf = basquin_fit(log_N, log_s)
        print(f"  {grp}: b = {b:.4f}, σ'f = {10**log_sf:.1f} MPa")

        # Fit line spans each group's own data range only
        N_min = failed["cycles"].min()
        N_max = RUNOUT_LIMIT
        N_fit = np.logspace(np.log10(N_min), np.log10(N_max), 300)
        s_fit = 10 ** (log_sf + b * np.log10(N_fit))
        ax.plot(N_fit, s_fit, color=color, lw=LINE_WIDTH, linestyle="-", zorder=2)

        basquin_records.append({
            "group": grp,
            "sigma_f_prime": 10 ** log_sf,
            "b": b,
        })

    # ---- Axes ----
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e3, 1e7)
    ax.set_ylim(200, 1100)

    ax.set_xlabel("Number of Cycles to Failure, $N$", fontsize=9)
    ax.set_ylabel(r"Maximum Stress, $\sigma_\mathrm{max}$ (MPa)", fontsize=9)

    # Y-axis: show actual MPa values rather than powers of 10
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_yticks([200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100])

    # Grid
    ax.grid(which="major", color="lightgrey", linestyle="-",  linewidth=0.6, zorder=0)
    ax.grid(which="minor", color="lightgrey", linestyle="--", linewidth=0.3, zorder=0)
    ax.minorticks_on()

    # ---- Annotations ----
    ax.text(
        0.97, 0.97,
        "$f$ = 30 Hz,  $R$ = 0.1",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.8),
    )

    # ---- Legend ----
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9, bbox_to_anchor=(0.99, 0.88))

    fig.tight_layout()
    out_path = SN_OUTPUT / "SN_curve.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")

    # ---- Export Basquin coefficients ----
    bdf = pd.DataFrame(basquin_records)
    bdf.to_csv(BASQUIN_CSV, index=False)
    print(f"Saved: {BASQUIN_CSV}")


if __name__ == "__main__":
    main()
