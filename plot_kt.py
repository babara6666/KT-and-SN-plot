"""
plot_kt.py — Kitagawa-Takahashi (K-T) Diagram Plotter
AUT Master's Dissertation: Fatigue Analysis of LPBF Co-29Cr-6Mo Alloy

Reads kt_input/kt_data.csv and produces 10 PNG files in kt_output/:
  - 4 individual condition plots × 2 Y values  = 8 files
  - 1 combined plot × 2 Y values               = 2 files
  Total = 10 files

Original bilinear K-T boundary: horizontal fatigue-limit line meeting the
LEFM slope at the intersection point, using ΔK_th values from Anuar et al. (2021).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.ticker

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
Y_VALUES     = [0.5, 0.65]
R_RATIO      = 0.1          # stress ratio R
RUNOUT_LIMIT = 1_000_000    # cycles — specimens at this count are run-outs

DPI = 300
FIG_W_MM, FIG_H_MM = 180, 130

ROOT      = Path(__file__).parent
KT_INPUT  = ROOT / "kt_input"  / "kt_data.csv"
KT_OUTPUT = ROOT / "kt_output"

REQUIRED_COLS = {"specimen_id", "group", "stress_max_MPa", "sqrt_area", "Nf_cycles", "runout"}

# ΔK_th per condition [MPa·m^0.5] — Anuar et al. (2021)
DELTA_K_TH = {
    "180W BD//CD": 6.2,
    "180W BD⊥CD": 4.9,
    "320W BD//CD": 6.7,
    "320W BD⊥CD": 5.2,
}

# Short filename suffix per condition
GROUP_CODES = {
    "180W BD//CD": "180H",
    "180W BD⊥CD": "180V",
    "320W BD//CD": "320H",
    "320W BD⊥CD": "320V",
}

# Visual styles (marker + line colour) — shared with S-N plot
GROUP_STYLES = {
    "180W BD//CD": {"color": "blue",    "marker": "D"},
    "180W BD⊥CD": {"color": "red",     "marker": "^"},
    "320W BD//CD": {"color": "green",   "marker": "s"},
    "320W BD⊥CD": {"color": "magenta", "marker": "o"},
}

MARKER_SIZE  = 6
LINE_WIDTH   = 1.5
CMAP_NAME    = "viridis"
SQRT_AREA_PLOT = np.logspace(np.log10(10), np.log10(100000), 400)  # μm, for curve


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
    # Drop rows with non-positive sqrt_area or Nf_cycles
    bad = df[(df["sqrt_area"] <= 0) | (df["Nf_cycles"] <= 0)]
    if not bad.empty:
        print(f"Warning: dropping {len(bad)} rows with non-positive sqrt_area or Nf_cycles: "
              f"{list(bad['specimen_id'])}")
        df = df[(df["sqrt_area"] > 0) & (df["Nf_cycles"] > 0)].copy()
    return df


def delta_sigma_e(group_df: pd.DataFrame) -> float:
    """
    Fatigue limit Δσ_e for a condition.
    Taken from the run-out specimen's stress_max_MPa × (1 − R).
    """
    ro = group_df[group_df["runout"] == 1]
    if ro.empty:
        # Fallback: minimum stress in group
        s = group_df["stress_max_MPa"].min()
        print(f"Warning: no run-out specimen found — using min stress {s} MPa as Δσ_e reference.")
    else:
        s = ro["stress_max_MPa"].iloc[0]
    return s * (1 - R_RATIO)


def kt_boundary(dk_th: float, Y: float, ds_e: float):
    """
    Compute the original bilinear Kitagawa-Takahashi boundary.

    Using Murakami's √area model (√area in μm on x-axis):
        LEFM line:  Δσ = ΔK_th / (Y · √(π · √area_m))
                    where √area_m = √area_μm × 1e-6  [m]
        Horizontal: Δσ = Δσ_e

    Intersection:
        √area* = (1/π) · (ΔK_th / (Y · Δσ_e))² × 1e6  [μm]

    Returns (x_horiz, y_horiz, x_lefm, y_lefm) — two segments to plot.
    """
    x_min = SQRT_AREA_PLOT[0]
    x_max = SQRT_AREA_PLOT[-1]

    sqrt_area_star = (1.0 / np.pi) * (dk_th / (Y * ds_e)) ** 2 * 1e6  # μm

    # Horizontal segment: from x_min to intersection (clipped to plot range)
    x_star_clipped = min(max(sqrt_area_star, x_min), x_max)
    x_horiz = np.array([x_min, x_star_clipped])
    y_horiz = np.array([ds_e, ds_e])

    # LEFM segment: from intersection to x_max (clipped to plot range)
    if sqrt_area_star < x_max:
        x_start = max(sqrt_area_star, x_min)
        x_lefm = np.logspace(np.log10(x_start), np.log10(x_max), 200)
        y_lefm = dk_th / (Y * np.sqrt(np.pi * x_lefm * 1e-6))
    else:
        x_lefm = np.array([])
        y_lefm = np.array([])

    return x_horiz, y_horiz, x_lefm, y_lefm


def abbreviated_id(specimen_id: str) -> str:
    """Strip the numeric power prefix from specimen ID. e.g. '180H3' → 'H3'."""
    # Remove leading digits (180, 320)
    import re
    m = re.match(r"^\d+([A-Za-z].*)", specimen_id)
    return m.group(1) if m else specimen_id


def make_fig():
    return plt.subplots(
        figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
        dpi=DPI,
        facecolor="white",
    )


def style_axes(ax, title_suffix=""):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, 100000)
    ax.set_ylim(10, 1000)
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([10, 20, 50, 100, 200, 500, 1000])
    ax.set_xlabel(r"$\sqrt{\mathrm{area}}$ (μm)", fontsize=9)
    ax.set_ylabel(r"$\Delta\sigma$ (MPa)", fontsize=9)
    ax.grid(which="major", color="lightgrey", linestyle="-",  linewidth=0.6, zorder=0)
    ax.grid(which="minor", color="lightgrey", linestyle="--", linewidth=0.3, zorder=0)
    ax.minorticks_on()


def add_colorbar(fig, ax, norm, label="$N_f$ cycles"):
    sm = cm.ScalarMappable(cmap=CMAP_NAME, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(label, fontsize=8)
    # Show actual cycle values on colourbar ticks
    tick_vals = [1e3, 1e4, 1e5, 1e6]
    cbar.set_ticks([np.log10(v) for v in tick_vals])
    cbar.set_ticklabels([f"$10^{int(np.log10(v))}$" for v in tick_vals])


# ---------------------------------------------------------------------------
# Plot builders
# ---------------------------------------------------------------------------
def plot_individual(grp: str, gdf: pd.DataFrame, Y: float,
                    norm: mcolors.Normalize, cmap, y_tag: str):
    """Generate one individual condition K-T plot and save it."""
    dk_th = DELTA_K_TH.get(grp)
    if dk_th is None:
        print(f"Warning: no ΔK_th for group '{grp}' — skipping individual plot.")
        return

    code   = GROUP_CODES.get(grp, grp.replace(" ", "_"))
    color  = GROUP_STYLES.get(grp, {}).get("color", "black")
    marker = GROUP_STYLES.get(grp, {}).get("marker", "o")

    ds_e   = delta_sigma_e(gdf)
    x_horiz, y_horiz, x_lefm, y_lefm = kt_boundary(dk_th, Y, ds_e)

    fig, ax = make_fig()

    # K-T bilinear boundary
    ax.plot(x_horiz, y_horiz, color=color, lw=LINE_WIDTH,
            linestyle="-", label=f"K-T boundary ({grp})", zorder=3)
    if x_lefm.size:
        ax.plot(x_lefm, y_lefm, color=color, lw=LINE_WIDTH,
                linestyle="-", zorder=3)

    # Data points
    for _, row in gdf.iterrows():
        c_val = cmap(norm(np.log10(row["Nf_cycles"])))
        ds    = row["stress_max_MPa"] * (1 - R_RATIO)
        is_ro = row["runout"] == 1

        ax.scatter(
            row["sqrt_area"], ds,
            color=c_val,
            marker=marker,
            s=MARKER_SIZE ** 2,
            edgecolors=color if is_ro else "none",
            linewidths=0.8,
            zorder=4,
        )
        # Abbreviated label
        ax.annotate(
            abbreviated_id(str(row["specimen_id"])),
            xy=(row["sqrt_area"], ds),
            xytext=(3, 3), textcoords="offset points",
            fontsize=6, color="dimgrey",
        )
        # Run-out arrow
        if is_ro:
            ax.annotate(
                "",
                xy=(row["sqrt_area"] * 1.15, ds),
                xytext=(row["sqrt_area"], ds),
                arrowprops=dict(arrowstyle="->", color=color, lw=0.9),
            )

    style_axes(ax)
    ax.text(0.04, 0.06, f"$Y$ = {Y:.2f}",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="grey", alpha=0.8))
    ax.legend(fontsize=7, loc="upper right")
    add_colorbar(fig, ax, norm)

    fig.tight_layout()
    out = KT_OUTPUT / f"KT_{y_tag}_{code}.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_combined(df: pd.DataFrame, groups: list, Y: float,
                  norm: mcolors.Normalize, cmap, y_tag: str):
    """Generate the combined (all-conditions) K-T plot and save it."""
    fig, ax = make_fig()

    for grp in groups:
        gdf   = df[df["group"] == grp]
        dk_th = DELTA_K_TH.get(grp)
        if dk_th is None:
            continue

        color  = GROUP_STYLES.get(grp, {}).get("color", "black")
        marker = GROUP_STYLES.get(grp, {}).get("marker", "o")
        ds_e   = delta_sigma_e(gdf)
        x_horiz, y_horiz, x_lefm, y_lefm = kt_boundary(dk_th, Y, ds_e)

        # K-T bilinear boundary
        ax.plot(x_horiz, y_horiz, color=color, lw=LINE_WIDTH,
                linestyle="-", label=grp, zorder=3)
        if x_lefm.size:
            ax.plot(x_lefm, y_lefm, color=color, lw=LINE_WIDTH,
                    linestyle="-", zorder=3)

        # Data points
        for _, row in gdf.iterrows():
            c_val = cmap(norm(np.log10(row["Nf_cycles"])))
            ds    = row["stress_max_MPa"] * (1 - R_RATIO)
            is_ro = row["runout"] == 1

            ax.scatter(
                row["sqrt_area"], ds,
                color=c_val,
                marker=marker,
                s=MARKER_SIZE ** 2,
                edgecolors=color if is_ro else "none",
                linewidths=0.8,
                zorder=4,
            )
            ax.annotate(
                abbreviated_id(str(row["specimen_id"])),
                xy=(row["sqrt_area"], ds),
                xytext=(3, 3), textcoords="offset points",
                fontsize=5, color="dimgrey",
            )
            if is_ro:
                ax.annotate(
                    "",
                    xy=(row["sqrt_area"] * 1.15, ds),
                    xytext=(row["sqrt_area"], ds),
                    arrowprops=dict(arrowstyle="->", color=color, lw=0.9),
                )

    style_axes(ax)
    ax.text(0.04, 0.06, f"$Y$ = {Y:.2f}",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="grey", alpha=0.8))
    ax.legend(fontsize=6, loc="upper right", framealpha=0.9)
    add_colorbar(fig, ax, norm)

    fig.tight_layout()
    out = KT_OUTPUT / f"KT_{y_tag}_combined.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    KT_OUTPUT.mkdir(parents=True, exist_ok=True)

    df = load_data(KT_INPUT)

    # Ordered group list
    groups_in_data = list(df["group"].unique())
    ordered_groups = [g for g in GROUP_STYLES if g in groups_in_data]
    unknown_groups = [g for g in groups_in_data if g not in GROUP_STYLES]
    ordered_groups += unknown_groups

    # Shared colourmap normalised to log10(Nf_cycles) range across all data
    log_nf_all = np.log10(df["Nf_cycles"].values.astype(float))
    norm = mcolors.Normalize(vmin=log_nf_all.min(), vmax=log_nf_all.max())
    cmap = matplotlib.colormaps[CMAP_NAME]

    for Y in Y_VALUES:
        y_tag = f"Y{int(Y * 100):03d}"   # 0.5 → 'Y050', 0.65 → 'Y065'

        # Individual plots
        for grp in ordered_groups:
            gdf = df[df["group"] == grp].copy()
            plot_individual(grp, gdf, Y, norm, cmap, y_tag)

        # Combined plot
        plot_combined(df, ordered_groups, Y, norm, cmap, y_tag)


if __name__ == "__main__":
    main()
