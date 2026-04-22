"""
plot_kt.py — Kitagawa-Takahashi Diagram Plotter
Reads 'D:/DownloadD/WAAM/LPBF/SN force.xlsx' and produces 6 PNG files in kt_output/:
  KT_Y050_combined.png  KT_Y065_combined.png  — all four groups
  KT_Y050_180W.png      KT_Y065_180W.png      — 180W groups overlaid
  KT_Y050_320W.png      KT_Y065_320W.png      — 320W groups overlaid
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RUNOUT_LIMIT = 10_000_000
R_RATIO      = 0.1
Y_VALUES     = [0.5, 0.65]
DPI          = 300
FIG_W_MM, FIG_H_MM = 180, 130
EXCEL_PATH   = Path("D:/DownloadD/WAAM/LPBF/SN force.xlsx")
OUT_DIR      = Path(__file__).parent / "kt_output"

# Fatigue limits Δσ_e [MPa] — user-specified
FATIGUE_LIMITS = {
    "180W BD//LD": 250,
    "180W BD⊥LD": 400,
    "320W BD//LD": 350,
    "320W BD⊥LD": 500,
}

# Threshold SIF range ΔK_th [MPa·√m] — Anuar et al. (2021)
DELTA_K_TH = {
    "180W BD//LD": 6.2,
    "180W BD⊥LD": 4.9,
    "320W BD//LD": 6.7,
    "320W BD⊥LD": 5.2,
}

GROUP_STYLES = {
    "180W BD//LD": {"color": "blue",    "marker": "D"},
    "180W BD⊥LD": {"color": "red",     "marker": "^"},
    "320W BD//LD": {"color": "green",   "marker": "s"},
    "320W BD⊥LD": {"color": "magenta", "marker": "o"},
}

WATTAGE_GROUPS = {
    "180W": ["180W BD//LD", "180W BD⊥LD"],
    "320W": ["320W BD//LD", "320W BD⊥LD"],
}

MARKER_SIZE  = 6
LINE_WIDTH   = 1.5
CMAP_NAME    = "RdYlBu"   # warm (red, low Nf) → cool (blue, high Nf)
CBAR_LOG_MIN = 4.0         # log10(1e4)
CBAR_LOG_MAX = 7.0         # log10(1e7) — matches RUNOUT_LIMIT
ARROW_SCALE  = 1.15
X_CURVE      = np.logspace(np.log10(50), np.log10(2000), 400)  # µm


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fix_group(g):
    if "//" in g:
        return g
    prefix = g.split("W")[0] + "W"
    return f"{prefix} BD⊥LD"


def shorten_id(specimen_id):
    m = re.match(r"^\d+([A-Za-z].*)", str(specimen_id))
    return m.group(1) if m else str(specimen_id)


def kt_boundary(dk_th, Y, ds_e):
    """Bilinear K-T boundary. Returns horizontal and LEFM segments."""
    x_min, x_max   = X_CURVE[0], X_CURVE[-1]
    sqrt_area_star  = (1.0 / np.pi) * (dk_th / (Y * ds_e)) ** 2 * 1e6  # µm
    x_star          = min(max(sqrt_area_star, x_min), x_max)

    x_horiz = np.array([x_min, x_star])
    y_horiz = np.array([ds_e,  ds_e])

    if sqrt_area_star < x_max:
        x_lefm = np.logspace(np.log10(max(sqrt_area_star, x_min)),
                              np.log10(x_max), 200)
        y_lefm = dk_th / (Y * np.sqrt(np.pi * x_lefm * 1e-6))
    else:
        x_lefm, y_lefm = np.array([]), np.array([])

    return x_horiz, y_horiz, x_lefm, y_lefm


def load_data():
    df = pd.read_excel(EXCEL_PATH)
    df = df[df["stress_max_MPa"] > 0].dropna(subset=["cycles"]).copy()
    df["group"]       = df["group"].apply(fix_group)
    df["runout"]      = (df["cycles"] >= RUNOUT_LIMIT).astype(int)
    df["sqrt_area"]   = pd.to_numeric(df["sqrt area"], errors="coerce")
    df["delta_sigma"] = df["stress_max_MPa"] * (1 - R_RATIO)
    df = df[df["sqrt_area"] > 0].copy()
    return df


def style_ax(ax, Y):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(50, 2000)
    ax.set_ylim(100, 1000)
    ax.set_xlabel(r"$\sqrt{\mathrm{area}}$ (μm)", fontsize=9)
    ax.set_ylabel(r"$\Delta\sigma$ (MPa)", fontsize=9)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    ax.grid(which="major", color="lightgrey", linestyle="-",  linewidth=0.6, zorder=0)
    ax.grid(which="minor", color="lightgrey", linestyle="--", linewidth=0.3, zorder=0)
    ax.minorticks_on()
    ax.text(0.04, 0.06, f"$Y$ = {Y:.2f},  $R$ = {R_RATIO}",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="grey", alpha=0.8))


def add_colorbar(fig, ax, norm):
    sm = cm.ScalarMappable(cmap=CMAP_NAME, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("$N_f$ (cycles)", fontsize=8)
    tick_vals = [1e4, 1e5, 1e6, 1e7]
    cbar.set_ticks([np.log10(v) for v in tick_vals])
    cbar.set_ticklabels(["$10^4$", "$10^5$", "$10^6$", "$10^7$"])


def draw_groups(ax, fig, df, groups, Y, norm, cmap, show_wattage=True):
    for grp in groups:
        if grp not in df["group"].unique():
            continue
        gdf    = df[df["group"] == grp]
        color  = GROUP_STYLES[grp]["color"]
        marker = GROUP_STYLES[grp]["marker"]
        ds_e   = FATIGUE_LIMITS[grp]
        dk_th  = DELTA_K_TH[grp]
        # Strip "180W " / "320W " prefix for wattage-only plots, append H/V hint
        if show_wattage:
            legend_label = grp
        else:
            short = grp.split(" ", 1)[1]  # e.g. "BD//LD" or "BD⊥LD"
            hint  = "(H)" if "//" in short else "(V)"
            legend_label = f"{short}{hint}"

        x_horiz, y_horiz, x_lefm, y_lefm = kt_boundary(dk_th, Y, ds_e)
        ax.plot(x_horiz, y_horiz, color=color, lw=LINE_WIDTH, label=legend_label, zorder=3)
        if x_lefm.size:
            ax.plot(x_lefm, y_lefm, color=color, lw=LINE_WIDTH, zorder=3)

        for _, row in gdf.iterrows():
            c_val = cmap(norm(np.log10(row["cycles"])))
            is_ro = row["runout"] == 1
            ax.scatter(row["sqrt_area"], row["delta_sigma"],
                       color=c_val, marker=marker, s=MARKER_SIZE**2,
                       edgecolors=color if is_ro else "none",
                       linewidths=0.8, zorder=4)
            pt_label = str(row["specimen_id"]) if show_wattage else shorten_id(row["specimen_id"])
            ax.annotate(pt_label,
                        xy=(row["sqrt_area"], row["delta_sigma"]),
                        xytext=(3, 3), textcoords="offset points",
                        fontsize=6, color="dimgrey", zorder=5)
            if is_ro:
                ax.annotate("",
                    xy=(row["sqrt_area"] * ARROW_SCALE, row["delta_sigma"]),
                    xytext=(row["sqrt_area"], row["delta_sigma"]),
                    arrowprops=dict(arrowstyle="->", color=color, lw=0.9))

    style_ax(ax, Y)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    add_colorbar(fig, ax, norm)


def save_fig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_data()

    print("KT data:")
    print(df[["specimen_id", "group", "stress_max_MPa", "delta_sigma",
               "sqrt_area", "cycles", "runout"]].to_string(index=False))

    norm = mcolors.Normalize(vmin=CBAR_LOG_MIN, vmax=CBAR_LOG_MAX)
    cmap = matplotlib.colormaps[CMAP_NAME]

    all_groups = [g for g in GROUP_STYLES if g in df["group"].unique()]

    for Y in Y_VALUES:
        y_tag = f"Y{int(Y * 100):03d}"   # 0.5 → 'Y050', 0.65 → 'Y065'

        # Combined: all groups — show full label with wattage
        fig, ax = plt.subplots(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
                               dpi=DPI, facecolor="white")
        draw_groups(ax, fig, df, all_groups, Y, norm, cmap, show_wattage=True)
        save_fig(fig, OUT_DIR / f"KT_{y_tag}_combined.png")

        # Per wattage — strip wattage prefix from legend
        for watt, groups in WATTAGE_GROUPS.items():
            fig, ax = plt.subplots(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
                                   dpi=DPI, facecolor="white")
            draw_groups(ax, fig, df, groups, Y, norm, cmap, show_wattage=False)
            save_fig(fig, OUT_DIR / f"KT_{y_tag}_{watt}.png")


if __name__ == "__main__":
    main()
