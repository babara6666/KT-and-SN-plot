"""
plot_sn_darea.py — S-N Curve with defect area (darea) labels
Same layout as plot_sn.py, but annotates each point with its darea value (µm²).
Points with darea = 0 or NaN are not labelled.
"""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

RUNOUT_LIMIT = 10_000_000
DPI          = 300
FIG_W_MM, FIG_H_MM = 180, 130
EXCEL_PATH   = Path("D:/DownloadD/WAAM/LPBF/SN force.xlsx")
OUT_DIR      = Path(__file__).parent / "sn_output"

GROUP_STYLES = {
    "180W BD//LD": {"color": "blue",    "marker": "D"},
    "180W BD⊥LD": {"color": "red",     "marker": "^"},
    "320W BD//LD": {"color": "green",   "marker": "s"},
    "320W BD⊥LD": {"color": "magenta", "marker": "o"},
}
MARKER_SIZE  = 6
LINE_WIDTH   = 1.2
ARROW_OFFSET = 0.15


def fix_group(g):
    if "//" in g:
        return g
    prefix = g.split("W")[0] + "W"
    return f"{prefix} BD⊥LD"


def load_data():
    df = pd.read_excel(EXCEL_PATH)
    df = df[df["stress_max_MPa"] > 0].dropna(subset=["cycles"]).copy()
    df["group"]  = df["group"].apply(fix_group)
    df["runout"] = (df["cycles"] >= RUNOUT_LIMIT).astype(int)
    df["darea"]  = pd.to_numeric(df["darea"], errors="coerce").fillna(0)
    return df


def annotate_point(ax, row, color):
    """Draw specimen ID above and darea (if valid) below the data point."""
    ax.annotate(row["specimen_id"][3:],
        xy=(row["cycles"], row["stress_max_MPa"]),
        xytext=(4, 4), textcoords="offset points",
        fontsize=6, color=color, zorder=4)
    if row["darea"] > 0:
        ax.annotate(f'{int(row["darea"])}',
            xy=(row["cycles"], row["stress_max_MPa"]),
            xytext=(4, -10), textcoords="offset points",
            fontsize=5.5, color="dimgrey", zorder=4)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_data()

    print("Data loaded:")
    print(df[["specimen_id", "group", "stress_max_MPa", "cycles", "runout", "darea"]].to_string(index=False))

    ordered_groups = [g for g in GROUP_STYLES if g in df["group"].unique()]

    fig, ax = plt.subplots(
        figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
        dpi=DPI, facecolor="white"
    )

    for grp in ordered_groups:
        gdf    = df[df["group"] == grp]
        color  = GROUP_STYLES[grp]["color"]
        marker = GROUP_STYLES[grp]["marker"]
        failed  = gdf[gdf["runout"] == 0]
        runouts = gdf[gdf["runout"] == 1]

        ax.scatter(failed["cycles"], failed["stress_max_MPa"],
                   color=color, marker=marker, s=MARKER_SIZE**2, zorder=3, label=grp)
        for _, row in failed.iterrows():
            annotate_point(ax, row, color)

        for _, row in runouts.iterrows():
            ax.scatter(row["cycles"], row["stress_max_MPa"],
                       color=color, marker=marker, s=MARKER_SIZE**2, zorder=3)
            annotate_point(ax, row, color)
            ax.annotate("",
                xy=(row["cycles"] * 10**ARROW_OFFSET, row["stress_max_MPa"]),
                xytext=(row["cycles"], row["stress_max_MPa"]),
                arrowprops=dict(arrowstyle="->", color=color, lw=LINE_WIDTH))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e4, 2e7)
    ax.set_ylim(200, 1100)
    ax.set_xlabel("Number of Cycles to Failure, $N$", fontsize=9)
    ax.set_ylabel(r"Maximum Stress, $\sigma_\mathrm{max}$ (MPa)", fontsize=9)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_yticks([200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100])
    ax.grid(which="major", color="lightgrey", linestyle="-",  linewidth=0.6, zorder=0)
    ax.grid(which="minor", color="lightgrey", linestyle="--", linewidth=0.3, zorder=0)
    ax.minorticks_on()
    ax.text(0.97, 0.97, "$f$ = 30 Hz,  $R$ = 0.1",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.8))
    ax.legend(fontsize=7, loc="lower left", framealpha=0.9)

    fig.tight_layout()
    out_path = OUT_DIR / "SN_curve_darea.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
