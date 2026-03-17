# Fatigue Data Visualisation Tool Suite

Two standalone Python scripts for generating publication-quality fatigue plots from tensile fatigue test data:

| Script | Output |
|--------|--------|
| `plot_sn.py` | `sn_output/SN_curve.png` — log-log S-N curve with Basquin regression |
| `plot_kt.py` | 10 PNG files in `kt_output/` — Kitagawa-Takahashi diagrams (El Haddad correction) |

---

## Requirements

```bash
pip install pandas numpy matplotlib
```

Python ≥ 3.8 required.

---

## Directory Structure

Create the following folders alongside the scripts before running:

```
project_root/
├── plot_sn.py
├── plot_kt.py
│
├── sn_input/
│   └── fatigue_data.csv      ← S-N input data
│
├── sn_output/                ← auto-created by plot_sn.py
│
├── kt_input/
│   └── kt_data.csv           ← K-T input data
│
└── kt_output/                ← auto-created by plot_kt.py
```

---

## Input CSV Format

### `sn_input/fatigue_data.csv` — S-N Curve Data

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `specimen_id` | string | — | Unique specimen label, e.g. `GroupA-1` |
| `group` | string | — | Condition label, e.g. `GroupA BD//CD` |
| `stress_max_MPa` | float | MPa | Maximum applied stress σ_max |
| `cycles` | integer | — | Cycles to failure (or run-out limit) |
| `runout` | integer | 0 or 1 | `1` = run-out (not failed); `0` = failed |

**Example:**

```csv
specimen_id,group,stress_max_MPa,cycles,runout
A-H1,GroupA BD//CD,950,2300,0
A-H2,GroupA BD//CD,820,18000,0
A-H7,GroupA BD//CD,430,780000,0
A-H8,GroupA BD//CD,380,1000000,1
A-V1,GroupA BD⊥CD,910,1900,0
...
```

**Supported group names and their plot styles:**

| Group name | Colour | Marker |
|------------|--------|--------|
| `GroupA BD//CD` | Blue | Diamond |
| `GroupA BD⊥CD` | Red | Triangle |
| `GroupB BD//CD` | Green | Square |
| `GroupB BD⊥CD` | Magenta | Circle |

> The script auto-detects groups from the CSV — additional groups will be assigned default styles automatically.

---

### `kt_input/kt_data.csv` — Kitagawa-Takahashi Diagram Data

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `specimen_id` | string | — | Specimen label matching S-N data |
| `group` | string | — | Condition label (must match S-N group names) |
| `stress_max_MPa` | float | MPa | Maximum applied stress σ_max |
| `sqrt_area` | float | μm | √(defect projected area) — Murakami √area parameter from SEM fractography |
| `Nf_cycles` | integer | — | Fatigue life in cycles (used for colourmap colouring) |
| `runout` | integer | 0 or 1 | `1` = run-out; `0` = failed |

**Example:**

```csv
specimen_id,group,stress_max_MPa,sqrt_area,Nf_cycles,runout
A-H1,GroupA BD//CD,950,210,2300,0
A-H2,GroupA BD//CD,820,170,18000,0
A-H8,GroupA BD//CD,380,160,1000000,1
A-V1,GroupA BD⊥CD,910,230,1900,0
...
```

> `sqrt_area` must be in **μm** (micrometres). The script converts internally to SI units for fracture mechanics calculations.

---

## Key Physical Parameters (hardcoded in `plot_kt.py`)

The ΔK_th values must be updated in `plot_kt.py` to match your material and test conditions before running. Locate the `DELTA_K_TH` dictionary near the top of the file:

```python
DELTA_K_TH = {
    "GroupA BD//CD": <value>,   # MPa·m½
    "GroupA BD⊥CD":  <value>,   # MPa·m½
    "GroupB BD//CD": <value>,   # MPa·m½
    "GroupB BD⊥CD":  <value>,   # MPa·m½
}
```

| Parameter | Description |
|-----------|-------------|
| Stress ratio *R* | Set in `R_RATIO` (default: 0.1) |
| Run-out limit | Set in `RUNOUT_LIMIT` (default: 10⁶ cycles) |
| ΔK_th | Threshold stress intensity range per condition — update `DELTA_K_TH` to match your literature/experimental values |

The smooth-specimen fatigue limit **Δσ_e** is derived automatically from each condition's run-out specimen: `Δσ_e = stress_max_MPa(run-out) × (1 − R)`.

---

## Usage

```bash
python plot_sn.py
python plot_kt.py
```

`plot_kt.py` outputs 10 files — 4 individual condition plots + 1 combined plot, for each of the two geometry correction factors Y = 0.50 and Y = 0.65:

```
kt_output/
├── KT_Y050_GroupAH.png      Y=0.50, GroupA BD//CD
├── KT_Y050_GroupAV.png      Y=0.50, GroupA BD⊥CD
├── KT_Y050_GroupBH.png      Y=0.50, GroupB BD//CD
├── KT_Y050_GroupBV.png      Y=0.50, GroupB BD⊥CD
├── KT_Y050_combined.png     Y=0.50, all conditions
├── KT_Y065_GroupAH.png      Y=0.65, GroupA BD//CD
├── KT_Y065_GroupAV.png      Y=0.65, GroupA BD⊥CD
├── KT_Y065_GroupBH.png      Y=0.65, GroupB BD//CD
├── KT_Y065_GroupBV.png      Y=0.65, GroupB BD⊥CD
└── KT_Y065_combined.png     Y=0.65, all conditions
```

---

*ΔK_th values should be sourced from literature or experimental measurements appropriate to your material system.*
