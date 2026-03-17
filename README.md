# Fatigue Data Visualisation Tool Suite

**AUT Master's Dissertation — Fatigue Analysis of LPBF Co-29Cr-6Mo Alloy**

Two standalone Python scripts for generating publication-quality fatigue plots:

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
| `specimen_id` | string | — | Unique specimen label, e.g. `180H1` |
| `group` | string | — | Condition label, e.g. `180W BD//CD` |
| `stress_max_MPa` | float | MPa | Maximum applied stress σ_max |
| `cycles` | integer | — | Cycles to failure (or run-out limit) |
| `runout` | integer | 0 or 1 | `1` = run-out (not failed); `0` = failed |

**Example:**

```csv
specimen_id,group,stress_max_MPa,cycles,runout
180H1,180W BD//CD,900,1000,0
180H2,180W BD//CD,800,10000,0
180H7,180W BD//CD,400,100000,0
180H8,180W BD//CD,350,1000000,1
180V1,180W BD⊥CD,900,900,0
...
```

**Supported group names and their plot styles:**

| Group name | Colour | Marker |
|------------|--------|--------|
| `180W BD//CD` | Blue | Diamond |
| `180W BD⊥CD` | Red | Triangle |
| `320W BD//CD` | Green | Square |
| `320W BD⊥CD` | Magenta | Circle |

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
180H1,180W BD//CD,900,220,1000,0
180H2,180W BD//CD,800,150,10000,0
180H8,180W BD//CD,350,150,1000000,1
180V1,180W BD⊥CD,900,220,900,0
...
```

> `sqrt_area` must be in **μm** (micrometres). The script converts internally to SI units for fracture mechanics calculations.

---

## Key Physical Parameters (hardcoded in `plot_kt.py`)

| Parameter | Value | Source |
|-----------|-------|--------|
| Stress ratio *R* | 0.1 | — |
| Run-out limit | 10⁶ cycles | — |
| ΔK_th — 180W BD//CD | 6.2 MPa·m½ | Anuar et al. (2021) |
| ΔK_th — 180W BD⊥CD | 4.9 MPa·m½ | Anuar et al. (2021) |
| ΔK_th — 320W BD//CD | 6.7 MPa·m½ | Anuar et al. (2021) |
| ΔK_th — 320W BD⊥CD | 5.2 MPa·m½ | Anuar et al. (2021) |

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
├── KT_Y050_180H.png      Y=0.50, 180W BD//CD
├── KT_Y050_180V.png      Y=0.50, 180W BD⊥CD
├── KT_Y050_320H.png      Y=0.50, 320W BD//CD
├── KT_Y050_320V.png      Y=0.50, 320W BD⊥CD
├── KT_Y050_combined.png  Y=0.50, all conditions
├── KT_Y065_180H.png      Y=0.65, 180W BD//CD
├── KT_Y065_180V.png      Y=0.65, 180W BD⊥CD
├── KT_Y065_320H.png      Y=0.65, 320W BD//CD
├── KT_Y065_320V.png      Y=0.65, 320W BD⊥CD
└── KT_Y065_combined.png  Y=0.65, all conditions
```

---

*Auckland University of Technology — Department of Mechanical Engineering*
*ΔK_th values: Anuar et al. (2021), J. Mech. Behav. Biomed. Mater. 123, 104741*
