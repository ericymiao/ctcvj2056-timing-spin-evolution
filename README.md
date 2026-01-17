# CTCV J2056-3014 X-ray Timing Analysis

X-ray timing analysis of the cataclysmic variable **CTCV J2056-3014** using NICER, XMM-Newton, and NuSTAR data.

## Overview

**CTCV J2056-3014** is a cataclysmic variable system with a white dwarf accreting material from a companion star. This project detects X-ray pulsations at the WD spin period (~29.6 s) and measures the spin evolution (period derivative, Ṗ) to constrain accretion physics and WD mass.

**Key Results:**
- Spin period: P = 29.60968584 ± 0.00000008 s
- Current Ṗ constraint: (0.0 ± 1.7) × 10⁻¹⁵ s/s
- Hardness ratio anticorrelated with flux (5.5σ in 2024 data)

## Project Structure

```
├── src/                        # Python modules
│   ├── timing_analysis.py      # Z² periodograms, phase folding
│   ├── statistics.py           # F-tests, correlation analysis
│   ├── plotting.py             # Visualization
│   └── demodulation.py         # Orbital motion correction
├── notebooks/
│   ├── 01_nicer_timing.ipynb   # NICER pulsation analysis
│   ├── 02_xmm_timing.ipynb     # XMM-Newton timing
│   ├── 03_nustar_timing.ipynb  # NuSTAR timing
│   ├── 04_pint_analysis.ipynb  # Precision timing with PINT
│   ├── 05_pdot_modeling.ipynb  # Accretion torque modeling
│   └── 06_demodulation.ipynb   # Orbital demodulation
├── data/                       # TOAs, timing models
├── figures/                    # Publication figures
└── wd_mass_radius/             # WD mass-radius tables
```

## Notebooks

1. **01-03**: Pulsation analysis for each instrument (Z² periodograms, phase folding, hardness ratios)
2. **04_pint_analysis**: Precision timing to measure spin frequency and derivative
3. **05_pdot_modeling**: Theoretical Ṗ from accretion torque theory
4. **06_demodulation**: Orbital motion correction to optimize pulsation detection

## Author

Eric Miao
