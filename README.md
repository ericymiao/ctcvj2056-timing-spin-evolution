# CTCV J2056-3014 X-ray Timing Analysis

X-ray timing analysis of the cataclysmic variable **CTCV J2056-3014**, investigating pulsations, spin evolution, and orbital characteristics using data from NICER, XMM-Newton, and NuSTAR observatories.

## Scientific Background

**CTCV J2056-3014** is a cataclysmic variable (CV) system consisting of:
- A **white dwarf (WD)** primary accreting material from a companion star
- An **accretion disk** channeling material onto the WD surface
- X-ray emission from the accretion shock near the WD poles

This project focuses on detecting and characterizing **X-ray pulsations** at the WD spin period (~29.6 seconds) and measuring the **spin evolution** (period derivative, Ṗ) to understand accretion physics.

### Key Scientific Questions

1. **Is the WD spinning up or down?** Accreting material carries angular momentum that should spin up the WD. Measuring Ṗ tests accretion torque theory.
2. **What is the WD mass?** The spin-up rate depends on the WD mass and moment of inertia through the mass-radius relation.
3. **How does the X-ray spectrum vary with spin phase?** Hardness ratio analysis reveals whether absorption or emission geometry drives the pulsations.

## Observations

| Observatory | ObsID | Date | Exposure | Energy Band |
|-------------|-------|------|----------|-------------|
| NICER | 4592010200 | 2021 | 47.6 ks | 0.3-7 keV |
| NICER | 4592010300 | 2021 | 18.5 ks | 0.3-10 keV |
| NICER | 7656010100 | 2024 | 40.1 ks | 0.3-10 keV |
| XMM-Newton | Multiple | 2019-2024 | Various | 0.15-15 keV |
| NuSTAR | Multiple | Various | Various | 3-79 keV |

## Analysis Methods

### Z² Periodogram

The Z² statistic is a Rayleigh-based test for periodicity in X-ray event data. For each trial frequency, photon arrival phases are computed and the statistic measures phase coherence:

```
Z²_n = (2/N) * Σ[k=1 to n] [(Σcos(kφ))² + (Σsin(kφ))²]
```

Under the null hypothesis (no pulsation), Z²_n follows a χ² distribution with 2n degrees of freedom. Higher harmonics (n > 1) can detect non-sinusoidal pulse shapes.

**Key Results:**
- Spin period: **P = 29.60968584 ± 0.00000008 s** (combined 2021 data)
- Spin frequency: **F0 = 0.03377273 Hz**
- Highly significant pulsations detected in all observations

### Phase Folding

Events are binned by arrival phase (modulo the spin period) to build pulse profiles. The analysis includes:
- **GTI-aware exposure correction** for accurate count rates
- **Background subtraction** using SCORPEON model estimates (NICER)
- **Energy-resolved profiles** in soft (0.3-2 keV) and hard (2-10 keV) bands

### Hardness Ratio Analysis

The hardness ratio HR = (H-S)/(H+S) traces spectral variations with pulse phase:
- **F-tests** compare constant vs. sinusoidal HR models
- **Phase correlation tests** determine if HR is correlated or anticorrelated with flux

**Key Finding:** The 2024 observation shows **significant HR variability** (5.5σ) that is **anticorrelated** with flux—the spectrum is harder when the source is fainter, suggesting absorption by intervening material (e.g., accretion curtain) or geometric occultation effects.

### PINT Pulsar Timing

The PINT package performs precision timing to measure F0 (spin frequency) and F1 (frequency derivative):
- **TOA fitting**: Minimizes residuals between observed and model-predicted pulse arrival times
- **Parameter correlation**: F0 and F1 are highly correlated (r ~ -0.92) over limited baselines
- **Future projections**: Simulations predict when F1 will become detectable (>3σ)

**Current Constraints:**
- F1 = (0.0 ± 2.0) × 10⁻¹⁸ Hz/s
- Equivalent Ṗ = (0.0 ± 1.7) × 10⁻¹⁵ s/s
- Longer baseline needed to break F0-F1 degeneracy

### Spin Evolution Modeling

The expected period derivative from accretion torque theory:

```
Ṗ = -Ṁ × R_m² × P / I
```

where:
- **Ṁ** = mass accretion rate (2.3-9.0 × 10¹⁴ g/s from Salcedo et al.)
- **R_m** = magnetospheric radius (assumed = corotation radius for spin equilibrium)
- **I** = WD moment of inertia (from mass-radius relation)

For the constrained WD mass range (0.7-1.0 M☉), expected |Ṗ| ~ 10⁻¹⁵ to 10⁻¹⁴ s/s.

## Project Structure

```
.
├── src/                        # Python analysis modules
│   ├── timing_analysis.py      # Z² periodograms, phase folding
│   ├── statistics.py           # F-tests, correlation analysis
│   ├── plotting.py             # Bokeh visualization
│   └── demodulation.py         # Orbital motion correction
├── notebooks/
│   ├── 01_nicer_timing.ipynb   # NICER pulsation analysis
│   ├── 02_xmm_timing.ipynb     # XMM-Newton timing
│   ├── 03_nustar_timing.ipynb  # NuSTAR high-energy timing
│   ├── 04_pint_analysis.ipynb  # Precision timing with PINT
│   ├── 05_pdot_modeling.ipynb  # Accretion torque modeling
│   └── 06_demodulation.ipynb   # Orbital demodulation search
├── data/                       # Processed lightcurves, TOAs
├── figures/                    # Publication figures
└── wd_mass_radius/             # WD equation of state tables
```

## Notebook Descriptions

### 01_nicer_timing.ipynb
NICER-specific timing analysis including:
- Loading barycentered event files processed with HENDRICS
- Energy filtering (0.3-7/10 keV) and GTI handling
- Z² periodogram computation and frequency uncertainty estimation
- Background-subtracted folded profiles with SCORPEON model
- Hardness ratio F-tests and phase correlation analysis

### 02_xmm_timing.ipynb & 03_nustar_timing.ipynb
Multi-instrument confirmation of pulsations with different energy coverage and systematics.

### 04_pint_analysis.ipynb
Precision timing model fitting:
- TOA generation from folded profiles
- Weighted least-squares fitting for F0, F1
- Parameter correlation and error ellipse visualization
- Monte Carlo simulations for future observation planning
- Detection threshold predictions for F1 measurement

### 05_pdot_modeling.ipynb
Theoretical framework for interpreting Ṗ measurements:
- Accretion torque calculations as function of WD mass
- Mass-radius relation from Hamada & Salpeter models
- Inferred B-field under spin equilibrium assumption
- Constraints on WD mass from future Ṗ detection

### 06_demodulation.ipynb
Orbital motion correction to optimize pulsation detection:
- Grid search over orbital period, amplitude, phase offset
- Timing delay model: δt = A × sin(2π(t-φ)/T)
- Z² power maximization in demodulated event times

## Installation

```bash
git clone https://github.com/ericymiao/ctcvj2056-timing-spin-evolution.git
cd ctcvj2056-timing-spin-evolution
pip install -r requirements.txt
```

### Dependencies

- **numpy, scipy, matplotlib** - Numerical computing
- **astropy** - Time handling, units, coordinates
- **stingray** - X-ray timing (Z² search, event lists, lightcurves)
- **hendrics** - High-energy event processing
- **pint-pulsar** - Precision pulsar timing
- **bokeh** - Interactive visualization
- **nifty-ls** - Fast Lomb-Scargle periodograms

## Usage

Run notebooks in sequence:

```bash
cd notebooks
jupyter notebook 01_nicer_timing.ipynb
```

Import modules for custom analysis:

```python
from src.timing_analysis import z2_periodogram, fold_events, get_frequency_uncertainty
from src.statistics import calculate_hardness_ratio, test_hardness_ratio_variability
from src.plotting import plot_folded_lightcurve, plot_hardness_ratio
```

## References

- Salcedo et al. - CTCV J2056-3014 system parameters
- Hamada & Salpeter (1961) - WD mass-radius relation
- Buccheri et al. (1983) - Z² periodogram method

## Author

Eric Miao

## License

Research code - please contact the author for usage permissions.
