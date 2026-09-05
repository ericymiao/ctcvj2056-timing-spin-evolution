# CTCV J2056-3014: X-ray timing

Code and data products behind *Precision X-ray Timing and Spin Evolution of the Rapidly Rotating
Cataclysmic Variable CTCV J2056-3014* (Miao et al.). CTCV J2056-3014 hosts the second
fastest-spinning white dwarf known (P = 29.6 s). Using XMM-Newton, NICER and NuSTAR data from
2019 to 2024, the paper measures the orbital period in the X-ray band from the Rømer delay of
the spin pulses and the spin period and its derivative by pulsar timing of template-matched
times of arrival (TOAs):

| quantity | value |
|---|---|
| orbital period | P_orb = 6223 ± 20 s |
| spin period | P = 29.60968584 ± 3.7e-8 s |
| spin period derivative | Ṗ = (0.0 ± 1.8) × 10⁻¹⁵ s s⁻¹ |

## Layout

```
notebooks/            one notebook per section of the paper; each cell says which figure or table it makes
  01_nicer_timing     Z^2 periodograms, folded profiles, hardness ratio, pulse template   (Figures 1-3, 6, 9)
  02_xmm_timing       XMM-Newton periodograms and folded profiles                        (Figures 1, 2)
  03_nustar_timing    NuSTAR periodogram and folded profiles                             (Figures 1, 2)
  04_orbital_search   Roemer-delay grid search and control study                         (Figure 4)
  05_pulsar_timing    initial frequency, PINT fit for f and f-dot                        (Figures 7-12)
  archive/            earlier exploratory notebooks (not maintained)
pipeline/             scripts that produce the TOAs (NICERsoft photon_toa.py)
src/                  library: data loaders, periodogram and folding tools, statistics, plotting
data/                 TOAs, timing models, pulse template, GTI table (see data/README.md)
paper_figures/        the figures of the paper, written by the notebooks
results/              the numbers quoted in the paper (JSON), written by the notebooks
external/NICERsoft/   git submodule (commit 98f0ff9) with a one-line patch for PINT 1.1.3
```

## Reproducing the paper

1. Environment: Python 3.13 with the packages in `requirements.txt` (PINT 1.1.3, Stingray 2.2.10,
   HENDRICS 8.1). `bash pipeline/setup_nicersoft.sh` fetches NICERsoft and applies the patch.
2. Data: cleaned, barycentered event files in the layout described in `src/config.py`
   (`CTCV_DATA_DIR`). Event screening and barycentering follow Section 2 of the paper
   (HEASoft/NICERDAS, SAS, NuSTARDAS); the HENDRICS event files are made with
   `HENreadevents` and `HENcalibrate`.
3. TOAs: `python pipeline/01_prepare_events.py` then `bash pipeline/02_make_toas.sh`
   (writes `data/toas.tim` and `data/toas_within_2021A.tim`; both are tracked, so this step can
   be skipped).
4. Figures and numbers: run the notebooks in order (`pipeline/run_all.sh` does all of the above).
   The orbital search caches its two grid searches in `data/*.npy` (about 20 minutes to rebuild).

## Method in brief

- **Periodograms and profiles** (`src/timing_analysis.py`): Z²₁ periodograms with the
  frequency-uncertainty heuristic of the paper's Appendix A; exposure-corrected phase folding.
- **Orbital search** (`src/demodulation.py`): each event time is corrected by
  −A sin(2πt/P_orb − φ) and Z²₁ at the spin period is maximised over a (P_orb, A, φ) grid; the
  false-alarm probability comes from the same search at control periods.
- **Pulsar timing** (`pipeline/`, `notebooks/05_pulsar_timing.ipynb`): photon phases from
  `photonphase`, a five-Gaussian template fitted with `nitemplate.py` to the 2021A profile
  (`data/template.gauss`), one TOA per observation from `photon_toa.py`, and a PINT fit for f
  and f-dot. The likelihood used for the TOAs is written out in the paper's Appendix B.

## Citation

Miao, E., Bhattacharya, S., Mori, K., Bhattacharyya, S., et al., in preparation.
