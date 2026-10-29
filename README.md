# CTCV J2056-3014: X-ray timing

Code and data products behind *Precision X-ray Timing Methods and Spin Evolution Modeling for the
Rapidly Rotating Cataclysmic Variable CTCV J2056-3014* (Miao et al.). CTCV J2056-3014 hosts the
second fastest-spinning white dwarf known (P = 29.6 s). Using XMM-Newton, NICER and NuSTAR data
from 2019 to 2024, the paper measures the spin period and its derivative by pulsar timing of
template-matched times of arrival (TOAs):

| quantity | value |
|---|---|
| spin period | P = 29.60968585 ± 3.6e-8 s |
| spin period derivative | Ṗ = (−0.2 ± 1.8) × 10⁻¹⁵ s s⁻¹ |

## Layout

```
notebooks/            one notebook per section of the paper; each cell says which figure or table it makes,
                      and the numbers quoted in the paper are printed at the end
  01_nicer_timing     Z^2 periodograms, folded profiles, hardness ratio, pulse template   (Figures 1-3, 5, 8)
  02_xmm_timing       XMM-Newton periodograms and folded profiles                        (Figures 1, 2)
  03_nustar_timing    NuSTAR periodogram and folded profiles                             (Figures 1, 2)
  05_pulsar_timing    initial frequency, PINT fit for f and f-dot                        (Figures 6-11)
  archive/            earlier exploratory notebooks (not maintained)
pipeline/             scripts that produce the TOAs (NICERsoft photon_toa.py)
src/                  library: data loaders, periodogram and folding tools, statistics, plotting
data/                 TOAs, timing models, pulse template, GTI table (see data/README.md)
paper_figures/        the figures of the paper, written by the notebooks
external/NICERsoft/   git submodule (commit 98f0ff9) with a one-line patch for PINT 1.1.3
```

## Reproducing the paper

1. Environment: Python 3.13 with the packages in `requirements.txt` (PINT 1.1.3, Stingray 2.2.10,
   HENDRICS 8.1). `bash pipeline/setup_nicersoft.sh` fetches NICERsoft and applies the patch.
2. Data: cleaned, barycentered event files in the layout described in `src/config.py`
   (`CTCV_DATA_DIR`). Event screening and barycentering follow Section 2 of the paper
   (HEASoft/NICERDAS, SAS, NuSTARDAS); the HENDRICS event files are made with
   `HENreadevents` and `HENcalibrate`.
3. TOAs: `python pipeline/01_prepare_events.py` (energy bands of Section 3 of the paper) then `bash pipeline/02_make_toas.sh`
   (writes `data/toas.tim` and `data/toas_within_2021A.tim`; both are tracked, so this step can
   be skipped).
4. Figures and numbers: run the notebooks in order (`pipeline/run_all.sh` does all of the above); the
   figures go to `paper_figures/` and each notebook prints the numbers it contributes to the paper.

## Method in brief

- **Periodograms and profiles** (`src/timing_analysis.py`): Z²₁ periodograms with the
  frequency-uncertainty heuristic of the paper's Appendix A; exposure-corrected phase folding.
- **Pulsar timing** (`pipeline/`, `notebooks/05_pulsar_timing.ipynb`): photon phases from
  `photonphase`, a five-Gaussian template fitted with `nitemplate.py` to the 2021A profile
  (`data/template.gauss`), one TOA per observation from `photon_toa.py`, and a PINT fit for f
  and f-dot. The likelihood used for the TOAs is written out in the paper's Appendix B.

## Citation

Miao, E., Bhattacharya, S., Mori, K., Bhattacharyya, S., et al., in preparation.
