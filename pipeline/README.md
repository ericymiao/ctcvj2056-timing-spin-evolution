# pipeline/

Scripts that produce the times of arrival (TOAs) used by `notebooks/05_pulsar_timing.ipynb`.

```
bash pipeline/setup_nicersoft.sh     # NICERsoft submodule (commit 98f0ff9) + one-line PINT 1.1.3 patch
python pipeline/01_prepare_events.py # 0.4-2 keV, count-rate-cut NICER 2021A file for the within-observation TOAs
bash pipeline/02_make_toas.sh        # photon_toa.py for every observation -> data/toas.tim, data/toas_within_2021A.tim
bash pipeline/run_all.sh             # the above plus every notebook
```

Inputs are the cleaned, barycentered event files (see `src/config.py` for the expected layout;
set `CTCV_DATA_DIR`). Outputs that are tracked in the repository: `data/toas.tim` and
`data/toas_within_2021A.tim`.
