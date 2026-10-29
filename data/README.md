# data/

| file | content |
|---|---|
| `toas.tim` | the five TOAs of the paper (one per observation, 0.3-10 keV; 0.3-7 keV for NICER 2021A), tempo2 format; made by `pipeline/02_make_toas.sh` |
| `toas_within_2021A.tim` | three TOAs within the NICER 2021A observation, used to refine the initial frequency |
| `template.gauss` | the five-Gaussian pulse template fitted with NICERsoft `nitemplate.py` to the folded 2021A 0.4-2 keV profile (60 profile bins) |
| `nicer2021A_template_phases.npy` | rotational phases of the 2021A 0.4-2 keV events (PINT `photonphase`, linear model f = 0.033772731 Hz, PEPOCH 59500) that the template was fitted to |
| `photon_toa_initial_model.par` | timing model used by `photon_toa.py` for phase assignment (f = 0.03377273 Hz, PEPOCH 59420, DE440) |
| `timing_model_z2.par`, `timing_model_resid_slope.par`, `timing_model_centerpeak.par` | initial models of the three frequency-refinement stages (Section 5 of the paper) |
| `timing_model_final.par`, `timing_model_final_f0only.par` | the fitted ephemeris of Table 4 (two-parameter fit) and the adopted f-only fit, written by notebook 05 |
| `test_timing_model.par`, `CTCV_leftpeak*.par`, `CTCV_rightpeak*.par` | central and off-centre alias frequencies tested with the within-observation TOAs |
| `gti/nicer_2021A_crcut.ecsv` | good time intervals of NICER 2021A after the count-rate cut (5 s bins below 9 c/s) |
| `events/` | event files written by `pipeline/01_prepare_events.py` (not tracked) |
