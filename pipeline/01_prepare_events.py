"""Prepare the event file used for the within-observation TOAs of NICER 2021A.

The five TOAs of the paper are made directly from the cleaned, barycentered source event files
(pipeline/02_make_toas.sh). The three within-observation TOAs that refine the initial frequency
use the 2021A observation restricted to 0.4-2 keV (PI 40-200) with a count-rate cut (5 s bins,
rate below 9 c/s; GTI table in data/gti/). This script writes that file to data/events/.

Run:  python pipeline/01_prepare_events.py
"""

import os
import sys
import datetime

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import config
from events import read_events, check_time_header, energy_filter_events, apply_gti, read_gti_table

os.makedirs(config.EVENTS_DIR, exist_ok=True)
today = datetime.date.today().isoformat()

# Every input must be barycentered (TDB, SOLARSYSTEM)
for key, path in config.EVENT_FILES.items():
    hdul, events, gti = read_events(path)
    info = check_time_header(events.header)
    print(f'{key:12s} {os.path.basename(path):32s} {events.header["NAXIS2"]:7d} events, GTI {np.sum(gti[:, 1] - gti[:, 0]):.1f} s, {info["TIMESYS"]}/{info["TIMEREF"]}')
    hdul.close()

# NICER 2021A, 0.4-2 keV
soft = os.path.join(config.EVENTS_DIR, 'nicer2021A_0.4-2keV.evt')
n = energy_filter_events(config.EVENT_FILES['nicer2021A'], soft, 40, 200,
                         [f'pipeline/01_prepare_events.py {today}: PI 40-200 (0.4-2 keV) selection'])
print(f'{os.path.basename(soft)}: {n} events')

# Count-rate cut (GTI table) for the within-observation TOAs
gti, meta = read_gti_table(config.CRCUT_GTI)
cut = os.path.join(config.EVENTS_DIR, 'nicer2021A_0.4-2keV_crcut.evt')
n = apply_gti(soft, cut, gti, [f'pipeline/01_prepare_events.py {today}: count-rate cut GTI from data/gti/nicer_2021A_crcut.ecsv '
                               f'({len(gti)} intervals, {np.sum(gti[:, 1] - gti[:, 0]):.1f} s)'])
print(f'{os.path.basename(cut)}: {n} events, {len(gti)} GTI rows, {np.sum(gti[:, 1] - gti[:, 0]):.1f} s')
