#!/bin/bash
# Times of arrival with NICERsoft photon_toa.py (submodule external/NICERsoft, PINT 1.1.3).
#
# One TOA per observation from the cleaned, barycentered source event files, using the Gaussian
# template data/template.gauss and the initial model data/photon_toa_initial_model.par; --tint
# is the good time accumulated (whole GTIs) into the TOA and --maxint the maximum time span.
# The three within-observation TOAs of NICER 2021A use the 0.4-2 keV, count-rate-cut file made
# by pipeline/01_prepare_events.py.
#
# Run:  bash pipeline/02_make_toas.sh        (about 10 minutes)
set -e
cd "$(dirname "$0")/.."
DATA_ROOT="${CTCV_DATA_DIR:-/Users/ericmiao/Documents/Kaya Research/CTCV_J2056}"
PY="${CTCV_PYTHON:-/opt/miniconda3/envs/CTCV_PINT/bin/python}"
PHOTON_TOA=external/NICERsoft/scripts/photon_toa.py
TEMPLATE=data/template.gauss
PAR=data/photon_toa_initial_model.par
NICER="$DATA_ROOT/analysis_files_strictfilter"
XMM="$DATA_ROOT/xmm"

# PINT reads the mission table $HEADAS/bin/xselect.mdb to find the XMM event extension; without a
# HEASoft installation, point HEADAS at a directory holding the copy that ships with stingray
STUB="$PWD/external/heasoft_stub"
mkdir -p "$STUB/bin"
[ -f "$STUB/bin/xselect.mdb" ] || cp "$($PY -c 'import os, stingray; print(os.path.join(os.path.dirname(stingray.__file__), "datasets", "xselect.mdb"))')" "$STUB/bin/"
export HEADAS="$STUB"

mkdir -p data/toas
toa() {   # name tint maxint minexp eventfile
    echo "photon_toa $1 (tint=$2 s, maxint=$3 s)"
    "$PY" $PHOTON_TOA --ephem DE440 --tint "$2" --maxint "$3" --minexp "$4" --outfile "data/toas/$1.tim" "$5" $TEMPLATE $PAR > "data/toas/$1.log" 2>&1
    grep '^photon_toa' "data/toas/$1.tim"
}
toa xmm2019_pn    15000  864000  0 "$XMM/0842570101/pn_clean_src_bary.evt"
toa nicer2021A    50000  864000  0 "$NICER/ni4592010200_0mpu7_cl_bary.evt"
toa nicer2021B    30000  864000  0 "$NICER/ni4592010300_0mpu7_cl_bary.evt"
toa xmm2022_pn    15000  864000  0 "$XMM/0902500101/pn_clean_src.evt"
toa nicer2024    100000 8640000  0 "$NICER/ni7656010100_0mpu7_cl_bary.evt"
toa within_2021A  15000  864000 2000 data/events/nicer2021A_0.4-2keV_crcut.evt

# Assemble data/toas.tim (five TOAs) and data/toas_within_2021A.tim (three TOAs, -group obs1)
"$PY" - <<'PYEOF'
import sys, datetime
sys.path.insert(0, 'src')
from events import parse_tim, write_tim
today = datetime.date.today().isoformat()
names = ['xmm2019_pn', 'nicer2021A', 'nicer2021B', 'xmm2022_pn', 'nicer2024']
toas = [t for n in names for t in parse_tim(f'data/toas/{n}.tim')]
write_tim('data/toas.tim', toas, comments=[f'One TOA per observation (photon_toa.py, PINT 1.1.3, DE440), pipeline/02_make_toas.sh {today}',
                                           'template data/template.gauss, model data/photon_toa_initial_model.par'] +
          [f'{n}: {t["flags"]["t"]}, exposure {float(t["flags"]["exposure"]):.1f} s' for n, t in zip(names, toas)])
within = parse_tim('data/toas/within_2021A.tim')
write_tim('data/toas_within_2021A.tim', within, comments=[f'Three TOAs within NICER 2021A (0.4-2 keV, count-rate cut, tint 15 ks), pipeline/02_make_toas.sh {today}'],
          extra_flags={'group': 'obs1'})
print(f'data/toas.tim: {len(toas)} TOAs; data/toas_within_2021A.tim: {len(within)} TOAs')
PYEOF
