#!/bin/bash
# Initialise the NICERsoft submodule (commit 98f0ff9) and apply the one-line patch that
# photon_toa.py needs with PINT 1.1.3 (the include_gps keyword was removed from
# pint.event_toas.get_XMM_TOAs). Idempotent: safe to run again.
set -e
cd "$(dirname "$0")/.."

git submodule update --init external/NICERsoft

# The unpatched file still passes include_gps to get_XMM_TOAs; patch only in that case
if grep -q "include_gps=args.use_gps" external/NICERsoft/scripts/photon_toa.py; then
    echo "Applying pipeline/nicersoft_xmm_pint113.patch"
    (cd external/NICERsoft && git apply ../../pipeline/nicersoft_xmm_pint113.patch)
else
    echo "NICERsoft photon_toa.py already patched"
fi

# Smoke test with the pipeline's Python
PY=$(/opt/miniconda3/envs/CTCV_PINT/bin/python -c "import sys; print(sys.executable)")
"$PY" external/NICERsoft/scripts/photon_toa.py --help > /dev/null && echo "photon_toa.py --help OK ($PY)"
