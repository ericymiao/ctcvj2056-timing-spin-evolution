#!/bin/bash
# Rebuild everything: TOAs, then the notebooks (figures in paper_figures/, numbers in results/).
set -e
cd "$(dirname "$0")/.."
PY="${CTCV_PYTHON:-/opt/miniconda3/envs/CTCV_PINT/bin/python}"
JX="$(dirname "$PY")/jupyter-execute"

bash pipeline/setup_nicersoft.sh
"$PY" pipeline/01_prepare_events.py
bash pipeline/02_make_toas.sh
for nb in 01_nicer_timing 02_xmm_timing 03_nustar_timing 04_orbital_search 05_pulsar_timing; do
    (cd notebooks && "$JX" --inplace --timeout=-1 --kernel_name=python3 $nb.ipynb)
done
echo "done"
