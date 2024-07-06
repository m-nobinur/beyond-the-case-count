#!/usr/bin/env bash
# Reproduce the whole project from raw data to the PDF report.
# Needs R (packages listed in README) and a Python venv: python3 -m venv .venv && .venv/bin/pip install -r python/requirements.txt
set -euo pipefail
cd "$(dirname "$0")"
[ -f data/raw/mental_health_2021/mh2021.sav ] || scripts/download_data.sh
rm -rf output && mkdir -p output
Rscript scripts/01_national.R
Rscript scripts/02_panel.R
Rscript scripts/03_mh_explore.R
Rscript scripts/04_mh_profiles.R
Rscript scripts/05_mh_ml.R
.venv/bin/python python/ml_crosscheck.py
Rscript scripts/06_ml_crosscheck.R
mkdir -p docs/figures   # figures shown in the README
cp output/figures/{fig_policy_mobility,fig_panel_quintile,fig_ml_repeated}.png docs/figures/
# Fixed build date (20 Jul 2024) for the PDF metadata, so rebuilds keep the same date
export SOURCE_DATE_EPOCH=1721481120 FORCE_SOURCE_DATE=1
(cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex > ../output/latex_build.log 2>&1) || { echo "LaTeX failed: see report/main.log"; exit 1; }
cp report/main.pdf Final_Report.pdf
(cd report && latexmk -c > /dev/null 2>&1 || true)
echo "Done: Final_Report.pdf"
