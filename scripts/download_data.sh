#!/usr/bin/env bash
# Download all raw data into data/raw/ (or into $1). Keeps only Bangladesh rows of large international files.
set -euo pipefail
OUT="${1:-$(cd "$(dirname "$0")/.." && pwd)/data/raw}"
mkdir -p "$OUT"/national "$OUT"/mental_health_2021 "$OUT"/women_panel_2019_2020

# --- National series ---------------------------------------------------------
curl -sL https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv \
  | awk -F, 'NR==1 || ($1=="BGD" && $4<="2022-12-31")' > "$OUT/national/owid_bgd.csv"   # study period ends 2022
curl -sL https://raw.githubusercontent.com/OxCGRT/covid-policy-dataset/main/data/OxCGRT_compact_national_v1.csv \
  | awk -F, 'NR==1 || $2=="BGD"' > "$OUT/national/oxcgrt_bgd.csv"
TMP=$(mktemp -d)
for y in 2020 2021 2022; do
  curl -sL "https://www.gstatic.com/covid19/mobility/${y}_BD_Region_Mobility_Report.csv" -o "$TMP/$y.csv"
done
{ cat "$TMP/2020.csv"; tail -n +2 "$TMP/2021.csv"; tail -n +2 "$TMP/2022.csv"; } \
  | awk -F, 'NR==1 || $3==""' > "$OUT/national/mobility_BD_national.csv"   # national rows only
rm -rf "$TMP"

# --- Mendeley Data (CC BY 4.0), via the public API -----------------------------
mendeley() {  # $1 dataset id, $2 target dir, $3 python regex of files to keep
  curl -s "https://data.mendeley.com/public-api/datasets/$1" | python3 -c "
import sys, json, re
for f in json.load(sys.stdin)['files']:
    if re.search(r'$3', f['filename']): print(f['content_details']['download_url'] + '\t' + f['filename'])" |
  while IFS=$'\t' read -r url fn; do curl -sL "$url" -o "$2/$fn"; done
}
mendeley 6szgbkpzn8 "$OUT/mental_health_2021" '\.sav$'
mv "$OUT"/mental_health_2021/*.sav "$OUT/mental_health_2021/mh2021.sav"
mendeley hxf7nb964d "$OUT/women_panel_2019_2020" '\.(csv|docx|xlsx)$'
echo "Raw data written to $OUT"
