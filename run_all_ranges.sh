#!/bin/bash
# Chay range_calc cho ca 6 vat lieu s1..s6 roi gop ket qua thanh 1 file CSV
# layout Energy_keV x s1..s6 (giong layout msp_results.csv / parse_msp_range.py)
#
# Dung: ./run_all_ranges.sh [energies_kev.txt] [nEvents] [output.csv]

set -e

ENERGY_FILE="${1:-energies_kev.txt}"
NEVENTS="${2:-3000}"
OUTPUT="${3:-projected_range_results.csv}"

MATERIALS="s1 s2 s3 s4 s5 s6"

echo "Dang chay range_calc cho ${MATERIALS} (nEvents=${NEVENTS} moi nang luong)..."
echo "CANH BAO: viec nay co the mat vai phut den vai chuc phut tuy so nang luong x vat lieu x nEvents."

for mat in $MATERIALS; do
    echo ""
    echo "=== Vat lieu: $mat ==="
    ./range_calc "$mat" "$ENERGY_FILE" "range_${mat}.csv" "$NEVENTS"
done

echo ""
echo "Dang gop ket qua vao $OUTPUT ..."

python3 - "$OUTPUT" $MATERIALS <<'PYEOF'
import sys, csv

output_path = sys.argv[1]
materials = sys.argv[2:]

data = {}   # energy_kev -> {mat: range_um}
energies = []

for mat in materials:
    path = f"range_{mat}.csv"
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            e = float(row["Energy_keV"])
            if e not in data:
                data[e] = {}
                energies.append(e)
            data[e][mat] = row["Range_um"]

energies = sorted(set(energies))

with open(output_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Energy_keV"] + materials)
    for e in energies:
        row = [e] + [data[e].get(m, "") for m in materials]
        writer.writerow(row)

print(f"Da gop xong: {output_path}")
PYEOF

echo "Hoan tat. Xem file: $OUTPUT"
