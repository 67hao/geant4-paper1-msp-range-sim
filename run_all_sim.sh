#!/bin/bash
# Chay sim_msp_range cho CA proton VA alpha x 6 vat lieu s1..s6, gop ket qua
# thanh 4 file CSV cuoi cung:
#   MSP_proton.csv, MSP_alpha.csv, Range_proton.csv, Range_alpha.csv
# Moi file layout: Energy_keV,s1,s2,s3,s4,s5,s6
#
# Dung: ./run_all_sim.sh [energies_kev.txt] [nEvents] [calibEvents]

set -e

ENERGY_FILE="${1:-energies_kev.txt}"
NEVENTS="${2:-3000}"
CALIB="${3:-200}"

MATERIALS="s1 s2 s3 s4 s5 s6"
PARTICLES="proton alpha"

echo "Dang chay sim_msp_range cho ${PARTICLES} x ${MATERIALS}"
echo "(nEvents=${NEVENTS}, calibEvents=${CALIB} moi nang luong)"
echo "CANH BAO: viec nay co the mat kha lau (2 hat x 6 vat lieu x so nang luong x nEvents)."

for particle in $PARTICLES; do
    for mat in $MATERIALS; do
        echo ""
        echo "=== $particle / $mat ==="
        ./sim_msp_range "$particle" "$mat" "$ENERGY_FILE" "sim_${particle}_${mat}" "$NEVENTS" "$CALIB"
    done
done

echo ""
echo "Dang gop ket qua..."

python3 - $PARTICLES <<'PYEOF'
import sys, csv

particles = sys.argv[1:]
materials = ["s1", "s2", "s3", "s4", "s5", "s6"]

for particle in particles:
    for kind, colname in [("msp", "MSP_MeV_cm2_per_g"), ("range", "Range_um")]:
        data = {}
        energies = []
        for mat in materials:
            path = f"sim_{particle}_{mat}_{kind}.csv"
            with open(path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    e = float(row["Energy_keV"])
                    if e not in data:
                        data[e] = {}
                        energies.append(e)
                    data[e][mat] = row[colname]
        energies = sorted(set(energies))

        outname = ("MSP_" if kind == "msp" else "Range_") + particle + ".csv"
        with open(outname, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Energy_keV"] + materials)
            for e in energies:
                writer.writerow([e] + [data[e].get(m, "") for m in materials])
        print(f"Da gop xong: {outname}")
PYEOF

echo ""
echo "Hoan tat. 4 file ket qua: MSP_proton.csv, MSP_alpha.csv, Range_proton.csv, Range_alpha.csv"
