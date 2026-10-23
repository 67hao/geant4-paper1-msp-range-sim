# msp_calc - Tinh Mass Stopping Power proton cho s1..s6 bang Geant4

## Build

```bash
cd msp_geant4
mkdir build && cd build
cmake ..
make -j$(nproc)
```

Se ra 2 chuong trinh: `msp_calc` (MSP + CSDA range) va `range_calc` (Projected
Range that qua tracking Monte Carlo).

Neu ban da co bien moi truong Geant4 (source thisgeant4.sh / geant4make.sh
truoc do), co the bo qua -DGeant4_DIR va chi can:

```bash
cmake ..
make -j$(nproc)
```

De tim dung Geant4_DIR neu can:
```bash
find / -name "Geant4Config.cmake" 2>/dev/null
```

## Chay

```bash
./msp_calc                          # dung energies_kev.txt + xuat msp_results.csv
./msp_calc ../energies_kev.txt out.csv   # tuy chinh duong dan file nang luong / output
```

## Tuy chinh dai nang luong

Sua truc tiep file `energies_kev.txt` (moi dong 1 gia tri keV, dong bat dau
bang `#` se bi bo qua) - khong can build lai code, chi can chay lai
`./msp_calc`.

## Ket qua

Chuong trinh xuat RA 2 FILE CSV, cung layout Energy x s1..s6:

1. `msp_results.csv` - Mass Stopping Power (MeV.cm2/g)
2. `msp_results_csda_range_um.csv` - CSDA Range (um)

```
Energy_keV,s1,s2,s3,s4,s5,s6
15,...,...,...,...,...,...
20,...,...,...,...,...,...
...
```

- MSP: co the doi chieu TRUC TIEP (1-1) voi cot MSP trong
  `msp_range_results.xlsx` (sheet MSP_MeV_cm2_per_g) sinh boi
  `parse_msp_range.py` tu ket qua PHITS.

- CSDA Range: **KHONG PHAI** Projected Range nhu PHITS tinh (sheet
  Projected_Range_um). CSDA range la tong quang duong doc theo TOAN BO
  quy dao cong (ke ca doan bi lech huong do tan xa goc), luon >= Projected
  Range (khoang cach hinh chieu theo truc z). O vung nang luong thap
  (duoi ~100 keV) hai gia tri nay co the lech dang ke do tan xa goc lon.
  Chi nen dung CSDA range lam CAN TREN THAM KHAO, khong so sanh 1-1 truc
  tiep voi cot Range cua parse_msp_range.py.

  Neu can Projected Range THAT su de so sanh 1-1, phai chay tracking Monte
  Carlo that (track nhieu proton, ghi lai toa do z luc dung, lay trung
  binh) - DA DUOC THEM vao project qua chuong trinh `range_calc` (xem ben
  duoi).

## ⭐ sim_msp_range - TINH CA MSP VA RANGE BANG SIMULATION (khong dung G4EmCalculator)

Day la chuong trinh CHINH duoc khuyen dung. Khac voi `msp_calc`/`range_calc`
(dung G4EmCalculator tra bang cho MSP/CSDA range), `sim_msp_range` tinh CA
HAI DAI LUONG TRUC TIEP TU TRACKING MONTE CARLO THAT, ho tro CA proton LAN
alpha:

- **MSP(entrance)**: cong don nang luong hap thu cua hat so cap trong 1 lop
  mong sat be mat (be day tu dong = 1% range uoc luong tu 1 luot calibration
  nhanh truoc do), chia cho be day lop va mat do -> MeV.cm2/g. Chinh la cach
  PHITS/`parse_msp_range.py` da lam, chi khac nguon du lieu (Geant4 track
  thay vi file .out).
- **Projected Range**: trung binh toa do z luc hat dung han qua nhieu hat -
  giong het `range_calc` (van dung nghia, khong doi).

### Chay tu dong ca 2 hat x 6 vat lieu:

```bash
cd build
../run_all_sim.sh ../energies_kev.txt 3000 200
```

(3000 = so hat/nang luong o main pass, 200 = so hat o calibration pass de
uoc luong be day lop entrance)

Ket qua: 4 file CSV, layout `Energy_keV,s1,...,s6`:
- `MSP_proton.csv`, `MSP_alpha.csv`
- `Range_proton.csv`, `Range_alpha.csv`

### Chay rieng le (vi du kiem tra nhanh 1 truong hop):

```bash
./sim_msp_range proton s3 ../energies_kev.txt test_proton_s3 3000 200
./sim_msp_range alpha  s6 ../energies_kev.txt test_alpha_s6  3000 200
```

### Luu y

- Alpha co range NGAN HON proton dang ke o cung nang luong (do khoi luong
  lon hon, dien tich z=2 -> stopping power cao hon nhieu) - `TARGET_LENGTH`
  (5mm) da du du cho ca 2 hat trong dai 15 keV - 15 MeV, nhung neu doi dai
  nang luong cao hon nhieu (vd >50 MeV) can kiem tra lai canh bao "escaped".
- MSP tinh theo phuong phap nay van la XAP XI (be day lop entrance huu han,
  khong phai gioi han Δz->0 that su) - giong het gioi han da neu trong
  `parse_msp_range.py` ban dau, chi khac la be day duoc TU DONG dieu chinh
  theo tung nang luong (thay vi co dinh nhu mesh PHITS truoc day) nen chinh
  xac hon o vung nang luong thap.
- Tang `calibEvents` (vd 500) neu thay range uoc luong dao dong nhieu giua
  cac lan chay (do thong ke it hat).
- Tang `nEvents` de giam sai so thong ke cho ca MSP lan Range.

## (Legacy) range_calc - chi tinh rieng Projected Range bang tracking (proton)

Van con giu trong project de tham khao/doi chieu - `sim_msp_range` da bao
gom ca chuc nang nay (va them ca alpha), nen thong thuong CHI CAN dung
`sim_msp_range` la du, khong bat buoc phai chay `range_calc` nua.

Chuong trinh `range_calc` ban NHIEU proton (mac dinh 3000/nang luong) vao
MOT khoi vat lieu duy nhat, ghi lai toa do z LUC HAT DUNG HAN, lay trung
binh -> day la Projected Range DUNG NGHIA, so sanh 1-1 duoc voi cot
`Projected_Range_um` trong `msp_range_results.xlsx` (PHITS).

Chay RIENG cho tung vat lieu (de tranh hieu ung bien khi nhieu vat lieu
dat canh nhau trong 1 geometry):

```bash
./range_calc s1 ../energies_kev.txt range_s1.csv 3000
./range_calc s2 ../energies_kev.txt range_s2.csv 3000
# ... hoac tu dong ca 6 vat lieu + gop ket qua bang script co san:
cd ..    # ve thu muc msp_geant4 (co chua run_all_ranges.sh)
cd build
../run_all_ranges.sh ../energies_kev.txt 3000 projected_range_results.csv
```

Ket qua `projected_range_results.csv` co layout `Energy_keV,s1,...,s6`,
don vi um.

### Tham so co the can chinh

- **nEvents** (mac dinh 3000): tang len (vd 10000-20000) de giam sai so
  thong ke, doi lai chay lau hon. Do lech chuan cua trung binh giam theo
  can bac hai N.
- **TARGET_LENGTH** trong `range_calc.cc` (mac dinh 2 mm): neu chuong
  trinh in canh bao "escaped" (hat thoat khoi khoi vat lieu truoc khi
  dung), tang gia tri nay len va build lai.
- **SetDefaultCutValue** trong `PhysicsList` cua `range_calc.cc` (mac dinh
  1 um): Geant4 co gioi han duoi ve energy cut (~990 eV), nen o E=15 keV
  (range du kien chi vai chuc nm) do phan giai co the khong hoan hao -
  day la gioi han vat ly cua phuong phap Monte Carlo o vung nang luong rat
  thap, khong phai loi code. Neu can chinh xac hon o vung nay, can giam
  cut xuong nua (nhung se lam cham dang ke).

## Ghi chu ky thuat

- Dung `G4EmStandardPhysics_option4` (models chinh xac nhat vung nang luong
  thap, dua tren Livermore/Penelope) - phu hop dai 15 keV - 15 MeV cua ban.
- `G4EmCalculator::ComputeTotalDEDX()` tra ve dE/dx TRUC TIEP tu physics
  model (khong chay Monte Carlo transport that su) -> khong co nhieu thong
  ke, chi can `BeamOn(1)` mot lan de build physics table.
- MSP = dE/dx (MeV/cm) / mat do (g/cm3), giong het cong thuc dung trong
  `parse_msp_range.py` (PHITS).
