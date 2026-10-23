// ---------------------------------------------------------------------------
// sim_msp_range : tinh CA MSP(entrance) VA Projected Range HOAN TOAN BANG
// MONTE CARLO SIMULATION (khong dung G4EmCalculator), cho PROTON hoac ALPHA,
// tren 1 vat lieu (s1..s6) x dai nang luong tuy chinh.
//
// *** PHIEN BAN 2 - DOI CACH TINH MSP(entrance) ***
// Phien ban truoc dung "track-boundary-crossing": bat thoi diem hat SO CAP
// vuot qua 1 mat phang ao z=layerThickness, noi suy DeltaE theo Z giua 2
// step. Cach do bi 1 loi he thong nghiem trong o vung nang luong THAP (hoac
// vat lieu nang, Z lon): khi tan xa Coulomb nhieu lan (MSC) / tan xa hat
// nhan dan hoi lam quy dao hat cong/ngoan ngoeo manh, PATH LENGTH thuc te
// hat di dai hon nhieu so voi khoang cach Z no tien duoc ("detour factor").
// Neu 1 step bi tan xa manh (di gan nhu ngang/lui) lai la step "vuot qua"
// nguong Z mong, TOAN BO nang luong mat cua ca step do (voi path length dai,
// nen mat nhieu nang luong that) bi gan het vao "mat tai bien" -> DeltaE bi
// thoi phong -> MSP OVERESTIMATE manh va rat nhieu nhieu (dac biet ro o vung
// 15-100 keV proton trong vat lieu giau La, Z=57).
//
// CACH MOI (giong huong PHITS: T-Deposit/T-Track tally theo mesh):
//   - Chia khoi vat lieu thanh N SLAB THAT (G4LogicalVolume rieng, xep lien
//     tiep theo z, tong do day ~ 1.6 x projected range uoc luong).
//   - Tally TONG EDEP CUA MOI HAT (ca hat so cap LAN delta-ray/secondary,
//     khong chi rieng track chinh) roi vao TUNG SLAB (theo copy number cua
//     slab ma step dang o trong, lay tu PreStepPoint touchable).
//   - MSP(entrance) = <Edep trong SLAB DAU TIEN> / nEvents / slabThickness / density.
//   Vi day la edep THAT SU (Geant4 tu tra step tai bien gioi hinh hoc that,
//   khong con phu thuoc quy dao hat so cap co ngoan ngoeo hay khong), loi
//   detour-bias noi tren bien mat hoan toan. Delta-ray sinh gan bien slab co
//   the mang nang luong sang slab ben canh, nhung do la mot HIEU UNG THAT
//   (finite bin resolution), khong phai bug - va co the kiem soat bang cach
//   chon slab du day (xem AUTO_TOTAL_LENGTH_FACTOR / AUTO_NSLABS ben duoi).
//
// Tham so tu dieu chinh theo nang luong - GIO GOM 2 PHA HOAN TOAN TACH BIET:
//   - RANGE-PASS (rieng): ban DAY DU nEvents hat vao 1 khoi lien duy nhat
//     (1 slab dai TARGET_LENGTH_CALIB) - DUNG Y HET geometry cua
//     range_calc.cc - de tinh Range (dinh nghia "z luc dung han" cua track
//     chinh). Day la NGUON DUY NHAT cho Range xuat ra CSV/console
//     (rangeFinalDedicated_um). Khong lien quan gi toi MSP/mesh slab. Pha
//     nay THAY THE hoan toan "Pha 1: CALIBRATION" (calibEvents) cua ban
//     truoc - gio dung thang nEvents thay vi mau nho rieng.
//   - MSP MAIN PASS (mesh slab): dung chinh ket qua Range-pass o tren
//     (rangeEst_um / maxRangeEst_um) de tu dinh kich thuoc mesh N-slab
//     (totalLen = rangeEst x AUTO_TOTAL_LENGTH_FACTOR, chia deu AUTO_NSLABS
//     slab, khong mong hon MIN_SLAB_UM), roi ban lai nEvents hat TREN MESH
//     nay, tally edep tung slab de tinh MSP(entrance) = Edep(slab dau)/
//     nEvents/slabThickness/density. Pha nay KHONG con duoc dung de bao cao
//     Range nua (chi con dung noi bo de auto-grow mesh khi con hat thoat).
//   - Xuat them 1 file profile MSP-theo-do-sau cho MOI muc nang luong de
//     kiem tra hinh dang duong cong Bragg (phai tron, 1 dinh) - dung de
//     validate ket qua, phat hien ngay neu con nhieu/hoi tu chua dat.
//
// Bien dich: da co san trong CMakeLists.txt cung project.
//
// Chay:
//   ./sim_msp_range <particle: proton|alpha> <material: s1..s6>
//                    [energies_kev.txt] [outputPrefix] [nEvents] [calibEvents]
//   vi du:
//   ./sim_msp_range proton s3 ../energies_kev.txt result_proton_s3 3000 200
//   ./sim_msp_range alpha  s3 ../energies_kev.txt result_alpha_s3  3000 200
//
// Dung script run_all_sim.sh de tu dong chay ca 2 hat x 6 vat lieu va gop
// ket qua thanh 4 file CSV (MSP_proton.csv, MSP_alpha.csv, Range_proton.csv,
// Range_alpha.csv), layout Energy_keV x s1..s6.
// ---------------------------------------------------------------------------

#include "Randomize.hh"
#include <ctime>
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHandle.hh"
#include "G4VTouchable.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// Thanh phan + mat do 6 vat lieu (giong het main.cc / range_calc.cc)
// ---------------------------------------------------------------------------
struct ElementFraction { std::string symbol; double wt; };

static const std::map<std::string, std::vector<ElementFraction>> COMPOSITION = {
    {"s1", {{"B",0.1629}, {"O",0.4648}, {"Si",0.0302}, {"Ca",0.0431}, {"La",0.2990}}},
    {"s2", {{"B",0.1329}, {"O",0.4086}, {"Si",0.0266}, {"Ca",0.0379}, {"La",0.3941}}},
    {"s3", {{"B",0.1095}, {"O",0.3644}, {"Si",0.0237}, {"Ca",0.0338}, {"La",0.4688}}},
    {"s4", {{"B",0.0905}, {"O",0.3289}, {"Si",0.0214}, {"Ca",0.0305}, {"La",0.5288}}},
    {"s5", {{"B",0.0750}, {"O",0.2996}, {"Si",0.0195}, {"Ca",0.0278}, {"La",0.5783}}},
    {"s6", {{"B",0.0620}, {"O",0.2752}, {"Si",0.0179}, {"Ca",0.0255}, {"La",0.6195}}},
};

static const std::map<std::string, double> DENSITY_G_CM3 = {
    {"s1", 3.0099}, {"s2", 3.1999}, {"s3", 3.3899},
    {"s4", 3.5799}, {"s5", 3.7699}, {"s6", 3.9599},
};

// Chieu dai khoi vat lieu (z) dung rieng cho PHA CALIBRATION (uoc luong tho
// projected range). PHAI DU LON de chua het range o 15 MeV ca proton lan
// alpha (5mm du du cho ca 2 hat trong dai 15 keV - 15 MeV).
static const G4double TARGET_LENGTH_CALIB = 5.0 * mm;
static const G4double TARGET_XY           = 2.0 * mm;

// --- Tham so mesh slab cho PHA MAIN (do MSP bang edep-tally) --------------
// Tong do day mesh = rangeEst x AUTO_TOTAL_LENGTH_FACTOR (du du 1.6x de
// chua ca straggling - mot so hat co range that dai hon mean range).
static const double AUTO_TOTAL_LENGTH_FACTOR = 1.6;
// He so an toan nhan them vao Z-stop LON NHAT quan sat duoc trong pha
// calibration (khong chi mean) - vi calibEvents chi la mau nho (vd 200 hat),
// straggling tail thuc te tren nEvents (vd 3000 hat) o main pass co the con
// dai hon gia tri max quan sat duoc trong mau nho do. Tang he so nay neu van
// con canh bao "escaped" sau khi da ap dung.
static const double MAX_STOP_SAFETY_MARGIN = 1.6;
// So lan toi da "phong to" mesh va CHAY LAI main pass neu van con hat thoat
// mesh (escaped) - moi lan phong to dua tren DU LIEU THUC cua chinh main
// pass vua chay (nEvents hat), khong con phu thuoc ngoai suy tu calibration
// (chi vai tram hat) nua - giai quyet dung nguyen nhan "extreme-value tail
// cua mau lon (nEvents) luon dai hon mau nho (calibEvents)".
static const int    MAX_MESH_RETRIES = 4;
static const double MESH_GROWTH_FACTOR = 1.5;
// So slab mac dinh muon chia (slab dau tien = vung "entrance" dung de tinh
// MSP). Se tu giam xuong neu slabThickness ra qua mong so voi MIN_SLAB_UM.
static const int    AUTO_NSLABS_DEFAULT = 400;
// San duoi be day 1 slab (um) - LUU Y QUAN TRONG SAU KHI DONG BO CUT XUONG
// 1 NM (tu 1 um truoc do): voi cut=1nm, Geant4 se TAO TRACK RIENG cho hau
// het delta-ray (chi nhung electron co range <1nm moi bi deposit tai cho),
// tuc la SO LUONG delta-ray duoc track rieng (co the "di" ra khoi slab truoc
// khi deposit het nang luong) TANG LEN DANG KE so voi khi cut=1um (luc do
// chi delta-ray range >1um moi duoc track rieng, con lai deposit tai cho).
// Nghia la NGUY CO delta-ray leak qua bien slab (lam MSP entrance bi
// UNDERESTIMATE) O DAY THUC SU CAO HON TRUOC, khac han voi range_calc.cc -
// file do KHONG lam edep-tally theo slab nen cut khong anh huong toi loai
// bug nay. BAT BUOC phai convergence-test lai MIN_SLAB_UM duoi cau hinh
// cut=1nm moi nay (quet vd 0.05/0.1/0.2/0.5 um) truoc khi tin ket qua MSP -
// gia tri 0.05 um ben duoi la ke thua tu thoi cut=1um, CO THE KHONG CON PHU
// HOP nua.
static const double MIN_SLAB_UM = 0.0001;

// ---------------------------------------------------------------------------
// Bien global de SteppingAction biet dang o CHE DO nao (calibration hay
// main/binning) va tong chieu dai khoi hien tai (de bat "thoat khoi vat
// lieu"). Duoc set truoc moi lan BeamOn().
// ---------------------------------------------------------------------------
static G4bool  g_binningActive              = false; // true = dang o main pass (co slab that)
static G4double g_targetTotalLength_internal = 0.0;   // don vi noi bo G4 (mm-scale)
// Vi Geant4 bat buoc World phai dat tai goc toa do, chong slab GIO nam doi
// xung qua goc (tu -totalLen/2 den +totalLen/2) thay vi tu 0 den totalLen
// nhu truoc. g_zEntrance_internal = toa do TUYET DOI (frame Geant4) cua MAT
// VAO slab dau tien (= -totalLen/2), dung de:
//   (1) dat vi tri ban hat cua gun ngay truoc mat vao (offset -1um) - xem
//       PrimaryGeneratorAction::SetGeometryOffset().
//   (2) doi zNow (toa do tuyet doi tu step->GetPosition()) VE LAI he quy
//       chieu "z tinh tu mat vao" (zRelative = zNow - g_zEntrance_internal)
//       truoc khi so sanh voi 0 (backscatter) hay voi
//       g_targetTotalLength_internal (escape forward) - GIU NGUYEN Y NGHIA
//       VAT LY nhu code goc, chi doi he quy chieu toa do.
static G4double g_zEntrance_internal = 0.0;

// ---------------------------------------------------------------------------
// Bo tich luy ket qua dung chung cho ca calibration va main pass.
// ---------------------------------------------------------------------------
class Accumulator {
public:
    static Accumulator& Instance() {
        static Accumulator inst;
        return inst;
    }
    void ResetForCalibration() {
        fSumZstop_um = 0.0; fNStopped = 0; fNEscapedForward = 0; fNEscapedBackward = 0;
        fMaxZstop_um = 0.0;
        fBinEdep_MeV.clear();
    }
    void ResetForMain(int nBins) {
        fSumZstop_um = 0.0; fNStopped = 0; fNEscapedForward = 0; fNEscapedBackward = 0;
        fMaxZstop_um = 0.0;
        fBinEdep_MeV.assign(nBins, 0.0);
    }
    void AddStopped(double z_um) {
        fSumZstop_um += z_um; ++fNStopped;
        if (z_um > fMaxZstop_um) fMaxZstop_um = z_um;
    }
    // QUAN TRONG: tach rieng 2 nguyen nhan "escaped" vi CHUNG CAN XU LY KHAC
    // NHAU HOAN TOAN:
    //   - Forward (hat di xuyen qua het chieu dai mesh z>totalLen): day la
    //     dau hieu mesh CHUA DU DAI -> retry-loop tang so slab se giai quyet
    //     duoc (giam dan ve 0 khi mesh du dai).
    //   - Backward/backscatter (hat tan xa nguoc ra khoi be mat z<0): day la
    //     HIEU UNG VAT LY THAT (proton nang luong thap bi tan xa hat nhan
    //     goc lon, phan xa nguoc khoi vat lieu nang nhu La) - KHONG lien
    //     quan gi den do dai mesh, tang mesh bao nhieu cung KHONG giam duoc
    //     con so nay. Neu gop chung 2 loai lai, retry-loop se retry vo ich
    //     (dung dung hien tuong ban dang thay: escaped khong giam du mesh
    //     tang tu 15 len 80 slab).
    void AddEscapedForward() { ++fNEscapedForward; }
    void AddEscapedBackward() { ++fNEscapedBackward; }
    // Range LON NHAT quan sat duoc trong so hat da dung han - dai dien cho
    // "duoi" cua phan bo straggling (thay vi chi dung mean nhu truoc), dung
    // de dinh co mesh du dai bao het ca cac hat straggle xa nhat.
    double MaxZstop_um() const { return fMaxZstop_um; }
    // edepMeV: nang luong deposit boi MOT step (cua BAT KY hat nao, so cap
    // hay secondary) roi vao slab co copy number = bin.
    void AddBinEdep(int bin, double edepMeV) {
        if (bin >= 0 && bin < (int)fBinEdep_MeV.size()) fBinEdep_MeV[bin] += edepMeV;
    }

    double MeanRange_um() const {
        return (fNStopped > 0) ? fSumZstop_um / fNStopped : 0.0;
    }
    int NStopped() const { return fNStopped; }
    int NEscapedForward() const { return fNEscapedForward; }
    int NEscapedBackward() const { return fNEscapedBackward; }
    int NEscapedTotal() const { return fNEscapedForward + fNEscapedBackward; }
    double BinEdep_MeV(int bin) const {
        return (bin >= 0 && bin < (int)fBinEdep_MeV.size()) ? fBinEdep_MeV[bin] : 0.0;
    }
    int NBins() const { return (int)fBinEdep_MeV.size(); }
private:
    double fSumZstop_um = 0.0;
    int fNStopped = 0;
    int fNEscapedForward = 0;
    int fNEscapedBackward = 0;
    double fMaxZstop_um = 0.0;
    std::vector<double> fBinEdep_MeV;
};

// ---------------------------------------------------------------------------
// DetectorConstruction: World (vacuum) chua khoi vat lieu, bat dau tu z=0.
// Co 2 che do hinh hoc, dung chung 1 vat lieu da build san (tranh tao lai
// G4Material nhieu lan trong bang vat lieu global):
//   - BuildCalibrationGeometry(): 1 slab dai TARGET_LENGTH_CALIB (uoc luong
//     range tho, khong tally edep theo bin).
//   - BuildSlabGeometry(slabThickness, nSlabs): N slab mong xep lien tiep
//     (dung cho main pass, tally MSP bang edep that su).
// ---------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const std::string& materialName) {
        G4NistManager* nist = G4NistManager::Instance();
        fVacuum = nist->FindOrBuildMaterial("G4_Galactic");
        fTarget = BuildMaterial(materialName, nist);
    }

    // Duoc goi 1 lan boi G4RunManager::Initialize() - tra ve geometry
    // CALIBRATION mac dinh ban dau.
    G4VPhysicalVolume* Construct() override {
        return BuildCalibrationGeometry();
    }

    G4VPhysicalVolume* BuildCalibrationGeometry() {
        return BuildWorldWithSlabs(TARGET_LENGTH_CALIB, 1);
    }

    G4VPhysicalVolume* BuildSlabGeometry(G4double slabThickness, int nSlabs) {
        return BuildWorldWithSlabs(slabThickness, nSlabs);
    }

private:
    // QUAN TRONG - BUG DA SUA (LAN 2): Geant4 BAT BUOC World physical volume
    // phai dat DUNG TAI GOC TOA DO (0,0,0) - G4Navigator::SetWorldVolume() se
    // NEM G4Exception "Volume must be centered on the origin" neu World bi
    // dich tam (day la ly do ban vua gap loi GeomNav0002 sau khi ap dung ban
    // sua lan 1, vi ban sua do dich tam World sang worldZCenter != 0, VI PHAM
    // rang buoc cung nay cua Geant4).
    //
    // Cach sua DUNG: GIU World tai goc toa do (0,0,0) nhu code goc, nhung
    // DICH CA CHONG SLAB sang de no nam GIUA World doi xung qua goc - tuc la
    // slab dau tien bat dau tu z=-totalLen/2 (thay vi z=0 nhu truoc) va slab
    // cuoi ket thuc tai z=+totalLen/2. World tinh kich thuoc doi xung binh
    // thuong (totalLen + margin) se LUON du bao trum ca chong slab bat ke
    // totalLen lon bao nhieu, vi ca 2 deu doi xung qua cung 1 goc.
    //
    // He qua: TOAN BO logic khac (vi tri ban hat cua gun, cac phep so sanh
    // zNow trong SteppingAction de tinh Range/escape) truoc gio dung he toa
    // do "z=0 la mat vao slab dau tien" - GIO PHAI DUOC DICH THEO CUNG 1
    // OFFSET (g_zOffset_internal = -totalLen/2) de giu NGUYEN Y NGHIA VAT LY
    // "z tinh tu mat vao" nhu truoc, chi khac he quy chieu toa do TUYET DOI
    // cua Geant4. Xem PrimaryGeneratorAction::SetGeometryOffset() va
    // SteppingAction (dung g_zOffset_internal) o duoi.
    G4VPhysicalVolume* BuildWorldWithSlabs(G4double slabThickness, int nSlabs) {
        G4double totalLen = slabThickness * nSlabs;
        G4double margin = 1.0 * mm;
        G4double worldZ = totalLen + 2.0 * margin;
        G4double worldXY = TARGET_XY + 1.0*mm;

        auto worldSolid = new G4Box("World", worldXY/2, worldXY/2, worldZ/2);
        auto worldLV = new G4LogicalVolume(worldSolid, fVacuum, "World");
        // World PHAI dat tai goc toa do (0,0,0) - bat buoc cua Geant4.
        auto worldPV = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                          worldLV, "World", nullptr, false, 0);

        auto slabSolid = new G4Box("Slab", TARGET_XY/2, TARGET_XY/2, slabThickness/2);
        auto slabLV = new G4LogicalVolume(slabSolid, fTarget, "SlabLV");
        G4double zStart = -totalLen / 2.0; // mat vao slab dau tien, DICH TAM
        for (int i = 0; i < nSlabs; ++i) {
            G4double zCenter = zStart + i * slabThickness + slabThickness/2;
            // copyNo = i  ->  dung lam "bin index" khi tally edep trong
            // SteppingAction (touchable->GetCopyNumber()).
            new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCenter), slabLV,
                               "SlabPV", worldLV, false, i);
        }
        return worldPV;
    }

    static G4Material* BuildMaterial(const std::string& name, G4NistManager* nist) {
        double density = DENSITY_G_CM3.at(name);
        const auto& comp = COMPOSITION.at(name);
        auto material = new G4Material(name, density * g/cm3, (G4int)comp.size());
        for (const auto& ef : comp) {
            G4Element* elem = nist->FindOrBuildElement(ef.symbol);
            material->AddElement(elem, ef.wt);
        }
        return material;
    }

    G4Material* fVacuum = nullptr;
    G4Material* fTarget = nullptr;
};

// ---------------------------------------------------------------------------
// PhysicsList - DONG BO voi range_calc.cc: dung option3 (thuc nghiem cho
// thay khop SRIM/TRIM tot hon option4 o vung nang luong thap cho proton/
// alpha trong bo vat lieu nay), cut = 1 nm (thay vi 1 um truoc do).
//
// LUU Y QUAN TRONG: day la 2 thay doi CUNG LUC (physics list VA cut) so voi
// ban truoc. Neu can biet chinh xac bien nao la nguyen nhan chinh khien ket
// qua khop SRIM/TRIM tot hon, nen test rieng tung bien (vd option4+1nm va
// option3+1um) truoc khi coi day la cau hinh "chuan" cuoi cung. O day chi
// dong bo cho GIONG HET range_calc.cc de 2 file so sanh duoc apples-to-apples
// (vi sim_msp_range.cc cung tu uoc luong range trong buoc calibration, dung
// khac physics se lam ca mesh sizing lan MSP bi lech theo dai voi
// range_calc.cc, du logic tally hoan toan dung).
// ---------------------------------------------------------------------------
class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList() {
        RegisterPhysics(new G4DecayPhysics());
        RegisterPhysics(new G4EmStandardPhysics_option4());
        SetDefaultCutValue(1.0 * nm);
    }
};

// ---------------------------------------------------------------------------
// PrimaryGeneratorAction - particle VA nang luong deu doi dong duoc.
// ---------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const std::string& particleName) {
        fGun = new G4ParticleGun(1);
        G4ParticleDefinition* particle =
            G4ParticleTable::GetParticleTable()->FindParticle(particleName);
        if (!particle) {
            G4cerr << "Khong tim thay particle: " << particleName << G4endl;
            std::exit(1);
        }
        fGun->SetParticleDefinition(particle);
        fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        fGun->SetParticleEnergy(1.0 * MeV);
        SetGeometryOffset(0.0); // mac dinh, se duoc goi lai dung truoc moi BeamOn
    }
    ~PrimaryGeneratorAction() override { delete fGun; }

    void SetEnergyMeV(G4double eMeV) { fGun->SetParticleEnergy(eMeV * MeV); }

    // PHAI goi truoc moi BeamOn(), NGAY SAU khi DefineWorldVolume() voi mot
    // world moi: zEntrance_internal la toa do TUYET DOI (frame Geant4, world
    // dat tai goc) cua mat vao slab dau tien (= -totalLen/2 vi chong slab gio
    // doi xung qua goc - xem BuildWorldWithSlabs). Gun duoc dat ngay TRUOC
    // mat vao 1 khoang nho (1um) theo huong ban (+z), GIONG HET y nghia vat
    // ly cu (truoc kia mat vao la z=0, gun dat tai z=-1um).
    void SetGeometryOffset(G4double zEntrance_internal) {
        fGun->SetParticlePosition(G4ThreeVector(0., 0., zEntrance_internal - 1.0*um));
    }

    void GeneratePrimaries(G4Event* anEvent) override {
        fGun->GeneratePrimaryVertex(anEvent);
    }
private:
    G4ParticleGun* fGun;
};

// ---------------------------------------------------------------------------
// SteppingAction: 2 nhiem vu HOAN TOAN DOC LAP nhau:
//   (A) Tally edep theo slab (MOI hat, chi hoat dong khi g_binningActive) -
//       day la co so tinh MSP(entrance) = BinEdep(0)/nEvents/slabThickness/density.
//   (B) Theo doi RIENG track so cap (ParentID==0, TrackID==1) de xac dinh z
//       luc dung han (Range) va phat hien thoat khoi vat lieu/tan xa nguoc -
//       khong lien quan gi toi (A), giu logic don gian, on dinh nhu ban goc.
// ---------------------------------------------------------------------------
class SteppingAction : public G4UserSteppingAction {
public:
    void UserSteppingAction(const G4Step* step) override {
        // --- (A) Edep-tally theo slab, cho MOI hat -------------------------
        if (g_binningActive) {
            G4double edep = step->GetTotalEnergyDeposit();
            if (edep > 0.0) {
                const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
                if (touch && touch->GetVolume() &&
                    touch->GetVolume()->GetName() == "SlabPV") {
                    G4int binIdx = touch->GetCopyNumber();
                    Accumulator::Instance().AddBinEdep(binIdx, edep / MeV);
                }
            }
        }

        // --- (B) Rieng track so cap: Range + phat hien thoat/tan xa nguoc --
        G4Track* track = step->GetTrack();
        if (track->GetParentID() != 0) return;
        if (track->GetTrackID() != 1) return;

        G4double zNow = step->GetPostStepPoint()->GetPosition().z()
                        - g_zEntrance_internal; // doi ve he "z tinh tu mat vao"
        G4double eKin = track->GetKineticEnergy();

        // QUAN TRONG: kiem tra BACKSCATTER (zNow < 0) PHAI DUNG TRUOC kiem
        // tra "da dung han" (eKin<=1eV). Ly do: mot hat co the tan xa nguoc
        // ra sau be mat VA dong thoi mat not chut dong nang cuoi cung NGAY
        // TRONG CUNG 1 STEP do (vi tan xa manh thuong di kem mat nang luong
        // lon). Neu kiem tra eKin truoc, hat nay se bi tinh nham thanh "dung
        // han tai z AM" -> cong z am vao Range -> Range trung binh co the ra
        // AM (bug da gap: Range=-0.083 um). Ve mat vat ly, mot khi track da
        // di nguoc qua be mat, no khong con "range tien" co nghia nua -
        // PHAI duoc tinh la backscatter bat ke luc do con dong nang hay
        // khong.
        static const G4double kBackscatterMargin = -0.01 * um;
        if (zNow < kBackscatterMargin) {
            Accumulator::Instance().AddEscapedBackward();
            track->SetTrackStatus(fStopAndKill);
            return;
        }

        if (eKin <= 1.0 * CLHEP::eV) {
            Accumulator::Instance().AddStopped(zNow / um);
            track->SetTrackStatus(fStopAndKill);
            return;
        }

        if (zNow > g_targetTotalLength_internal) {
            Accumulator::Instance().AddEscapedForward();
            track->SetTrackStatus(fStopAndKill);
        }
    }
};

// ---------------------------------------------------------------------------
// Doc danh sach nang luong keV tu file text
// ---------------------------------------------------------------------------
static std::vector<double> ReadEnergyListKeV(const std::string& path) {
    std::vector<double> energies;
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Khong mo duoc file nang luong: " << path << std::endl;
        return energies;
    }
    std::string line;
    while (std::getline(f, line)) {
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        if (line[start] == '#') continue;
        std::istringstream iss(line);
        double val;
        if (iss >> val) energies.push_back(val);
    }
    return energies;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Cach dung: ./sim_msp_range <particle: proton|alpha> "
                     "<material s1..s6> [energies_kev.txt] [outputPrefix] "
                     "[nEvents] [calibEvents]\n";
        return 1;
    }
    std::string particleName = argv[1];
    std::string material     = argv[2];
    std::string energyFile   = (argc > 3) ? argv[3] : "energies_kev.txt";
    std::string outPrefix    = (argc > 4) ? argv[4] : (particleName + "_" + material);
    int nEvents               = (argc > 5) ? std::stoi(argv[5]) : 3000;
    // calibEvents: GIU LAI cho tuong thich dong lenh cu (vd script/wrapper
    // hien co van truyen tham so thu 6 nay), nhung KHONG con duoc dung nua -
    // Range-pass moi dung thang nEvents day du tren khoi lien (giong het
    // range_calc.cc) thay vi ban rieng calibEvents cho 1 pha calibration mau
    // nho nhu truoc. (void) de tranh canh bao "unused variable".
    int calibEvents           = (argc > 6) ? std::stoi(argv[6]) : 200;
    (void)calibEvents;

    if (particleName != "proton" && particleName != "alpha") {
        std::cerr << "Particle khong hop le: " << particleName
                  << " (phai la proton hoac alpha)\n";
        return 1;
    }
    if (COMPOSITION.find(material) == COMPOSITION.end()) {
        std::cerr << "Vat lieu khong hop le: " << material << " (phai la s1..s6)\n";
        return 1;
    }

    std::vector<double> energiesKeV = ReadEnergyListKeV(energyFile);
    if (energiesKeV.empty()) {
        std::cerr << "Danh sach nang luong rong. Kiem tra file: " << energyFile << std::endl;
        return 1;
    }
    CLHEP::HepRandom::setTheSeed(static_cast<long>(std::time(nullptr)));
    auto runManager = new G4RunManager();
    auto detector = new DetectorConstruction(material);
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList());

    auto primaryGen = new PrimaryGeneratorAction(particleName);
    runManager->SetUserAction(primaryGen);
    auto steppingAction = new SteppingAction();
    runManager->SetUserAction(steppingAction);

    runManager->Initialize();   // xay dung geometry CALIBRATION lan dau (Construct())

    double density = DENSITY_G_CM3.at(material);

    std::ofstream csvMsp(outPrefix + "_msp.csv");
    csvMsp << "Energy_keV,MSP_MeV_cm2_per_g\n";
    std::ofstream csvRange(outPrefix + "_range.csv");
    csvRange << "Energy_keV,Range_um,NStopped,NEscapedForward,NEscapedBackward\n";

    for (double eKeV : energiesKeV) {
        double eMeV = eKeV / 1000.0;

        // === RANGE-PASS (rieng biet, TACH HOAN TOAN khoi MSP main pass):
        //     dung DUNG geometry 1-khoi-lien (giong het range_calc.cc) VA
        //     DUNG nEvents day du (khong con dung calibEvents mau nho nua) -
        //     day la nguon "su that" duy nhat cho Range bao cao ra CSV/
        //     console (rangeFinalDedicated_um). MSP main pass (mesh slab)
        //     phia duoi CHI dung de tinh MSP(entrance), khong con dung de
        //     bao cao Range nua - no van can 1 uoc luong range tho de tu
        //     dinh kich thuoc mesh, va gio lay luon tu chinh Range-pass nay
        //     (rangeEst_um/maxRangeEst_um) thay vi tu 1 pha calibration mau
        //     nho rieng - vua chinh xac hon (nEvents day du) vua khong con
        //     ban trung 2 lan cho cung 1 muc dich. ==========================
        primaryGen->SetEnergyMeV(eMeV);

        G4VPhysicalVolume* rangeWorld = detector->BuildCalibrationGeometry();
        runManager->DefineWorldVolume(rangeWorld, true);
        g_binningActive = false;
        g_targetTotalLength_internal = TARGET_LENGTH_CALIB;
        // Khoi calibration/Range-pass la 1 slab duy nhat (nSlabs=1), nen mat
        // vao cua no (theo BuildWorldWithSlabs moi) nam tai -TARGET_LENGTH_
        // CALIB/2 trong he toa do tuyet doi cua Geant4 (World o goc).
        g_zEntrance_internal = -TARGET_LENGTH_CALIB / 2.0;
        primaryGen->SetGeometryOffset(g_zEntrance_internal);

        Accumulator::Instance().ResetForCalibration();
        runManager->BeamOn(nEvents);
        Accumulator* rangeAcc = &Accumulator::Instance();

        double rangeFinalDedicated_um    = rangeAcc->MeanRange_um();
        int    rangePassNStopped         = rangeAcc->NStopped();
        int    rangePassNEscapedForward  = rangeAcc->NEscapedForward();
        int    rangePassNEscapedBackward = rangeAcc->NEscapedBackward();

        double rangeEst_um = rangeFinalDedicated_um;
        double maxRangeEst_um = rangeAcc->MaxZstop_um();
        if (rangeEst_um <= 0.0) rangeEst_um = MIN_SLAB_UM * 100;
        if (maxRangeEst_um <= 0.0) maxRangeEst_um = rangeEst_um;

        // === Tinh tham so mesh slab cho pha MAIN ==========================
        // Dung ca (a) mean range x he so co dinh VA (b) max range QUAN SAT
        // DUOC trong calibration x safety margin - lay gia tri LON HON trong
        // 2 cai, de mesh tu dong du dai hon khi straggling manh (nang luong
        // thap / vat lieu nang), thay vi chi dua vao 1 he so co dinh 1.6x
        // mean (co the khong du khi phan bo range co duoi dai).
        double totalLen_um = std::max(rangeEst_um * AUTO_TOTAL_LENGTH_FACTOR,
                                       maxRangeEst_um * MAX_STOP_SAFETY_MARGIN);

        // QUAN TRONG: be day slab (do phan giai vung entrance) duoc CHOT MOT
        // LAN DUY NHAT o day va GIU CO DINH xuyen suot moi lan retry ben
        // duoi. Neu de slabThickness = totalLen/nSlabs voi nSlabs co dinh,
        // moi lan retry "phong to" totalLen se lam slabThickness PHINH TO
        // THEO -> do phan giai entrance cang retry cang te (dung bug vua
        // gap). Cach dung: retry chi duoc phep THEM slab moi o phia sau (giu
        // slabThickness khong doi), khong duoc "gian" slab hien co.
        double slabThickness_um = std::max(MIN_SLAB_UM,
                                            totalLen_um / AUTO_NSLABS_DEFAULT);
        int nSlabs = std::max(1, (int)std::lround(totalLen_um / slabThickness_um));
        totalLen_um = slabThickness_um * nSlabs;

        // Tran an toan cho so slab (tranh mesh phinh qua lon ve so luong khi
        // retry nhieu lan, gay cham may / ton bo nho).
        static const int MAX_NSLABS_HARD_CAP = 500;

        // === Pha 2: MAIN PASS - mesh N slab, tally edep THAT SU cho MOI hat,
        //     Range van tinh tu track so cap nhu cu (khong doi). Neu con hat
        //     thoat mesh (escaped) VUOT NGUONG CHAP NHAN DUOC, TU DONG THEM
        //     SLAB MOI o phia sau (KHONG doi slabThickness) dua tren chinh du
        //     lieu 3000 hat nay (khong con phu thuoc uoc luong tu 200 hat
        //     calibration nua) roi chay lai, toi khi du it escaped hoac het
        //     so lan thu cho phep. Nguong chap nhan (thay vi doi hoi ==0)
        //     vi escaped con lai la nhieu thong ke tu nhien (extreme-value
        //     tail), khong the/khong can dua ve dung 0 tuyet doi. ===========
        static const double ESCAPE_FRACTION_TOLERANCE = 0.002; // 0.2%
        Accumulator* acc = nullptr;
        int retryCount = 0;
        for (;;) {
            G4VPhysicalVolume* mainWorld =
                detector->BuildSlabGeometry(slabThickness_um * um, nSlabs);
            runManager->DefineWorldVolume(mainWorld, true);
            g_binningActive = true;
            g_targetTotalLength_internal = totalLen_um * um;
            // Mesh N-slab gio doi xung qua goc (BuildWorldWithSlabs moi), nen
            // mat vao slab dau tien nam tai -totalLen/2 trong he toa do tuyet
            // doi. PHAI cap nhat lai moi lan mesh doi kich thuoc (retry).
            g_zEntrance_internal = -(totalLen_um * um) / 2.0;
            primaryGen->SetGeometryOffset(g_zEntrance_internal);

            Accumulator::Instance().ResetForMain(nSlabs);
            runManager->BeamOn(nEvents);
            acc = &Accumulator::Instance();

            double escapedForwardFrac = (double)acc->NEscapedForward() / (double)nEvents;
            if (escapedForwardFrac <= ESCAPE_FRACTION_TOLERANCE ||
                retryCount >= MAX_MESH_RETRIES ||
                nSlabs >= MAX_NSLABS_HARD_CAP) {
                break;
            }

            // THEM slab moi o phia sau dua tren max Z-stop THAT SU quan sat
            // duoc trong chinh nEvents hat vua chay (tail cua mau lon, dang
            // tin cay hon nhieu so voi uoc luong tu calibEvents). KHONG doi
            // slabThickness - chi tang nSlabs de keo dai mesh. LUU Y: CHI
            // retry dua tren FORWARD escape - backward (backscatter) la
            // hieu ung vat ly that (proton nang luong thap tan xa hat nhan
            // goc lon, phan xa nguoc khoi vat lieu nang), KHONG lien quan gi
            // toi do dai mesh, khong the/khong nen "fix" bang cach retry.
            double observedMax_um = acc->MaxZstop_um();
            double newTotalLen_um = std::max(totalLen_um * MESH_GROWTH_FACTOR,
                                              observedMax_um * MESH_GROWTH_FACTOR);
            int newNSlabs = std::max(nSlabs + 1,
                                      (int)std::lround(newTotalLen_um / slabThickness_um));
            newNSlabs = std::min(newNSlabs, MAX_NSLABS_HARD_CAP);
            std::cerr << "  [TU DONG THEM SLAB] " << acc->NEscapedForward() << "/"
                      << nEvents << " hat thoat mesh PHIA TRUOC (totalLen="
                      << totalLen_um << " um, " << nSlabs << " slab) -> tang len "
                      << newNSlabs << " slab (totalLen=" << (slabThickness_um*newNSlabs)
                      << " um, slabThickness KHONG DOI=" << slabThickness_um
                      << " um) (lan thu " << (retryCount+1) << "/"
                      << MAX_MESH_RETRIES << ")\n";
            nSlabs = newNSlabs;
            totalLen_um = slabThickness_um * nSlabs;
            ++retryCount;
        }

        // Range bao cao ra CSV/console GIO LAY TU RANGE-PASS RIENG (khoi lien,
        // nEvents day du) - KHONG con lay tu acc cua MSP main pass (mesh
        // slab) nua. acc (mesh) chi con dung cho MSP(entrance) ben duoi.
        double slabThickness_cm = slabThickness_um * 1.0e-4;
        double entranceEdep_MeV = acc->BinEdep_MeV(0);
        double msp = (entranceEdep_MeV / nEvents) / slabThickness_cm / density;

        csvMsp << eKeV << "," << msp << "\n";
        csvRange << eKeV << "," << rangeFinalDedicated_um << "," << rangePassNStopped
                  << "," << rangePassNEscapedForward << "," << rangePassNEscapedBackward << "\n";

        std::cout << "[" << particleName << "/" << material << "] E=" << eKeV
                  << " keV -> MSP=" << msp << " MeV.cm2/g, Range=" << rangeFinalDedicated_um
                  << " um (RangePass: khoi lien " << TARGET_LENGTH_CALIB/mm << " mm, stopped="
                  << rangePassNStopped << ", escaped_forward=" << rangePassNEscapedForward
                  << ", escaped_backward(backscatter)=" << rangePassNEscapedBackward
                  << ", /" << nEvents
                  << "; MSP mesh: slab=" << slabThickness_um << " um x " << nSlabs
                  << " = " << totalLen_um << " um, escaped_forward(mesh)="
                  << acc->NEscapedForward() << "/" << nEvents << ")"
                  << std::endl;

        double rangePassEscapedForwardFrac = (double)rangePassNEscapedForward / (double)nEvents;
        if (rangePassEscapedForwardFrac > ESCAPE_FRACTION_TOLERANCE) {
            std::cerr << "  [CANH BAO] Range-pass: " << rangePassNEscapedForward
                      << "/" << nEvents << " (" << (rangePassEscapedForwardFrac*100)
                      << "%) hat thoat khoi khoi lien TARGET_LENGTH_CALIB="
                      << (TARGET_LENGTH_CALIB/mm) << " mm truoc khi dung han -> "
                         "Range o nang luong nay co the bi hut nhe. Tang "
                         "TARGET_LENGTH_CALIB neu can chinh xac Range cao hon "
                         "(luu y: file nay dung chung 1 khoi TARGET_LENGTH_CALIB "
                         "cho moi muc nang luong, khong tu phong to nhu mesh).\n";
        }
        double finalEscapedForwardFrac = (double)acc->NEscapedForward() / (double)nEvents;
        if (finalEscapedForwardFrac > ESCAPE_FRACTION_TOLERANCE) {
            std::cerr << "  [CANH BAO] MSP mesh: Sau " << MAX_MESH_RETRIES
                      << " lan thu them slab, van con " << acc->NEscapedForward()
                      << "/" << nEvents << " (" << (finalEscapedForwardFrac*100)
                      << "%) hat thoat mesh PHIA TRUOC -> khong anh huong Range "
                         "(da tach rieng ve Range-pass), nhung co the anh huong "
                         "MSP profile o cac slab sau. Tang MAX_MESH_RETRIES, "
                         "MESH_GROWTH_FACTOR, hoac MAX_NSLABS_HARD_CAP neu can.\n";
        }
        double backscatterFrac = (double)rangePassNEscapedBackward / (double)nEvents;
        if (backscatterFrac > 0.002) {
            std::cerr << "  [GHI CHU] Range-pass: " << rangePassNEscapedBackward << "/"
                      << nEvents << " (" << (backscatterFrac*100)
                      << "%) hat tan xa NGUOC ra khoi be mat (backscatter) - "
                         "day la HIEU UNG VAT LY THAT o nang luong thap trong "
                         "vat lieu nang (La), KHONG phai loi mesh, KHONG the "
                         "sua bang cach tang do dai mesh. Cac hat nay khong co "
                         "'range tien' co nghia nen duoc loai khoi Range trung "
                         "binh mot cach hop ly.\n";
        }
        if (rangeFinalDedicated_um < slabThickness_um * 2.0) {
            std::cerr << "  [CANH BAO] Range (" << rangeFinalDedicated_um
                      << " um) qua gan be day 1 slab MSP-mesh (" << slabThickness_um
                      << " um) -> slab dau tien gan nhu chua het toan bo range, "
                         "MSP(entrance) o diem nay KHONG CON dai dien cho vung "
                         "entrance nua. Giam MIN_SLAB_UM hoac tang nEvents/"
                         "AUTO_NSLABS_DEFAULT cho vung nang luong nay.\n";
        }

        // === Xuat them profile MSP theo do sau (kiem tra hinh dang Bragg,
        //     dung de validate hoi tu / phat hien bin qua tho hay qua mong) ==
        std::ostringstream profName;
        profName << outPrefix << "_profile_" << (long)eKeV << "keV.csv";
        std::ofstream profFile(profName.str());
        profFile << "bin,z_center_um,MSP_MeV_cm2_per_g\n";
        for (int b = 0; b < nSlabs; ++b) {
            double zc = (b + 0.5) * slabThickness_um;
            double m = (acc->BinEdep_MeV(b) / nEvents) / slabThickness_cm / density;
            profFile << b << "," << zc << "," << m << "\n";
        }
    }
    csvMsp.close();
    csvRange.close();

    std::cout << "\nDa luu: " << outPrefix << "_msp.csv, "
              << outPrefix << "_range.csv, va cac file "
              << outPrefix << "_profile_<E>keV.csv (profile MSP theo do sau, "
                 "dung de kiem tra hoi tu / hinh dang duong cong Bragg)"
              << std::endl;

    delete runManager;
    return 0;
}
