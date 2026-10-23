// ---------------------------------------------------------------------------
// range_calc : tinh PROJECTED RANGE THAT (khong phai CSDA) cho proton HOAC
// alpha trong MOT vat lieu (s1..s6) tren dai nang luong tuy chinh, bang
// tracking Monte Carlo thuc su:
//   - Ban N hat (proton/alpha) vao 1 khoi vat lieu don (bat dau tu z=0)
//   - Ghi lai toa do z LUC HAT DUNG HAN (kinetic energy ~ 0)
//   - Trung binh z_stop qua N hat -> Projected Range
//
// Chay RIENG cho tung to hop (particle, vat lieu) - dung run_all_ranges.bat
// de tu dong chay het 2 hat x 6 vat lieu = 12 file CSV, roi dung
// merge_to_excel.py de gop lai thanh 1 file .xlsx voi 2 sheet (Proton/Alpha).
//
// Bien dich: da duoc them vao CMakeLists.txt cung project, chi can
//   cd build && cmake --build . --config Release
//
// Chay:
//   range_calc.exe <particle proton|alpha> <material s1..s6> \
//                   [energies_kev.txt] [output.csv] [nEvents] [seed]
//   vi du: range_calc.exe alpha s3 ../energies_kev.txt range_alpha_s3.csv 3000 111
// ---------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmParameters.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Thanh phan + mat do 6 vat lieu (giong het main.cc)
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

// Chieu dai khoi vat lieu (z) - PHAI DU LON de chua het range o nang luong
// cao nhat trong danh sach. Alpha o cung nang luong (keV) co range NGAN HON
// proton nhieu (do dE/dx cao hon ~4x do dien tich Z=2 va khoi luong lon hon),
// nen 2mm van thua du cho ca 2 loai hat trong dai nang luong thong thuong
// (~15 keV - 15 MeV). Neu chay alpha o nang luong rat cao (>50 MeV) can tang
// gia tri nay len va build lai.
static const G4double TARGET_LENGTH = 2.0 * mm;
static const G4double TARGET_XY     = 2.0 * mm;   // be rong ngang (x,y)

// ---------------------------------------------------------------------------
// Bo tich luy ket qua (don gian vi chay single-thread)
// ---------------------------------------------------------------------------
class RangeAccumulator {
public:
    static RangeAccumulator& Instance() {
        static RangeAccumulator inst;
        return inst;
    }
    void Reset() { fSumZ_um = 0.0; fNStopped = 0; fNEscapedForward = 0; fNEscapedBackward = 0; }
    void AddStopped(G4double z_um) { fSumZ_um += z_um; ++fNStopped; }
    void AddEscapedForward() { ++fNEscapedForward; }
    void AddEscapedBackward() { ++fNEscapedBackward; }
    G4double MeanRange_um() const {
        return (fNStopped > 0) ? fSumZ_um / fNStopped : 0.0;
    }
    int NStopped() const { return fNStopped; }
    int NEscapedForward() const { return fNEscapedForward; }
    int NEscapedBackward() const { return fNEscapedBackward; }
private:
    double fSumZ_um = 0.0;
    int fNStopped = 0;
    int fNEscapedForward = 0;
    int fNEscapedBackward = 0;
};

// ---------------------------------------------------------------------------
// DetectorConstruction: World (vacuum) chua 1 KHOI VAT LIEU DUY NHAT
// ---------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const std::string& materialName)
        : fMaterialName(materialName) {}

    G4VPhysicalVolume* Construct() override {
        G4NistManager* nist = G4NistManager::Instance();
        G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
        G4Material* target = BuildMaterial(fMaterialName, nist);

        G4double worldZ = TARGET_LENGTH + 1.0*mm;
        G4double worldXY = TARGET_XY + 1.0*mm;

        auto worldSolid = new G4Box("World", worldXY/2, worldXY/2, worldZ/2);
        auto worldLV = new G4LogicalVolume(worldSolid, vacuum, "World");
        auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);

        auto targetSolid = new G4Box("Target", TARGET_XY/2, TARGET_XY/2, TARGET_LENGTH/2);
        auto targetLV = new G4LogicalVolume(targetSolid, target, "TargetLV");
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, TARGET_LENGTH/2),
                           targetLV, "TargetPV", worldLV, false, 0);

        return worldPV;
    }

private:
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

    std::string fMaterialName;
};

// ---------------------------------------------------------------------------
// PhysicsList - dung option3 (khop SRIM/TRIM tot hon option4 o vung nang
// luong thap, da kiem chung qua doi chieu voi TRIM va PHITS).
// ---------------------------------------------------------------------------
class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList() {
        RegisterPhysics(new G4DecayPhysics());
        RegisterPhysics(new G4EmStandardPhysics_option3());
        SetDefaultCutValue(1.0 * nm);
    }
};

// ---------------------------------------------------------------------------
// PrimaryGeneratorAction - nhan ten hat (proton/alpha) qua constructor.
// ---------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const std::string& particleName) {
        fGun = new G4ParticleGun(1);
        G4ParticleDefinition* particle =
            G4ParticleTable::GetParticleTable()->FindParticle(particleName);
        if (!particle) {
            G4cerr << "LOI: khong tim thay dinh nghia hat '" << particleName
                   << "' trong G4ParticleTable." << G4endl;
            std::exit(1);
        }
        fGun->SetParticleDefinition(particle);
        fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        // bat dau NGAY TRUOC be mat khoi vat lieu (z=0), di chuyen vao +z
        fGun->SetParticlePosition(G4ThreeVector(0., 0., -1.0*um));
        fGun->SetParticleEnergy(1.0 * MeV);
    }
    ~PrimaryGeneratorAction() override { delete fGun; }

    void SetEnergyMeV(G4double eMeV) { fGun->SetParticleEnergy(eMeV * MeV); }

    void GeneratePrimaries(G4Event* anEvent) override {
        fGun->GeneratePrimaryVertex(anEvent);
    }
private:
    G4ParticleGun* fGun;
};

// ---------------------------------------------------------------------------
// SteppingAction: theo doi hat so cap (primary), ghi lai z luc dung han.
// ---------------------------------------------------------------------------
class SteppingAction : public G4UserSteppingAction {
public:
    void UserSteppingAction(const G4Step* step) override {
        G4Track* track = step->GetTrack();
        if (track->GetParentID() != 0) return;   // chi quan tam hat so cap
        if (track->GetTrackID() != 1) return;

        G4double eKin = track->GetKineticEnergy();
        G4StepPoint* post = step->GetPostStepPoint();
        G4double zNow = post->GetPosition().z();

        // Backscatter PHAI duoc kiem tra TRUOC "da dung han" - xem giai
        // thich chi tiet trong ban goc cua file nay.
        static const G4double kBackscatterMargin = -0.01 * um;
        if (zNow < kBackscatterMargin) {
            RangeAccumulator::Instance().AddEscapedBackward();
            track->SetTrackStatus(fStopAndKill);
            return;
        }

        if (eKin <= 1.0 * CLHEP::eV) {
            RangeAccumulator::Instance().AddStopped(zNow / um);
            track->SetTrackStatus(fStopAndKill);
            return;
        }

        if (zNow > TARGET_LENGTH) {
            RangeAccumulator::Instance().AddEscapedForward();
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
        std::cerr << "Cach dung: range_calc.exe <particle proton|alpha> "
                     "<material s1..s6> [energies_kev.txt] [output.csv] "
                     "[nEvents] [seed]\n"
                     "Vi du: range_calc.exe alpha s3 ../energies_kev.txt "
                     "range_alpha_s3.csv 3000 111\n";
        return 1;
    }

    // Ep dang ky dinh nghia hat NGAY TU DAU, truoc khi PhysicsList/RunManager
    // khoi tao - tranh truong hop FindParticle("alpha") tra ve nullptr do
    // G4Alpha chua duoc dang ky (proton thuong da san co, nhung alpha thi
    // phu thuoc thu tu khoi tao physics list nen can force o day cho chac).
    G4Proton::ProtonDefinition();
    G4Alpha::AlphaDefinition();

    std::string particleName = argv[1];
    std::string material     = argv[2];
    std::string energyFile = (argc > 3) ? argv[3] : "energies_kev.txt";
    std::string outFile    = (argc > 4) ? argv[4]
                              : ("range_" + particleName + "_" + material + ".csv");
    int nEvents             = (argc > 5) ? std::stoi(argv[5]) : 3000;
    long seed                = (argc > 6) ? std::stol(argv[6]) : 12345;

    if (particleName != "proton" && particleName != "alpha") {
        std::cerr << "Particle khong hop le: " << particleName
                  << " (phai la 'proton' hoac 'alpha')\n";
        return 1;
    }
    if (COMPOSITION.find(material) == COMPOSITION.end()) {
        std::cerr << "Vat lieu khong hop le: " << material
                  << " (phai la s1..s6)\n";
        return 1;
    }

    CLHEP::HepRandom::setTheSeed(seed);

    std::vector<double> energiesKeV = ReadEnergyListKeV(energyFile);
    if (energiesKeV.empty()) {
        std::cerr << "Danh sach nang luong rong. Kiem tra file: "
                  << energyFile << std::endl;
        return 1;
    }

    auto runManager = new G4RunManager();
    auto detector = new DetectorConstruction(material);
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList());

    auto primaryGen = new PrimaryGeneratorAction(particleName);
    runManager->SetUserAction(primaryGen);
    runManager->SetUserAction(new SteppingAction());

    runManager->Initialize();

    std::ofstream csv(outFile);
    csv << "Particle,Material,Energy_keV,Range_um,NStopped,NEscapedForward,NEscapedBackward\n";

    for (double eKeV : energiesKeV) {
        RangeAccumulator::Instance().Reset();
        primaryGen->SetEnergyMeV(eKeV / 1000.0);
        runManager->BeamOn(nEvents);

        auto& acc = RangeAccumulator::Instance();
        csv << particleName << "," << material << "," << eKeV << ","
            << acc.MeanRange_um() << "," << acc.NStopped() << ","
            << acc.NEscapedForward() << "," << acc.NEscapedBackward() << "\n";

        std::cout << "[" << particleName << "/" << material << "] E=" << eKeV
                  << " keV -> Range=" << acc.MeanRange_um() << " um "
                  << "(stopped=" << acc.NStopped()
                  << ", escaped_forward=" << acc.NEscapedForward()
                  << ", escaped_backward(backscatter)=" << acc.NEscapedBackward()
                  << ")" << std::endl;

        if (acc.NEscapedForward() > 0) {
            std::cerr << "  [CANH BAO] " << acc.NEscapedForward()
                      << " hat thoat khoi khoi vat lieu PHIA TRUOC ma chua "
                         "dung han -> TARGET_LENGTH dang qua nho cho nang "
                         "luong nay. Tang TARGET_LENGTH trong code va build lai.\n";
        }
        double backFrac = (double)acc.NEscapedBackward() / (double)nEvents;
        if (backFrac > 0.002) {
            std::cerr << "  [GHI CHU] " << acc.NEscapedBackward() << "/"
                      << nEvents << " (" << (backFrac*100)
                      << "%) hat tan xa NGUOC ra khoi be mat (backscatter) - "
                         "hieu ung vat ly that o nang luong thap.\n";
        }
        int totalAccounted = acc.NStopped() + acc.NEscapedForward() + acc.NEscapedBackward();
        if (totalAccounted != nEvents) {
            std::cerr << "  [CANH BAO NGHIEM TRONG] Chi hach toan duoc "
                      << totalAccounted << "/" << nEvents << " hat -> co "
                      << (nEvents - totalAccounted) << " hat BI MAT AM THAM.\n";
        }
    }
    csv.close();

    std::cout << "\nDa luu ket qua vao '" << outFile << "'" << std::endl;

    delete runManager;
    return 0;
}
