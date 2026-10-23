// ---------------------------------------------------------------------------
// msp_calc : tinh Mass Stopping Power (MeV.cm2/g) cho PROTON trong 6 vat lieu
// s1..s6 (thuy tinh B-O-Si-Ca-La), tren mot dai nang luong TUY CHINH doc tu
// file text (energies_kev.txt, moi dong 1 gia tri keV).
//
// Dung G4EmCalculator::ComputeTotalDEDX(), khong can chay simulation/tracking
// hat that su - chi can Physics Table duoc build (qua BeamOn(1) toi thieu).
//
// Bien dich:
//   mkdir build && cd build
//   cmake -DGeant4_DIR=<path_to_Geant4Config.cmake> ..
//   make -j$(nproc)
//
// Chay:
//   ./msp_calc [energies_kev.txt] [output.csv]
//   (mac dinh: energies_kev.txt va msp_results.csv trong thu muc hien tai)
// ---------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
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
#include "G4Event.hh"

#include "G4EmCalculator.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// Thanh phan (wt fraction, DUONG - khac quy uoc am cua PHITS/MCNP) va mat do
// cho 6 vat lieu s1..s6 - lay dung tu generate_msp_range_inputs.py, chi doi
// dau am -> duong va ten ZAID -> ky hieu nguyen to.
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

// Thu tu co dinh de xuat CSV luon dung mot layout Energy x s1..s6
static const std::vector<std::string> MATERIAL_ORDER = {"s1","s2","s3","s4","s5","s6"};

// ---------------------------------------------------------------------------
// DetectorConstruction toi gian: World (vacuum) chua 6 khoi vat lieu s1..s6
// dat canh nhau (khong chong lan) chi de dang ky material vao geometry ->
// Geant4 se build physics table (dE/dx) cho tung material nay.
// ---------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        G4NistManager* nist = G4NistManager::Instance();

        G4double worldSize = 20.0 * cm;
        G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

        auto worldSolid = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
        auto worldLV = new G4LogicalVolume(worldSolid, vacuum, "World");
        auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);

        G4double boxSize = 1.0 * cm;
        G4double spacing = 2.0 * cm;
        int i = 0;
        for (const auto& mat : MATERIAL_ORDER) {
            G4Material* material = BuildMaterial(mat, nist);
            auto solid = new G4Box("Box_" + mat, boxSize/2, boxSize/2, boxSize/2);
            auto lv = new G4LogicalVolume(solid, material, "LV_" + mat);
            G4double x = (i - (int)MATERIAL_ORDER.size()/2) * spacing;
            new G4PVPlacement(nullptr, G4ThreeVector(x, 0, 0), lv, "PV_" + mat, worldLV, false, 0);
            ++i;
        }
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
};

// ---------------------------------------------------------------------------
// PhysicsList toi gian: chi can EM (option4 - chinh xac nhat vung nang luong
// thap, dung Livermore/Penelope-based models, phu hop 15 keV - 15 MeV) + Decay
// (de cac particle duoc dang ky day du).
// ---------------------------------------------------------------------------
class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList() {
        RegisterPhysics(new G4DecayPhysics());
        RegisterPhysics(new G4EmStandardPhysics_option4());
    }
};

// ---------------------------------------------------------------------------
// PrimaryGeneratorAction toi thieu - chi de RunManager co du mandatory action,
// khong dung de tracking that (chung ta khong doc ket qua tu day).
// ---------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction() {
        fGun = new G4ParticleGun(1);
        G4ParticleDefinition* proton =
            G4ParticleTable::GetParticleTable()->FindParticle("proton");
        fGun->SetParticleDefinition(proton);
        fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        fGun->SetParticleEnergy(1.0 * MeV);
        fGun->SetParticlePosition(G4ThreeVector(0., 0., -9.0 * cm));
    }
    ~PrimaryGeneratorAction() override { delete fGun; }
    void GeneratePrimaries(G4Event* anEvent) override {
        fGun->GeneratePrimaryVertex(anEvent);
    }
private:
    G4ParticleGun* fGun;
};

// ---------------------------------------------------------------------------
// Doc danh sach nang luong (keV) tu file text, moi dong 1 gia tri.
// Bo qua dong trong / dong bat dau bang '#'.
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
        // bo whitespace
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
    std::string energyFile = (argc > 1) ? argv[1] : "energies_kev.txt";
    std::string outFile    = (argc > 2) ? argv[2] : "msp_results.csv";

    std::vector<double> energiesKeV = ReadEnergyListKeV(energyFile);
    if (energiesKeV.empty()) {
        std::cerr << "Danh sach nang luong rong, dung lai. Kiem tra file: "
                  << energyFile << std::endl;
        return 1;
    }

    // --- Bat co build bang CSDA range TRUOC khi Initialize() ----------------
    // MAC DINH Geant4 KHONG build bang nay (de tiet kiem thoi gian/bo nho),
    // neu quen buoc nay thi GetCSDARange() luon tra ve 0 nhu da gap phai.
    G4EmParameters::Instance()->SetBuildCSDARange(true);

    // --- Khoi tao Geant4 run manager toi gian --------------------------------
    auto runManager = new G4RunManager();
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserAction(new PrimaryGeneratorAction());

    runManager->Initialize();
    runManager->BeamOn(1);   // build physics table (khong dung ket qua event nay)

    // --- Tinh MSP cho tung vat lieu x tung nang luong ------------------------
    G4EmCalculator emCal;

    // CSV MSP (MeV.cm2/g)
    std::ofstream csvMsp(outFile);
    csvMsp << "Energy_keV";
    for (const auto& mat : MATERIAL_ORDER) csvMsp << "," << mat;
    csvMsp << "\n";

    // CSV CSDA Range (um) - LUU Y: day la CSDA range, KHONG PHAI Projected
    // Range nhu PHITS tinh (xem README). Chi dung de THAM KHAO / can tren,
    // khong nen so sanh 1-1 truc tiep voi cot Range cua parse_msp_range.py,
    // dac biet o vung nang luong thap (<100 keV) noi tan xa goc (angular
    // straggling) lam Projected Range < CSDA Range dang ke.
    std::string rangeFile = outFile;
    size_t dotPos = rangeFile.find_last_of('.');
    if (dotPos != std::string::npos) rangeFile = rangeFile.substr(0, dotPos);
    rangeFile += "_csda_range_um.csv";
    std::ofstream csvRange(rangeFile);
    csvRange << "Energy_keV";
    for (const auto& mat : MATERIAL_ORDER) csvRange << "," << mat;
    csvRange << "\n";

    for (double eKeV : energiesKeV) {
        double eMeV = eKeV / 1000.0;
        G4double kinE = eMeV * MeV;

        csvMsp << eKeV;
        csvRange << eKeV;
        for (const auto& mat : MATERIAL_ORDER) {
            G4double dedx = emCal.ComputeTotalDEDX(kinE, "proton", mat);
            G4double dedx_MeV_per_cm = dedx / (MeV/cm);
            G4double density = DENSITY_G_CM3.at(mat);
            G4double msp = dedx_MeV_per_cm / density;   // MeV.cm2/g
            csvMsp << "," << msp;

            G4double csdaRange = emCal.GetCSDARange(kinE, "proton", mat);
            G4double csdaRange_um = csdaRange / um;      // doi ve micromet
            csvRange << "," << csdaRange_um;
        }
        csvMsp << "\n";
        csvRange << "\n";
        std::cout << "Da tinh xong E = " << eKeV << " keV" << std::endl;
    }
    csvMsp.close();
    csvRange.close();
    std::cout << "Da luu CSDA range vao '" << rangeFile << "'" << std::endl;

    std::cout << "\nDa luu ket qua vao '" << outFile << "'" << std::endl;

    delete runManager;
    return 0;
}
