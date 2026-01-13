//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"

#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"
#include "G4PSEnergyDeposit.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),
    fCheckOverlaps(true),
    fConfig("water")
{
    fMessenger = new G4GenericMessenger(this, "/phantom/", "Phantom geometry control");
    fMessenger->DeclareMethod("select", &DetectorConstruction::SelectPhantom)
  .SetCandidates("water slab insert")
  .SetGuidance("Select phantom configuration.");

    DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    //clean out old geometry before updating
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    // Define volumes
    return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SelectPhantom(const G4String& name)
{
    fConfig = name;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4cout << "Phantom selected: " << fConfig << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    // standard material defined using NIST Manager
    G4NistManager *nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_WATER");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_B-100_BONE");

    // build low density lung equivalent
    // Elements
    auto* elH  = nistManager->FindOrBuildElement("H");
    auto* elC  = nistManager->FindOrBuildElement("C");
    auto* elN  = nistManager->FindOrBuildElement("N");
    auto* elO  = nistManager->FindOrBuildElement("O");
    auto* elNa = nistManager->FindOrBuildElement("Na");
    auto* elP  = nistManager->FindOrBuildElement("P");
    auto* elS  = nistManager->FindOrBuildElement("S");
    auto* elCl = nistManager->FindOrBuildElement("Cl");
    auto* elK = nistManager->FindOrBuildElement("K");

    auto compressedLung = new G4Material("compressedLung",0.385 *g/cm3,9);
    compressedLung -> AddElement(elH,0.103);
    compressedLung -> AddElement(elC,0.107);
    compressedLung -> AddElement(elN,0.032);
    compressedLung -> AddElement(elO,0.746);
    compressedLung -> AddElement(elNa,0.002);
    compressedLung -> AddElement(elP,0.002);
    compressedLung -> AddElement(elS,0.003);
    compressedLung -> AddElement(elCl,0.003);
    compressedLung -> AddElement(elK,0.002);
    compressedLung->GetIonisation()->SetMeanExcitationEnergy(75.3 *eV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildConfig_Water(G4LogicalVolume *worldLV, G4double phantomSize) {
    G4Material *waterMaterial = G4Material::GetMaterial("G4_WATER");
    auto* phantSolid   = new G4Box("phantom", phantomSize/2., phantomSize/2., phantomSize/2.);
    auto* phantLogical = new G4LogicalVolume(phantSolid, waterMaterial, "phantom");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), phantLogical, "phantom",
        worldLV, false, 0);

    // Set step limit in phantom
    auto stepLimit = new G4UserLimits(0.5*mm);
    phantLogical->SetUserLimits(stepLimit);

    auto *waterVisAtt = new G4VisAttributes(G4Colour(0, 0.63, 1));
    waterVisAtt->SetVisibility(true);
    waterVisAtt->SetDaughtersInvisible(false);
    phantLogical->SetVisAttributes(waterVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildConfig_Slab(G4LogicalVolume *worldLV, G4double phantomSize, G4double waterThickness,
    G4double boneThickness, G4double lungThickness) {

    // Set step limit in phantom
    auto stepLimit = new G4UserLimits(0.5*mm);

    // Materials in use
    G4Material *waterMaterial = G4Material::GetMaterial("G4_WATER");
    auto* boneMaterial  = G4Material::GetMaterial("G4_B-100_BONE");
    auto lungMaterial = G4Material::GetMaterial("compressedLung");

    // Definitions
    // Water
    auto* phantSolid   = new G4Box("phantom", phantomSize/2., phantomSize/2., phantomSize/2.);
    auto* phantLogical = new G4LogicalVolume(phantSolid, waterMaterial, "phantom");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), phantLogical, "phantom",
        worldLV, false, 0);

    phantLogical->SetUserLimits(stepLimit);

    auto *waterVisAtt = new G4VisAttributes(G4Colour(0, 0.63, 1));
    waterVisAtt->SetVisibility(true);
    waterVisAtt->SetDaughtersInvisible(false);
    phantLogical->SetVisAttributes(waterVisAtt);

    // Bone
    auto* boneSolid   = new G4Box("bone", boneThickness/2., phantomSize/2., phantomSize/2.);
    auto* boneLogical = new G4LogicalVolume(boneSolid, boneMaterial, "bone");
    new G4PVPlacement(nullptr,
      G4ThreeVector(-phantomSize/2 + waterThickness + boneThickness/2, 0, 0),
      boneLogical, "bone", phantLogical, false, 0, fCheckOverlaps);

    boneLogical->SetUserLimits(stepLimit);

    auto* boneVisAtt = new G4VisAttributes(G4Colour(0.94, 0.95, 0.75));
    boneVisAtt->SetVisibility(true);
    boneLogical->SetVisAttributes(boneVisAtt);

    // Lung
    auto* lungSolid   = new G4Box("lung", lungThickness/2., phantomSize/2., phantomSize/2.);
    auto* lungLogical = new G4LogicalVolume(lungSolid, lungMaterial, "lung");
    new G4PVPlacement(nullptr,
      G4ThreeVector(-phantomSize/2 + waterThickness + boneThickness + lungThickness/2, 0, 0),
      lungLogical, "lung", phantLogical, false, 0, fCheckOverlaps);

    lungLogical->SetUserLimits(stepLimit);

    auto* lungVisAtt = new G4VisAttributes(G4Colour(0.94, 0.85, 0.94));
    lungVisAtt->SetVisibility(true);
    lungLogical->SetVisAttributes(lungVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildConfig_Insert(G4LogicalVolume *worldLV, G4double phantomSize,  G4ThreeVector bonePos,
            G4double boneThickness) {
    // Set step limit in phantom
    auto stepLimit = new G4UserLimits(0.5*mm);

    // Materials in use
    G4Material *waterMaterial = G4Material::GetMaterial("G4_WATER");
    auto* boneMaterial  = G4Material::GetMaterial("G4_B-100_BONE");

    // Definitions
    // Water
    auto* phantSolid   = new G4Box("phantom", phantomSize/2., phantomSize/2., phantomSize/2.);
    auto* phantLogical = new G4LogicalVolume(phantSolid, waterMaterial, "phantom");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), phantLogical, "phantom",
        worldLV, false, 0);

    phantLogical->SetUserLimits(stepLimit);

    auto *waterVisAtt = new G4VisAttributes(G4Colour(0, 0.63, 1));
    waterVisAtt->SetVisibility(true);
    waterVisAtt->SetDaughtersInvisible(false);
    phantLogical->SetVisAttributes(waterVisAtt);

    // Bone
    auto* boneSolid   = new G4Box("bone", boneThickness/2., phantomSize/4., phantomSize/2.);
    auto* boneLogical = new G4LogicalVolume(boneSolid, boneMaterial, "bone");
    new G4PVPlacement(nullptr, bonePos,
      boneLogical, "bone", phantLogical, false, 0, fCheckOverlaps);

    boneLogical->SetUserLimits(stepLimit);

    auto* boneVisAtt = new G4VisAttributes(G4Colour(0.94, 0.95, 0.75));
    boneVisAtt->SetVisibility(true);
    boneLogical->SetVisAttributes(boneVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{
    G4cout << "\n### DefineVolumes() called. fConfig = " << fConfig << " ###\n" << G4endl;
    // Geometry parameters
    G4double phantomSize = 20*cm;
    G4double worldSizeXY = phantomSize*5;
    G4double worldSizeZ = phantomSize*5;

    // For slab phantom
    G4double waterSlabThickness = 2*cm;
    G4double boneSlabThickness = 1*cm;
    G4double lungSlabThickness = 2*cm;

    // For insert phantom
    G4double boneInsertThickness = 2*cm;
    G4double depth = 3*cm;
    G4ThreeVector boneInsertPos = G4ThreeVector(-phantomSize/2 + depth + boneInsertThickness/2,
        phantomSize/4., 0.);


    G4Material *defaultMaterial = G4Material::GetMaterial("G4_AIR");

    //
    // WORLD
    //
    G4VSolid *worldSolid
        = new G4Box("World", worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2);

    auto *worldLogical = new G4LogicalVolume(worldSolid, defaultMaterial, "World");

    G4VPhysicalVolume *worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogical,
            "World", nullptr, false, 0, fCheckOverlaps);

    //
    // PHANTOM - Call config builder
    //
    if (fConfig == "water") {
        BuildConfig_Water(worldLogical, phantomSize);
    }
    else if (fConfig == "slab") {
        BuildConfig_Slab(worldLogical, phantomSize, waterSlabThickness, boneSlabThickness, lungSlabThickness);
    }
    else if (fConfig == "insert") {
        BuildConfig_Insert(worldLogical, phantomSize, boneInsertPos, boneInsertThickness);
    }
    else {
        G4Exception("DetectorConstruction::DefineVolumes", "BadConfig", FatalException,
                    ("Unknown /phantom/select option: " + fConfig).c_str());
    }

    //
    // Visualization attributes
    //
    auto *worldVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    worldVisAtt->SetVisibility(false);
    worldLogical->SetVisAttributes(worldVisAtt);


    //
    // Always return the physical World
    //
    return worldPV;
}
