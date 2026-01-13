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
/// \file include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4GenericMessenger.hh"

class G4Box;
class G4Material;
class G4VPhysicalVolume;


class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
    	DetectorConstruction();
    	~DetectorConstruction() override;

    	G4VPhysicalVolume* Construct() override;
		G4bool      fCheckOverlaps; // option to activate checking of volumes overlaps

  	private:
    	static void DefineMaterials();
    	G4VPhysicalVolume* DefineVolumes();
		void BuildConfig_Water(G4LogicalVolume* worldLV, G4double phantomSize);
		void BuildConfig_Slab(G4LogicalVolume* worldLV, G4double phantomSize, G4double waterThickness,
			G4double boneThickness, G4double lungThickness);
		void BuildConfig_Insert(G4LogicalVolume* worldLV, G4double phantomSize, G4ThreeVector bonePos,
			G4double boneThickness);
		void SelectPhantom(const G4String& name);
		G4String fConfig;
		G4GenericMessenger* fMessenger;

};

#endif
