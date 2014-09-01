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
// $Id: ExN02DetectorConstruction.hh,v 1.10 2008-09-22 16:41:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN02DetectorConstruction_h
#define ExN02DetectorConstruction_h 1


#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include <vector>
using namespace std;

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "ExN02MagneticField.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;
class ExN02DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef enum {
	__SETUP_RESERVED = 0,
	__SETUP_STRIPS,
	__SETUP_SILICON1
} experiment_type;

class ExN02DetectorConstruction : public G4VUserDetectorConstruction
{
public:

	ExN02DetectorConstruction();
	ExN02DetectorConstruction(experiment_type);
	~ExN02DetectorConstruction();

public:

	G4VPhysicalVolume * Construct();
	G4VPhysicalVolume * Construct_SetupGas();
	G4VPhysicalVolume * Construct_SetupSilicon();

	void BuildAGasDetector(G4LogicalVolume *, G4String name);
	void BuildASiDetector(G4LogicalVolume *, G4String name);

	G4double GetWorldFullLength()   {return gD->fWorldLength;};

	void setTargetMaterial (G4String);
	void setChamberMaterial(G4String);
	void SetMagField(G4double);
	void SetMaxStep (G4double);

private:

	G4Box*             solidWorld;    //
	G4LogicalVolume*   logicWorld;    //
	G4VPhysicalVolume* physiWorld;    //

	G4Box*             solidDMEContainer;    //
	G4LogicalVolume*   logicDMEContainer;    //
	G4VPhysicalVolume* physiDMEContainer;    //
	G4Box*             solidDMEContainer2;    //
	G4LogicalVolume*   logicDMEContainer2;    //
	G4VPhysicalVolume* physiDMEContainer2;    //

	// For the silicon detectors
	G4Box*             solidSiBox;   //
	G4LogicalVolume*   logicSiBox;   //
	G4VPhysicalVolume* physiSiBox;   //
	G4Box*             solidSiBox2;   //
	G4LogicalVolume*   logicSiBox2;   //
	G4VPhysicalVolume* physiSiBox2;   //

	vector<G4Box*>             solidUpStreamDetector;   //
	vector<G4LogicalVolume*>   logicUpStreamDetector;   //
	vector<G4VPhysicalVolume*> physiUpStreamDetector;   //

	vector<G4Box*>             solidStripsEnv;   //
	vector<G4LogicalVolume*>   logicUSStripsEnv; //
	vector<G4VPhysicalVolume*> physiUSStripsEnv; //
	vector<G4LogicalVolume*>   logicDSStripsEnv; //
	vector<G4VPhysicalVolume*> physiDSStripsEnv; //

	vector<G4Box*>             solidDownStreamDetector;  //
	vector<G4LogicalVolume*>   logicDownStreamDetector;  //
	vector<G4VPhysicalVolume*> physiDownStreamDetector;  //

	vector<G4Box*>             solidKapton;    //
	vector<G4LogicalVolume*>   logicKaptonUS;  //
	vector<G4VPhysicalVolume*> physiKaptonUS;  //
	vector<G4LogicalVolume*>   logicKaptonDS;  //
	vector<G4VPhysicalVolume*> physiKaptonDS;  //

	vector<G4Box*>             solidStrips;  //
	vector<G4LogicalVolume*>   logicStripsUS;  //
	vector<G4VPhysicalVolume*> physiStripsUS;  //
	vector<G4LogicalVolume*>   logicStripsDS;  //
	vector<G4VPhysicalVolume*> physiStripsDS;  //

	vector<G4VPVParameterisation*> stripsParam; // pointer to chamber parameterisation

	G4Tubs*             solidWaterPhantom;   //
	G4LogicalVolume*   logicWaterPhantom;   //
	G4VPhysicalVolume* physiWaterPhantom;   //

	G4UserLimits* stepLimit;             // pointer to user step limits=
	ExN02MagneticField* fpMagField;   // pointer to the magnetic field
	ExN02DetectorMessenger* detectorMessenger;  // pointer to the Messenger

	//G4double fWorldLength;            // Full length of the world volume
	//G4int    NbOfChambers;            // Nb of chambers in the tracker region
	//G4double stripWidth;            // width of the chambers
	//G4double stripsPitch;	       // distance between chambers

	// Decides which experiment to consider
	experiment_type exp_type;

	typedef struct {

		G4double detectorEnvelopeHx;
		G4double detectorEnvelopeHy;
		G4double detectorEnvelopeHz;

		G4double kaptonHx;
		G4double kaptonHy;
		G4double kaptonHz;

		G4double stripsHx;
		G4double stripsHy;
		G4double stripsHz;

		G4double SiBoxHx;
		G4double SiBoxHy;
		G4double SiBoxHz;

		G4double DMEContainerHx;
		G4double DMEContainerHy;
		G4double DMEContainerHz;

		G4double phantomRadius;
		G4double phantomLengthHz;
		G4double posPhantomZ;

		G4double extraAfterBeamCenter;
		G4double HalfWorldLength;

		G4double posUpStreamDetectorZ;
		G4double posDownStreamDetectorZ;
		G4double posDMEContainer1Z;
		G4double posDMEContainer2Z;
		G4double posSiBox1Z;
		G4double posSiBox2Z;

		G4double fWorldLength;
		G4int    NbOfChambers;
		G4double stripWidth;
		G4double stripsPitch;

	} geobits;

	geobits * gD;

	typedef struct {
		G4Material * Air;
		G4Material * Kapton;
		G4Material * Al;
		G4Material * Water;
		G4Material * CO2;
		G4Material * DME;
		G4Material * GasMixture;
		G4Material * Si;
	} materialbits;

	materialbits * mD;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
