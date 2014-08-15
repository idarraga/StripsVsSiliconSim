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
// $Id: ExN02DetectorConstruction.cc,v 1.22 2010-01-22 11:57:03 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02DetectorConstruction.hh"
#include "ExN02DetectorMessenger.hh"
#include "ExN02ChamberParameterisation.hh"
#include "ExN02MagneticField.hh"
#include "ExN02TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02DetectorConstruction::ExN02DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 stepLimit(0), fpMagField(0),
 fWorldLength(0.),
 NbOfChambers(0)
{
	fpMagField = new ExN02MagneticField();
	detectorMessenger = new ExN02DetectorMessenger(this);
}

ExN02DetectorConstruction::ExN02DetectorConstruction(experiment_type et)
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 stepLimit(0), fpMagField(0),
 fWorldLength(0.),
 NbOfChambers(0)
{
	fpMagField = new ExN02MagneticField();
	detectorMessenger = new ExN02DetectorMessenger(this);
	exp_type = et;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02DetectorConstruction::~ExN02DetectorConstruction()
{
	delete fpMagField;
	delete stepLimit;
	delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN02DetectorConstruction::Construct(){

	if(exp_type == __SETUP_STRIPS) {

		physiWorld = Construct_SetupStrips();

	} else if(exp_type == __SETUP_SILICON1) {

		physiWorld =  Construct_SetupSilicon1();

	} else {
		G4cout << "[ERR ] Setup requested [" << (G4int)exp_type << "] can't be found.  Giving up." << G4endl;
	}

	return physiWorld;
}

G4VPhysicalVolume* ExN02DetectorConstruction::Construct_SetupSilicon1() {

	//--------- Material definition ---------
	G4NistManager * nistman = G4NistManager::Instance();
	//nistman->ListMaterials("all");
	// Air
	G4Material * Air = nistman->FindOrBuildMaterial("G4_AIR");
	// Si
	G4Material * Si = nistman->FindOrBuildMaterial("G4_Si");

	// Print all the materials defined.
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	//--------- Sizes of the principal geometrical components (solids)  ---------

	// WARNING ! Kapton and Strip have to fit in the envelope !
	G4double detectorEnvelopeHx = 30*mm;
	G4double detectorEnvelopeHy = 30*mm;
	G4double detectorEnvelopeHz = 100*um;//

	// 10cm of air in the middle, the two Silicon planes + 5cm or air on each side
	fWorldLength = 10*cm + 4*detectorEnvelopeHz + 2*5*cm;

	G4double posUpStreamDetectorZ = 50*mm + detectorEnvelopeHz;
	G4double posDownStreamDetectorZ = -posUpStreamDetectorZ;


	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

	//------------------------------
	// World
	//------------------------------

	G4double HalfWorldLength = 0.5*fWorldLength;

	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
	G4cout << "Computed tolerance = "
			<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
			<< " mm" << G4endl;

	solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	physiWorld = new G4PVPlacement(0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			logicWorld,      // its logical volume
			"World",         // its name
			0,               // its mother  volume
			false,           // no boolean operations
			0);              // copy number

	//------------------------------
	// Upstream detector
	//------------------------------
	G4ThreeVector positionUpStreamDetector = G4ThreeVector(0, 0, posUpStreamDetectorZ);

	solidUpStreamDetector = new G4Box("solidUpStreamDetector", detectorEnvelopeHx, detectorEnvelopeHy, detectorEnvelopeHz);
	logicUpStreamDetector = new G4LogicalVolume(solidUpStreamDetector, Si, "logicUpStreamDetector",0 ,0 ,0);
	physiUpStreamDetector = new G4PVPlacement(0,          // no rotation
			positionUpStreamDetector,  // at (x,y,z)
			logicUpStreamDetector,     // its logical volume
			"UpStreamDetector",        // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	//------------------------------
	// Downstream detector
	//------------------------------
	G4ThreeVector positionDownStreamDetector = G4ThreeVector(0, 0, posDownStreamDetectorZ);
	G4RotationMatrix * pRot = new G4RotationMatrix();
	pRot->rotateY( pi );
	solidDownStreamDetector = new G4Box("solidDownStreamDetector", detectorEnvelopeHx, detectorEnvelopeHy, detectorEnvelopeHz);
	logicDownStreamDetector = new G4LogicalVolume(solidDownStreamDetector, Si, "logicDownStreamDetector",0 ,0 ,0);
	physiDownStreamDetector = new G4PVPlacement(pRot,          // no rotation
			positionDownStreamDetector,  // at (x,y,z)
			logicDownStreamDetector,     // its logical volume
			"DownStreamDetector",        // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number


	G4VisAttributes* detVisAtt = new G4VisAttributes( G4Colour::Blue() );
	detVisAtt->SetForceSolid(true);
	logicUpStreamDetector->SetVisAttributes(detVisAtt);
	logicDownStreamDetector->SetVisAttributes(detVisAtt);


	return physiWorld;
}

G4VPhysicalVolume* ExN02DetectorConstruction::Construct_SetupStrips()
{

	//--------- Material definition ---------
	G4NistManager * nistman = G4NistManager::Instance();
	//nistman->ListMaterials("all");
	// Air
	G4Material * Air = nistman->FindOrBuildMaterial("G4_AIR");
	// KAPTON
	G4Material * Kapton = nistman->FindOrBuildMaterial("G4_KAPTON");
	// Al
	G4Material * Al = nistman->FindOrBuildMaterial("G4_Al");

	// DME, gas --> C2H6O
	G4double a, z;            // z=mean number of protons;
	a = 1.01*g/mole;
	G4Element * elH  = new G4Element("Hydrogen", "H" , z= 1., a);
	a = 12.01*g/mole;
	G4Element * elC  = new G4Element("Carbon"  , "C" , z= 6., a);
	a = 16.00*g/mole;
	G4Element * elO  = new G4Element("Oxygen"  , "O" , z= 8., a);

	G4double DME_density     = 0.0021146*g/cm3;
	G4int DME_ncomponents = 3;
	// temperature and pressure by default are STP
	G4Material * DME = new G4Material("DME", DME_density, DME_ncomponents, kStateGas);
	DME->AddElement(elC, 2);
	DME->AddElement(elH, 6);
	DME->AddElement(elO, 1);

	// Print all the materials defined.
	//
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	//--------- Sizes of the principal geometrical components (solids)  ---------


	// WARNING ! Kapton and Strip have to fit in the envelope !
	G4double detectorEnvelopeHx = 30*mm;
	G4double detectorEnvelopeHy = 30*mm;
	G4double detectorEnvelopeHz = 40*um;//*100;

	G4double kaptonHx = 30*mm;
	G4double kaptonHy = 30*mm;
	G4double kaptonHz = 25*um;//*100;

	G4double stripsHx = 30*mm;
	G4double stripsHy = 30*mm;
	G4double stripsHz = 15*um;//*100;

	// The two detectors fit in the DME container leaving an
	//  effective distance of 4cm in between them.
	G4double DMEContainerHx = 30*mm;
	G4double DMEContainerHy = 30*mm;
	G4double DMEContainerHz = 20*mm + kaptonHz*2 + stripsHz*2;

	// The DME + 5cm on each side
	fWorldLength = 2*DMEContainerHz + 2*50*mm;

	G4double posUpStreamDetectorZ = 20*mm + kaptonHz+stripsHz;
	G4double posDownStreamDetectorZ = -posUpStreamDetectorZ;

	NbOfChambers = 15;
	stripWidth = 2*mm;
	stripsPitch = 2*stripWidth;
	// Make enough chambers to cover 50% of the area
	NbOfChambers = (G4int)(detectorEnvelopeHy*2 / stripsPitch);

	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

	//------------------------------
	// World
	//------------------------------

	G4double HalfWorldLength = 0.5*fWorldLength;

	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
	G4cout << "Computed tolerance = "
			<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
			<< " mm" << G4endl;

	solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	physiWorld = new G4PVPlacement(0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			logicWorld,      // its logical volume
			"World",         // its name
			0,               // its mother  volume
			false,           // no boolean operations
			0);              // copy number

	//------------------------------
	// DME container
	//------------------------------
	G4ThreeVector positionDMEContainer = G4ThreeVector(0, 0, 0);

	solidDMEContainer = new G4Box("solidDMEContainer", DMEContainerHx, DMEContainerHy, DMEContainerHz);
	logicDMEContainer = new G4LogicalVolume(solidDMEContainer, DME, "logicDMEContainer",0 ,0 ,0);
	physiDMEContainer = new G4PVPlacement(0,          // no rotation
			positionDMEContainer,  // at (x,y,z)
			logicDMEContainer,     // its logical volume
			"DMEContainer",        // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	//------------------------------
	// Upstream detector
	//------------------------------
	G4ThreeVector positionUpStreamDetector = G4ThreeVector(0, 0, posUpStreamDetectorZ);

	solidUpStreamDetector = new G4Box("solidUpStreamDetector", detectorEnvelopeHx, detectorEnvelopeHy, detectorEnvelopeHz);
	logicUpStreamDetector = new G4LogicalVolume(solidUpStreamDetector, Air, "logicUpStreamDetector",0 ,0 ,0);
	physiUpStreamDetector = new G4PVPlacement(0,          // no rotation
			positionUpStreamDetector,  // at (x,y,z)
			logicUpStreamDetector,     // its logical volume
			"UpStreamDetector",        // its name
			logicDMEContainer,      // its mother  volume
			false,           // no boolean operations
			0);              // copy number

	//------------------------------
	// Downstream detector
	//------------------------------

	G4ThreeVector positionDownStreamDetector = G4ThreeVector(0, 0, posDownStreamDetectorZ);
	G4RotationMatrix * pRot = new G4RotationMatrix();
	pRot->rotateY( pi );
	solidDownStreamDetector = new G4Box("solidDownStreamDetector", detectorEnvelopeHx, detectorEnvelopeHy, detectorEnvelopeHz);
	logicDownStreamDetector = new G4LogicalVolume(solidDownStreamDetector, Air, "logicDownStreamDetector",0 ,0 ,0);
	physiDownStreamDetector = new G4PVPlacement(pRot,          // no rotation
			positionDownStreamDetector,  // at (x,y,z)
			logicDownStreamDetector,     // its logical volume
			"DownStreamDetector",        // its name
			logicDMEContainer,      // its mother  volume
			false,           // no boolean operations
			0);              // copy number

	//------------------------------
	// Kapton Upstream
	//------------------------------
	G4ThreeVector positionKaptonUS = G4ThreeVector( 0, 0, stripsHz );

	solidKapton = new G4Box("solidKapton", kaptonHx, kaptonHy, kaptonHz);
	logicKaptonUS = new G4LogicalVolume(solidKapton, Kapton, "logicKaptonUS", 0, 0, 0);
	physiKaptonUS = new G4PVPlacement(0,              // no rotation
			positionKaptonUS, // at (x,y,z)
			logicKaptonUS,    // its logical volume
			"KaptonUS",       // its name
			logicUpStreamDetector,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number
	//------------------------------
	// Strips envelope Upstream
	//------------------------------
	G4ThreeVector positionStripsEnvelopeUS = G4ThreeVector( 0, 0, -kaptonHz );
	solidStripsEnv = new G4Box("stripsEnv", stripsHx, stripsHy, stripsHz);
	logicUSStripsEnv = new G4LogicalVolume(solidStripsEnv, Air, "logicUSStripsEnv", 0, 0, 0);
	physiUSStripsEnv = new G4PVPlacement(0,              // no rotation
			positionStripsEnvelopeUS, // at (x,y,z)
			logicUSStripsEnv,    // its logical volume
			"StripsUSEnv",       // its name
			logicUpStreamDetector,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	//------------------------------
	// Kapton Downstream
	//------------------------------
	G4ThreeVector positionKaptonDS = G4ThreeVector( 0, 0, stripsHz );

	logicKaptonDS = new G4LogicalVolume(solidKapton, Kapton, "logicKaptonDS", 0, 0, 0);
	physiKaptonDS = new G4PVPlacement(0,              // no rotation
			positionKaptonDS, // at (x,y,z)
			logicKaptonDS,    // its logical volume
			"KaptonDS",       // its name
			logicDownStreamDetector,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number
	//------------------------------
	// Strips envelope Downstream
	//------------------------------
	G4ThreeVector positionStripsEnvelopeDS = G4ThreeVector( 0, 0, -kaptonHz);
	logicDSStripsEnv = new G4LogicalVolume(solidStripsEnv, Air, "logicDSStripsEnv", 0, 0, 0);
	physiDSStripsEnv = new G4PVPlacement(0,              // no rotation
			positionStripsEnvelopeDS, // at (x,y,z)
			logicDSStripsEnv,    // its logical volume
			"StripsDSEnv",       // its name
			logicDownStreamDetector,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	//------------------------------
	// Strip segments US
	//------------------------------
	solidStrips = new G4Box("solidStrips", stripsHx, stripsHy, stripsHz);
	logicStripsUS = new G4LogicalVolume(solidStrips, Al, "logicStrips", 0, 0, 0);

	G4double firstPosition = -detectorEnvelopeHy + 0.5*stripWidth;
	stripsParam = new ExN02ChamberParameterisation(
			NbOfChambers,          // NoChambers
			firstPosition,         // Y of center of first
			stripsPitch,           // Y spacing of centers
			stripWidth,            // Width Chamber
			stripsHx,
			stripsHz);

	physiStripsUS = new G4PVParameterised(
			"StripsUS",       // their name
			logicStripsUS,      // their logical volume
			logicUSStripsEnv,    // Mother logical volume
			kYAxis,           // Are placed along this axis
			NbOfChambers,     // Number of chambers
			stripsParam, true);     // The parametrisation

	G4cout << "There are " << NbOfChambers << " chambers in the Strips US region" << G4endl;

	//------------------------------
	// Strip segments DS
	//------------------------------
	logicStripsDS = new G4LogicalVolume(solidStrips, Al, "logicStrips", 0, 0, 0);
	physiStripsDS = new G4PVParameterised(
			"StripsDS",       // their name
			logicStripsDS,      // their logical volume
			logicDSStripsEnv,    // Mother logical volume
			kYAxis,           // Are placed along this axis
			NbOfChambers,     // Number of chambers
			stripsParam, true);     // The parametrisation

	G4cout << "There are " << NbOfChambers << " chambers in the Strips DS region" << G4endl;


	//------------------------------------------------
	// Sensitive detectors
	//------------------------------------------------
	/*
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4String trackerChamberSDname = "ExN02/TrackerChamberSD";
	ExN02TrackerSD* aTrackerSD = new ExN02TrackerSD( trackerChamberSDname );
	SDman->AddNewDetector( aTrackerSD );
	logicChamber->SetSensitiveDetector( aTrackerSD );

	//--------- Visualization attributes -------------------------------

	G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	logicWorld  ->SetVisAttributes(BoxVisAtt);
	logicTarget ->SetVisAttributes(BoxVisAtt);
	logicTracker->SetVisAttributes(BoxVisAtt);

	G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	logicChamber->SetVisAttributes(ChamberVisAtt);

	//--------- example of User Limits -------------------------------

	// below is an example of how to set tracking constraints in a given
	// logical volume(see also in N02PhysicsList how to setup the processes
	// G4StepLimiter or G4UserSpecialCuts).

	// Sets a max Step length in the tracker region, with G4StepLimiter
	//
	G4double maxStep = 0.5*ChamberWidth;
	stepLimit = new G4UserLimits(maxStep);
	logicTracker->SetUserLimits(stepLimit);

	// Set additional contraints on the track, with G4UserSpecialCuts
	//
	// G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
	// logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
	//                                               minEkin));
	 */

	G4VisAttributes* kaptonVisAtt = new G4VisAttributes( G4Colour::Yellow() );
	kaptonVisAtt->SetForceSolid(true);
	logicKaptonUS->SetVisAttributes(kaptonVisAtt);
	logicKaptonDS->SetVisAttributes(kaptonVisAtt);

	G4VisAttributes* stripsVisAtt = new G4VisAttributes( G4Colour::Red() );
	stripsVisAtt->SetForceSolid(true);
	logicStripsUS->SetVisAttributes(stripsVisAtt);
	logicStripsDS->SetVisAttributes(stripsVisAtt);

	//logicUSStripsEnv->Set

	//G4VisAttributes* UpStreamDetectorVisAtt = new G4VisAttributes( G4Colour::Red() );
	//UpStreamDetectorVisAtt->SetForceSolid(true);
	//logicUpStreamDetector->SetVisAttributes(UpStreamDetectorVisAtt);

	return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::setTargetMaterial(G4String)
{

	/*
	// search the material by its name
	G4Material* pttoMaterial = G4Material::GetMaterial(materialName);
	if (pttoMaterial)
	{TargetMater = pttoMaterial;
	logicTarget->SetMaterial(pttoMaterial);
	G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
			<< materialName << G4endl;
	}
	 */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::setChamberMaterial(G4String)
{
	/*
	// search the material by its name
	G4Material* pttoMaterial = G4Material::GetMaterial(materialName);
	if (pttoMaterial)
	{ChamberMater = pttoMaterial;
	logicChamber->SetMaterial(pttoMaterial);
	G4cout << "\n----> The chambers are " << ChamberWidth/cm << " cm of "
			<< materialName << G4endl;
	}
	 */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::SetMagField(G4double fieldValue)
{
	fpMagField->SetMagFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::SetMaxStep(G4double maxStep)
{
	if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
