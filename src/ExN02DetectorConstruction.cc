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
#include "G4Tubs.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02DetectorConstruction::ExN02DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 stepLimit(0), fpMagField(0)
{
	fpMagField = new ExN02MagneticField();
	detectorMessenger = new ExN02DetectorMessenger(this);
}

ExN02DetectorConstruction::ExN02DetectorConstruction(experiment_type et,  TTree * T, OutPhantomData * od)
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 stepLimit(0), fpMagField(0)
{

	fpMagField = new ExN02MagneticField();
	detectorMessenger = new ExN02DetectorMessenger(this);
	exp_type = et;

	// Output for the SD
	_T = T;
	_od = od;

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

		physiWorld = Construct_SetupGas();

	} else if(exp_type == __SETUP_SILICON1) {

		physiWorld =  Construct_SetupSilicon();

	} else {
		G4cout << "[ERR ] Setup requested [" << (G4int)exp_type << "] can't be found.  Giving up." << G4endl;
	}

	return physiWorld;
}

void ExN02DetectorConstruction::BuildASiDetector(G4LogicalVolume * log_dme, G4String name){

	// Determine the correct id for this detector set
	// I'll rely on one of the vectors
	G4int id = (G4int) solidUpStreamDetector.size();

	//------------------------------
	// Upstream detector
	//------------------------------
	G4ThreeVector positionUpStreamDetector = G4ThreeVector(0, 0, gD->posUpStreamDetectorZ);

	solidUpStreamDetector.push_back( new G4Box("solidUpStreamDetector"+name, gD->detectorEnvelopeHx, gD->detectorEnvelopeHy, gD->detectorEnvelopeHz) );
	logicUpStreamDetector.push_back( new G4LogicalVolume(solidUpStreamDetector[id], mD->Si, "logicUpStreamDetector"+name,0 ,0 ,0) );
	physiUpStreamDetector.push_back( new G4PVPlacement(0,          // no rotation
			positionUpStreamDetector,  // at (x,y,z)
			logicUpStreamDetector[id],     // its logical volume
			"UpStreamDetector"+name,        // its name
			log_dme,      // its mother  volume
			false,           // no boolean operations
			0, true) );              // copy number

	//------------------------------
	// Downstream detector
	//------------------------------
	G4ThreeVector positionDownStreamDetector = G4ThreeVector(0, 0, gD->posDownStreamDetectorZ);
	G4RotationMatrix * pRot = new G4RotationMatrix();
	pRot->rotateY( pi );
	solidDownStreamDetector.push_back( new G4Box("solidDownStreamDetector"+name, gD->detectorEnvelopeHx, gD->detectorEnvelopeHy, gD->detectorEnvelopeHz) );
	logicDownStreamDetector.push_back( new G4LogicalVolume(solidDownStreamDetector[id], mD->Si, "logicDownStreamDetector"+name,0 ,0 ,0) );
	physiDownStreamDetector.push_back( new G4PVPlacement(pRot,          // no rotation
			positionDownStreamDetector,  // at (x,y,z)
			logicDownStreamDetector[id],     // its logical volume
			"DownStreamDetector"+name,        // its name
			log_dme,      // its mother  volume
			false,           // no boolean operations
			0, true) );              // copy number

}

G4VPhysicalVolume* ExN02DetectorConstruction::Construct_SetupSilicon() {

	mD = new materialbits;

	//--------- Material definition ---------
	G4NistManager * nistman = G4NistManager::Instance();
	//nistman->ListMaterials("all");
	// Air
	mD->Air = nistman->FindOrBuildMaterial("G4_AIR");
	// Si
	mD->Si = nistman->FindOrBuildMaterial("G4_Si");
	// Water
	mD->Water = nistman->FindOrBuildMaterial("G4_WATER");

	// Print all the materials defined.
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	//--------- Sizes of the principal geometrical components (solids)  ---------
	gD = new geobits;

	//
	gD->detectorEnvelopeHx = 30*mm;
	gD->detectorEnvelopeHy = 30*mm;
	gD->detectorEnvelopeHz = 100*um;//

	gD->SiBoxHx = 30*mm;
	gD->SiBoxHy = 30*mm;
	gD->SiBoxHz = 5*cm + 200*um;

	//
	gD->extraAfterBeamCenter = gD->SiBoxHz + 5*cm + 20*cm + 5*cm + 2*gD->SiBoxHz + 5*cm;
	gD->fWorldLength = 300*cm + gD->extraAfterBeamCenter;
	gD->HalfWorldLength = 0.5*gD->fWorldLength;

	gD->posSiBox1Z = -gD->HalfWorldLength + 5*cm + 2*gD->SiBoxHz + 5*cm + 20*cm + 5*cm + gD->SiBoxHz;
	gD->posSiBox2Z = -gD->HalfWorldLength + 5*cm + gD->SiBoxHz;

	gD->posUpStreamDetectorZ = 50*mm + gD->detectorEnvelopeHz;
	gD->posDownStreamDetectorZ = -gD->posUpStreamDetectorZ;

	// Phantom
	gD->posPhantomZ = -gD->HalfWorldLength + 5*cm + 2*gD->SiBoxHz + 5*cm + 10*cm;
	gD->phantomRadius = 10*cm;
	gD->phantomLengthHz = 20*cm;

	G4cout << "Half World Length = " << gD->HalfWorldLength/mm << " mm" << G4endl;

	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

	//------------------------------
	// World
	//------------------------------


	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(gD->fWorldLength);
	G4cout << "Computed tolerance = "
			<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
			<< " mm" << G4endl;

	solidWorld= new G4Box("world",gD->HalfWorldLength,gD->HalfWorldLength,gD->HalfWorldLength);
	logicWorld= new G4LogicalVolume( solidWorld, mD->Air, "World", 0, 0, 0);

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
	// Si detector 1
	//------------------------------
	G4ThreeVector positionSiBox1 = G4ThreeVector(0, 0, gD->posSiBox1Z);
	solidSiBox = new G4Box("solidSiBox_Si1", gD->SiBoxHx, gD->SiBoxHy, gD->SiBoxHz);
	logicSiBox = new G4LogicalVolume(solidSiBox, mD->Air, "logicSiBox_Si1", 0,0,0);
	physiSiBox = new G4PVPlacement(0,          // no rotation
			positionSiBox1,  // at (x,y,z)
			logicSiBox,     // its logical volume
			"physiSiBox_Si1",        // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number


	// Build a Si detector
	BuildASiDetector(logicSiBox, "_Si1");

	//------------------------------
	// Si detector 2
	//------------------------------
	G4ThreeVector positionSiBox2 = G4ThreeVector(0, 0, gD->posSiBox2Z);
	solidSiBox2 = new G4Box("solidSiBox_Si2", gD->SiBoxHx, gD->SiBoxHy, gD->SiBoxHz);
	logicSiBox2 = new G4LogicalVolume(solidSiBox2, mD->Air, "logicSiBox_Si2", 0,0,0);
	physiSiBox2 = new G4PVPlacement(0,          // no rotation
			positionSiBox2,  // at (x,y,z)
			logicSiBox2,     // its logical volume
			"physiSiBox_Si2",        // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	// Build a Gas detector
	BuildASiDetector(logicSiBox2, "_Si2");

	//------------------------------
	// Water Phantom
	//------------------------------
	G4ThreeVector positionWaterPhantom = G4ThreeVector( 0, 0, gD->posPhantomZ);
	G4RotationMatrix * pRot = new G4RotationMatrix();
	pRot->rotateY( pi/2. );
	solidWaterPhantom = new G4Tubs("solidWaterPhantom", 0., gD->phantomRadius, gD->phantomLengthHz, 0., 2*pi);
	logicWaterPhantom = new G4LogicalVolume(solidWaterPhantom, mD->Water, "logicWaterPhantom", 0, 0, 0);
	physiWaterPhantom = new G4PVPlacement(pRot,
			positionWaterPhantom, // at (x,y,z)
			logicWaterPhantom,    // its logical volume
			"WaterPhantom",       // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number



	G4VisAttributes* detVisAtt = new G4VisAttributes( G4Colour::Blue() );
	detVisAtt->SetForceSolid(true);
	G4int nVols = (G4int) solidUpStreamDetector.size();
	for(G4int i = 0 ; i < nVols ; i++) {
		logicUpStreamDetector[i]->SetVisAttributes(detVisAtt);
		logicDownStreamDetector[i]->SetVisAttributes(detVisAtt);
	}

	G4VisAttributes* phanVisAtt = new G4VisAttributes( G4Colour::Blue() );
	phanVisAtt->SetForceSolid(true);
	phanVisAtt->SetForceWireframe(true);
	phanVisAtt->SetForceAuxEdgeVisible(true);
	logicWaterPhantom->SetVisAttributes(phanVisAtt);

	return physiWorld;
}

G4VPhysicalVolume* ExN02DetectorConstruction::Construct_SetupGas()
{

	mD = new materialbits;

	//--------- Material definition ---------
	G4NistManager * nistman = G4NistManager::Instance();
	//nistman->ListMaterials("all");
	// Air
	mD->Air = nistman->FindOrBuildMaterial("G4_AIR");
	// KAPTON
	mD->Kapton = nistman->FindOrBuildMaterial("G4_KAPTON");
	// Al
	mD->Al = nistman->FindOrBuildMaterial("G4_Al");
	// Water
	mD->Water = nistman->FindOrBuildMaterial("G4_WATER");
	// C02
	mD->CO2 = nistman->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

	//G4_CARBON_DIOXIDE

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
	mD->DME = new G4Material("DME", DME_density, DME_ncomponents, kStateGas);
	mD->DME->AddElement(elC, 2);
	mD->DME->AddElement(elH, 6);
	mD->DME->AddElement(elO, 1);

	// Mixture for chamber.  Considering an ideal gas.
	G4double Mixture_density = 0.5* ( DME_density + mD->CO2->GetDensity() );
	mD->GasMixture = new G4Material("GasMixture", Mixture_density, 2);
	mD->GasMixture->AddMaterial(mD->DME, 50*perCent);
	mD->GasMixture->AddMaterial(mD->CO2, 50*perCent);

	// Print all the materials defined.
	//
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	//--------- Sizes of the principal geometrical components (solids)  ---------

	gD = new geobits;

	// WARNING ! Kapton and Strip have to fit in the envelope !
	gD->detectorEnvelopeHx = 30*mm;
	gD->detectorEnvelopeHy = 30*mm;
	gD->detectorEnvelopeHz = 40*um;//*100;

	gD->kaptonHx = 30*mm;
	gD->kaptonHy = 30*mm;
	gD->kaptonHz = 25*um;//*100;

	gD->stripsHx = 30*mm;
	gD->stripsHy = 30*mm;
	gD->stripsHz = 15*um;//*100;

	// The two detectors fit in the DME container leaving an
	//  effective distance of 4cm in between them.
	gD->DMEContainerHx = 30*mm;
	gD->DMEContainerHy = 30*mm;
	gD->DMEContainerHz = 20*mm + gD->kaptonHz*2 + gD->stripsHz*2;

	// Phantom
	gD->phantomRadius = 10*cm;
	gD->phantomLengthHz = 20*cm;

	// The DME + 5cm on each side
	//fWorldLength = 2*DMEContainerHz + 2*50*mm;
	// World Length = 3m down to the center of the DME box + half of the DME box +
	// 				  5cm up to the water phantom + 20cm of water phantom +
	//				  5cm after phantom + DME box + 5cm downstream
	gD->extraAfterBeamCenter = gD->DMEContainerHz + 5*cm + 20*cm + 5*cm + 2*gD->DMEContainerHz + 5*cm;
	gD->fWorldLength = 300*cm + gD->extraAfterBeamCenter;
	gD->HalfWorldLength = 0.5*gD->fWorldLength;
	gD->posPhantomZ = -gD->HalfWorldLength + 5*cm + 2*gD->DMEContainerHz + 5*cm + 10*cm;

	G4cout << "Half World Length = " << gD->HalfWorldLength/mm << " mm" << G4endl;

	gD->posUpStreamDetectorZ = 20*mm + gD->kaptonHz+gD->stripsHz;
	gD->posDownStreamDetectorZ = -gD->posUpStreamDetectorZ;
	gD->posDMEContainer1Z = -gD->HalfWorldLength + gD->extraAfterBeamCenter;
	gD->posDMEContainer2Z = -gD->HalfWorldLength + 5*cm + gD->DMEContainerHz;

	gD->stripWidth = 2*mm;
	gD->stripsPitch = 2*gD->stripWidth;
	// Make enough chambers to cover 50% of the area
	gD->NbOfChambers = (G4int)(gD->detectorEnvelopeHy*2 / gD->stripsPitch);

	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

	//------------------------------
	// World
	//------------------------------
	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(gD->fWorldLength);
	G4cout << "Computed tolerance = "
			<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
			<< " mm" << G4endl;

	solidWorld= new G4Box("world",gD->HalfWorldLength,gD->HalfWorldLength,gD->HalfWorldLength);
	logicWorld= new G4LogicalVolume( solidWorld, mD->Air, "World", 0, 0, 0);

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
	// DME container 1
	//------------------------------
	G4ThreeVector positionDMEContainer1 = G4ThreeVector(0, 0, gD->posDMEContainer1Z);

	solidDMEContainer = new G4Box("solidDMEContainer_GAS1", gD->DMEContainerHx, gD->DMEContainerHy, gD->DMEContainerHz);
	logicDMEContainer = new G4LogicalVolume(solidDMEContainer, mD->GasMixture, "logicDMEContainer_GAS1",0 ,0 ,0);
	physiDMEContainer = new G4PVPlacement(0,          // no rotation
			positionDMEContainer1,  // at (x,y,z)
			logicDMEContainer,     // its logical volume
			"DMEContainer_GAS1",        // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	// Build a Gas detector
	BuildAGasDetector(logicDMEContainer, "_GAS1");

	//------------------------------
	// DME container 2
	//------------------------------
	G4ThreeVector positionDMEContainer2 = G4ThreeVector(0, 0, gD->posDMEContainer2Z);
	solidDMEContainer2 = new G4Box("solidDMEContainer_GAS2", gD->DMEContainerHx, gD->DMEContainerHy, gD->DMEContainerHz);
	logicDMEContainer2 = new G4LogicalVolume(solidDMEContainer2, mD->GasMixture, "logicDMEContainer_GAS2", 0 ,0 ,0);
	physiDMEContainer2 = new G4PVPlacement(0,          // no rotation
			positionDMEContainer2,  // at (x,y,z)
			logicDMEContainer2,     // its logical volume
			"DMEContainer_GAS2",        // its name
			logicWorld,      	// its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	// Build a Gas detector
	BuildAGasDetector(logicDMEContainer2, "_GAS2");


	//------------------------------
	// Water Phantom
	//------------------------------
	G4ThreeVector positionWaterPhantom = G4ThreeVector( 0, 0, gD->posPhantomZ);
	G4RotationMatrix * pRot = new G4RotationMatrix();
	pRot->rotateY( pi/2. );
	solidWaterPhantom = new G4Tubs("solidWaterPhantom", 0., gD->phantomRadius, gD->phantomLengthHz, 0., 2*pi);
	logicWaterPhantom = new G4LogicalVolume(solidWaterPhantom, mD->Water, "logicWaterPhantom", 0, 0, 0);
	physiWaterPhantom = new G4PVPlacement(pRot,
			positionWaterPhantom, // at (x,y,z)
			logicWaterPhantom,    // its logical volume
			"WaterPhantom",       // its name
			logicWorld,      // its mother  volume
			false,           // no boolean operations
			0, true);              // copy number

	//------------------------------------------------
	// Sensitive detectors
	//------------------------------------------------
	/*
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4String trackerChamberSDname = "ExN02/TrackerChamberSD";
	ExN02TrackerSD* aTrackerSD = new ExN02TrackerSD( trackerChamberSDname );
	SDman->AddNewDetector( aTrackerSD );
	logicChamber->SetSensitiveDetector( aTrackerSD );


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

	G4int nVols = (G4int) solidUpStreamDetector.size();
	G4VisAttributes* kaptonVisAtt = new G4VisAttributes( G4Colour::Yellow() );
	kaptonVisAtt->SetForceSolid(true);
	G4VisAttributes* stripsVisAtt = new G4VisAttributes( G4Colour::Red() );
	stripsVisAtt->SetForceSolid(true);

	for(G4int i = 0 ; i < nVols ; i++) {
		logicKaptonUS[i]->SetVisAttributes(kaptonVisAtt);
		logicKaptonDS[i]->SetVisAttributes(kaptonVisAtt);
		logicStripsUS[i]->SetVisAttributes(stripsVisAtt);
		logicStripsDS[i]->SetVisAttributes(stripsVisAtt);
	}
	G4VisAttributes* phanVisAtt = new G4VisAttributes( G4Colour::Blue() );
	phanVisAtt->SetForceSolid(true);
	phanVisAtt->SetForceWireframe(true);
	phanVisAtt->SetForceAuxEdgeVisible(true);
	logicWaterPhantom->SetVisAttributes(phanVisAtt);

	return physiWorld;
}

void ExN02DetectorConstruction::BuildAGasDetector(G4LogicalVolume * log_dme, G4String name){

	// Determine the correct id for this detector set
	// I'll rely on one of the vectors
	G4int id = (G4int) solidUpStreamDetector.size();


	//------------------------------
	// Upstream detector
	//------------------------------
	G4ThreeVector positionUpStreamDetector = G4ThreeVector(0, 0, gD->posUpStreamDetectorZ);

	solidUpStreamDetector.push_back( new G4Box("solidUpStreamDetector"+name, gD->detectorEnvelopeHx, gD->detectorEnvelopeHy, gD->detectorEnvelopeHz) );
	logicUpStreamDetector.push_back( new G4LogicalVolume(solidUpStreamDetector[id], mD->Air, "logicUpStreamDetector"+name,0 ,0 ,0) );
	physiUpStreamDetector.push_back( new G4PVPlacement(0,          // no rotation
			positionUpStreamDetector,  // at (x,y,z)
			logicUpStreamDetector[id],     // its logical volume
			"UpStreamDetector"+name,        // its name
			log_dme,      // its mother  volume
			false,           // no boolean operations
			0) );              // copy number

	//------------------------------
	// Downstream detector
	//------------------------------

	G4ThreeVector positionDownStreamDetector = G4ThreeVector(0, 0, gD->posDownStreamDetectorZ);
	G4RotationMatrix * pRot = new G4RotationMatrix();
	pRot->rotateY( pi );
	solidDownStreamDetector.push_back( new G4Box("solidDownStreamDetector"+name, gD->detectorEnvelopeHx, gD->detectorEnvelopeHy, gD->detectorEnvelopeHz) );
	logicDownStreamDetector.push_back( new G4LogicalVolume(solidDownStreamDetector[id], mD->Air, "logicDownStreamDetector"+name,0 ,0 ,0) );
	physiDownStreamDetector.push_back( new G4PVPlacement(pRot,          // no rotation
			positionDownStreamDetector,  // at (x,y,z)
			logicDownStreamDetector[id],     // its logical volume
			"DownStreamDetector"+name,        // its name
			log_dme,      // its mother  volume
			false,           // no boolean operations
			0) );              // copy number

	//------------------------------
	// Kapton Upstream
	//------------------------------
	G4ThreeVector positionKaptonUS = G4ThreeVector( 0, 0, gD->stripsHz );

	solidKapton.push_back( new G4Box("solidKapton"+name, gD->kaptonHx, gD->kaptonHy, gD->kaptonHz) );
	logicKaptonUS.push_back( new G4LogicalVolume(solidKapton[id], mD->Kapton, "logicKaptonUS"+name, 0, 0, 0) );
	physiKaptonUS.push_back( new G4PVPlacement(0,              // no rotation
			positionKaptonUS, // at (x,y,z)
			logicKaptonUS[id],    // its logical volume
			"KaptonUS"+name,       // its name
			logicUpStreamDetector[id],      // its mother  volume
			false,           // no boolean operations
			0, true) );              // copy number
	//------------------------------
	// Strips envelope Upstream
	//------------------------------
	G4ThreeVector positionStripsEnvelopeUS = G4ThreeVector( 0, 0, -gD->kaptonHz );
	solidStripsEnv.push_back( new G4Box("stripsEnv"+name, gD->stripsHx, gD->stripsHy, gD->stripsHz) );
	logicUSStripsEnv.push_back( new G4LogicalVolume(solidStripsEnv[id], mD->Air, "logicUSStripsEnv"+name, 0, 0, 0) );
	physiUSStripsEnv.push_back( new G4PVPlacement(0,              // no rotation
			positionStripsEnvelopeUS, // at (x,y,z)
			logicUSStripsEnv[id],    // its logical volume
			"StripsUSEnv"+name,       // its name
			logicUpStreamDetector[id],      // its mother  volume
			false,           // no boolean operations
			0, true) );              // copy number

	//------------------------------
	// Kapton Downstream
	//------------------------------
	G4ThreeVector positionKaptonDS = G4ThreeVector( 0, 0, gD->stripsHz );

	logicKaptonDS.push_back( new G4LogicalVolume(solidKapton[id], mD->Kapton, "logicKaptonDS"+name, 0, 0, 0) );
	physiKaptonDS.push_back( new G4PVPlacement(0,              // no rotation
			positionKaptonDS, // at (x,y,z)
			logicKaptonDS[id],    // its logical volume
			"KaptonDS"+name,       // its name
			logicDownStreamDetector[id],      // its mother  volume
			false,           // no boolean operations
			0, true) );              // copy number

	//------------------------------
	// Strips envelope Downstream
	//------------------------------
	G4ThreeVector positionStripsEnvelopeDS = G4ThreeVector( 0, 0, -gD->kaptonHz);
	logicDSStripsEnv.push_back( new G4LogicalVolume(solidStripsEnv[id], mD->Air, "logicDSStripsEnv"+name, 0, 0, 0) );
	physiDSStripsEnv.push_back( new G4PVPlacement(0,              // no rotation
			positionStripsEnvelopeDS, // at (x,y,z)
			logicDSStripsEnv[id],    // its logical volume
			"StripsDSEnv"+name,       // its name
			logicDownStreamDetector[id],      // its mother  volume
			false,           // no boolean operations
			0, true) );              // copy number

	//------------------------------
	// Strip segments US
	//------------------------------
	solidStrips.push_back( new G4Box("solidStrips"+name, gD->stripsHx, gD->stripsHy, gD->stripsHz) );
	logicStripsUS.push_back( new G4LogicalVolume(solidStrips[id], mD->Al, "logicStrips"+name, 0, 0, 0) );

	G4double firstPosition = -gD->detectorEnvelopeHy + 0.5*gD->stripWidth;
	stripsParam.push_back( new ExN02ChamberParameterisation(
			gD->NbOfChambers,          // NoChambers
			firstPosition,         // Y of center of first
			gD->stripsPitch,           // Y spacing of centers
			gD->stripWidth,            // Width Chamber
			gD->stripsHx,
			gD->stripsHz) );

	physiStripsUS.push_back( new G4PVParameterised(
			"StripsUS"+name,       // their name
			logicStripsUS[id],      // their logical volume
			logicUSStripsEnv[id],    // Mother logical volume
			kYAxis,           // Are placed along this axis
			gD->NbOfChambers,     // Number of chambers
			stripsParam[id], true) );     // The parametrisation

	G4cout << "There are " << gD->NbOfChambers << " chambers in the Strips US region" << G4endl;

	//------------------------------
	// Strip segments DS
	//------------------------------
	logicStripsDS.push_back( new G4LogicalVolume(solidStrips[id], mD->Al, "logicStrips"+name, 0, 0, 0) );
	physiStripsDS.push_back( new G4PVParameterised(
			"StripsDS"+name,       // their name
			logicStripsDS[id],      // their logical volume
			logicDSStripsEnv[id],    // Mother logical volume
			kYAxis,           // Are placed along this axis
			gD->NbOfChambers,     // Number of chambers
			stripsParam[id], true) );     // The parametrisation

	G4cout << "There are " << gD->NbOfChambers << " chambers in the Strips DS region" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::ConstructSDandField()
{

	// Sensitive NaI
	// _NaILogic
	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	ExN02TrackerSD * Phantom_SD = new ExN02TrackerSD("Phantom_SD");
	Phantom_SD->SetOutput(_T, _od);
	SDman->AddNewDetector( Phantom_SD );
	logicWaterPhantom->SetSensitiveDetector( Phantom_SD );

}
