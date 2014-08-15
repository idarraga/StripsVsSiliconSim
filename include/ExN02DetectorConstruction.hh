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
     G4VPhysicalVolume * Construct_SetupStrips();
     G4VPhysicalVolume * Construct_SetupSilicon1();

     G4double GetWorldFullLength()   {return fWorldLength;}; 
     
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

     G4Box*             solidUpStreamDetector;   //
     G4LogicalVolume*   logicUpStreamDetector;   //
     G4VPhysicalVolume* physiUpStreamDetector;   //

     G4Box*             solidStripsEnv;   //
     G4LogicalVolume*   logicUSStripsEnv; //
     G4VPhysicalVolume* physiUSStripsEnv; //
     G4LogicalVolume*   logicDSStripsEnv; //
     G4VPhysicalVolume* physiDSStripsEnv; //

     G4Box*             solidDownStreamDetector;  //
     G4LogicalVolume*   logicDownStreamDetector;  //
     G4VPhysicalVolume* physiDownStreamDetector;  //
     
     G4Box*             solidKapton;    //
     G4LogicalVolume*   logicKaptonUS;  //
     G4VPhysicalVolume* physiKaptonUS;  //
     G4LogicalVolume*   logicKaptonDS;  //
     G4VPhysicalVolume* physiKaptonDS;  //
     
     G4Box*             solidStrips;  //
     G4LogicalVolume*   logicStripsUS;  //
     G4VPhysicalVolume* physiStripsUS;  //
     G4LogicalVolume*   logicStripsDS;  //
     G4VPhysicalVolume* physiStripsDS;  //
     

     G4VPVParameterisation* stripsParam; // pointer to chamber parameterisation
     G4UserLimits* stepLimit;             // pointer to user step limits

     ExN02MagneticField* fpMagField;   // pointer to the magnetic field 
     
     ExN02DetectorMessenger* detectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume
     G4int    NbOfChambers;            // Nb of chambers in the tracker region
     G4double stripWidth;            // width of the chambers
     G4double stripsPitch;	       // distance between chambers

     // Decides which experiment to consider
     experiment_type exp_type;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
