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
// $Id: ExN02ChamberParameterisation.cc,v 1.9 2006-06-29 17:47:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02ChamberParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02ChamberParameterisation::ExN02ChamberParameterisation(  
		G4int    NoChambers,
		G4double startY,          //  Z of center of first
		G4double spacingY,        //  Z spacing of centers
		G4double widthChamber,
		G4double hx,
		G4double hz )
{
	fNoChambers =  NoChambers;
	fStartY     =  startY;
	fSpacing    =  spacingY;
	fHalfWidth  =  widthChamber*0.5;
	fhx = hx;
	fhz = hz;

	if( NoChambers > 0 ){

		if (spacingY < widthChamber) {
			G4Exception("ExN02ChamberParameterisation::ExN02ChamberParameterisation()",
					"InvalidSetup", FatalException,
					"Width>Spacing");
		}
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02ChamberParameterisation::~ExN02ChamberParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02ChamberParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
	//G4cout << "Compute transformation on copy No: " << copyNo << G4endl;
	G4double      Yposition = fStartY + (copyNo) * fSpacing;
	G4ThreeVector origin(0, Yposition, 0);
	physVol->SetTranslation(origin);
	//physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02ChamberParameterisation::ComputeDimensions
(G4Box& trackerChamber, const G4int /*copyNo*/, const G4VPhysicalVolume*) const
{
	//G4cout << "Compute dimensions on copy No: " << copyNo << G4endl;
	trackerChamber.SetXHalfLength(fhx);
	trackerChamber.SetYHalfLength(fHalfWidth);
	trackerChamber.SetZHalfLength(fhz);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
