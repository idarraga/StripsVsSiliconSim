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
// $Id: ExN02SteppingAction.cc,v 1.9 2006-06-29 17:48:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02SteppingAction.hh"
#include "G4SteppingManager.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02SteppingAction::ExN02SteppingAction()
{

}

ExN02SteppingAction::ExN02SteppingAction(Outdata * od)
{
	m_outputData = od;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02SteppingAction::UserSteppingAction(const G4Step * step)
{ 

	G4Track * track = step->GetTrack();

	// Consider only the primary here
	if(track->GetTrackID() == 1) {

		G4ThreeVector deltaPos = step->GetDeltaPosition();
		// Theta
		G4double theta = deltaPos.angle();
		G4double phi = deltaPos.phi();
		// This is the interesting angle
		m_outputData->scatteredAngle.push_back( CLHEP::pi - theta ); // comes in radians
		// This direction is not reall important as the system is symetric in this direction
		//m_outputData->scatteredPhi.push_back( phi );

		//
		m_outputData->edep.push_back( step->GetTotalEnergyDeposit()/keV );
		//
		//m_outputData->noniedep.push_back( step->GetNonIonizingEnergyDeposit()/keV );

		//track->Get !!!
		G4StepPoint * postStep = step->GetPostStepPoint();
		G4String procName = " NoProc";
		if( postStep ) {
			const G4VProcess * process = step->GetPostStepPoint()->GetProcessDefinedStep();
			if ( process ) {
				procName = " UserLimit";
				if (process) procName = process->GetProcessName();
				if ( fpSteppingManager->GetfStepStatus() == fWorldBoundary ) procName = "OutOfWorld";
			}
		}
		m_outputData->processName.push_back( TString( procName.c_str() ) );

	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

