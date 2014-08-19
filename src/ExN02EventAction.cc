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
// $Id: ExN02EventAction.cc,v 1.11 2006-06-29 17:48:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02EventAction::ExN02EventAction() :
m_outputTree(0), m_outputData(0)
{

}

ExN02EventAction::ExN02EventAction(TTree * T, Outdata * od)
{

	m_outputTree = T;
	m_outputData = od;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02EventAction::~ExN02EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02EventAction::EndOfEventAction(const G4Event* evt)
{
	G4int event_id = evt->GetEventID();

	// get number of stored trajectories
	//
	G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
	G4int n_trajectories = 0;
	if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

	//G4cout << "N primary vertex = " << evt->GetNumberOfPrimaryVertex() << G4endl;
	//G4PrimaryVertex * pV = evt->GetPrimaryVertex();
	//G4cout << "next ==> " << pV->GetNext() << G4endl;

	// periodic printing
	//
	if (event_id < 100 || event_id%100 == 0) {
		G4cout << ">>> Event " << evt->GetEventID() << G4endl;
		G4cout << "    " << n_trajectories
				<< " trajectories stored in this event." << G4endl;
	}

	/*
  // Before clearing make a resume here
  // Check the process and extract only the name checking on what repeats
  // map<Process name, number of repetitions>
  map<TString, int> processMap;
  map<TString, int>::iterator it;

  vector<TString>::iterator itr = m_outputData->processName.begin();
  vector<TString>::iterator itrE = m_outputData->processName.end();
  for( ; itr != itrE ; itr++ ) {
	  it = processMap.find(*itr);
	  if ( it == processMap.end() ) { // if not found start it up
		  processMap[*itr] = 1;
	  } else {
		  processMap[*itr]++;
	  }
  }
  // Now shorten the output vector.  No 1 to 1 correspondance anymore !
  // We need to sort by value and not by key.
  // Flip the map first (becomes a multimap !)
  std::multimap<int, TString> flippedProcessMap = flip_map(processMap);
  std::multimap<int, TString>::iterator itrF = flippedProcessMap.begin();
  // Clear info now and refill it in order
  m_outputData->processName.clear();
  // Fill
  for ( ; itrF != flippedProcessMap.end() ; itrF++ ) {
	  m_outputData->processName.push_back( itrF->second );
  }
	 */

	// First store only the resulting angle
	/*
	unsigned int maxAL = m_outputData->scatteredAngle.size();
	double angleSum = 0.;
	for(unsigned int iAL = 0 ; iAL < maxAL ; iAL++ ) {

		if ( m_outputData->scatteredPhi[iAL] >= 0. ) {
			angleSum += m_outputData->scatteredAngle[iAL];
		} else {
			angleSum -= m_outputData->scatteredAngle[iAL];
		}
		//G4cout << m_outputData->scatteredAngle[iAL] << " --> phi : " << m_outputData->scatteredPhi[iAL] << G4endl;
	}
	m_outputData->scatteredAngle.clear();
	*/
	unsigned int nStepsSaved = m_outputData->steps.size();
	G4ThreeVector firstPoint =  m_outputData->steps[ 0 ];
	G4ThreeVector lastPoint  =  m_outputData->steps[ nStepsSaved - 1 ];
	G4ThreeVector diffFirstLast = lastPoint - firstPoint;

	m_outputData->scatteredAngle.clear();
	m_outputData->scatteredPhi.clear();

	m_outputData->scatteredAngle.push_back( CLHEP::pi - (diffFirstLast.theta()/radian) );

	// Now store process names by order in edep.  Map will make the key unique
	map<TString, double> processEnergyMap;
	map<TString, double>::iterator pEM_itr;

	unsigned int iEPMax = (unsigned int) m_outputData->edep.size();
	unsigned int iEP = 0;
	for( ; iEP < iEPMax ; iEP++ ) {
		pEM_itr = processEnergyMap.find( m_outputData->processName[iEP] );
		if ( pEM_itr == processEnergyMap.end() ) {	// If not found start it up
			processEnergyMap[ m_outputData->processName[iEP] ]  = m_outputData->edep[iEP];
		} else { 								// Otherwise add
			processEnergyMap[ m_outputData->processName[iEP] ] += m_outputData->edep[iEP];
		}
	}
	// Now flip it
	multimap<double, TString> processEnergyMap_flipped;
	map<TString, double>::iterator itrEM = processEnergyMap.begin();
	map<TString, double>::iterator itrEM_E = processEnergyMap.end();
	for ( ; itrEM != itrEM_E ; itrEM++ ) {
		processEnergyMap_flipped.insert ( pair<double, TString>(itrEM->second, itrEM->first ) );
	}
	// This map should be sorted by energy now.  In ascendent order !
	// Copy it to the corresponding vectors
	m_outputData->edep.clear();
	m_outputData->processName.clear();

	multimap<double, TString>::iterator i = processEnergyMap_flipped.begin();
	multimap<double, TString>::iterator iE = processEnergyMap_flipped.end();
	for( ; i != iE ; i++) {
		m_outputData->edep.push_back( i->first );
		m_outputData->processName.push_back( i->second );
	}
	// Reverse -> Make it descendent order.
	reverse( m_outputData->edep.begin(), m_outputData->edep.end() );
	reverse( m_outputData->processName.begin(), m_outputData->processName.end() );

	// End of event, fill the tree and clear variables
	m_outputTree->Fill();
	m_outputData->scatteredAngle.clear();
	m_outputData->edep.clear();
	m_outputData->processName.clear();
	// Used but not persistified
	m_outputData->steps.clear();
	// Not used
	m_outputData->totalScatteringAngle = 0.;
	m_outputData->scatteredPhi.clear();
	m_outputData->noniedep.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
