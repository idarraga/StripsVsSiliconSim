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
// $Id: exampleN02.cc,v 1.16 2009-10-30 14:59:59 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02DetectorConstruction.hh"
#include "ExN02PhysicsList.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "ExN02SteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "CLHEP/Random/RanecuEngine.h"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include <string>
#include <vector>

using namespace std;

#include <time.h>
#include <sys/types.h>
#include <unistd.h>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SaveRandomSeed(long int seed, TString name);
bool CheckPars(int argc, char ** argv);

int main(int argc, char ** argv)
{

	// Pars
	if ( ! CheckPars(argc, argv) ) return 1;

	// Seed the random number generator manually
	time_t rawtime;
	time(&rawtime);
	G4long myseed = G4long(rawtime);
	G4cout << "[INFO] The local time : " << myseed << G4endl;
	// Get the process id and use it to offset the time
	G4long pid = (G4long) getpid();
	G4cout << "[INFO] The pid { getpid() } : " << pid << G4endl;
	// Offset the local time
	myseed += pid;
	// Finally this will be the seed
	G4cout << "[INFO] The random seed (local time + pid): " << myseed << G4endl;

	//choose the Random engine
	RanecuEngine * reng = new CLHEP::RanecuEngine();
	reng->setSeed(myseed);
	G4Random::setTheEngine(reng);

	// Check the experiment type entered by the user
	experiment_type expt = (experiment_type) atoi(argv[1]);

	// Output root file and save the random seed
	TString outputfn =  "output_";
	TString output_ph_fn =  "output_";

	if(expt == __SETUP_STRIPS) {
		outputfn += "Strips.root";
		output_ph_fn += "Phantom_Strips.root";
		SaveRandomSeed(myseed, "Strips");
	} else if(expt == __SETUP_SILICON1) {
		outputfn += "Si.root";
		output_ph_fn += "Phantom_Si.root";
		SaveRandomSeed(myseed, "Si");
	} else {
		G4cout << "[ERR ] Setup requested [" << (G4int)expt << "] can't be found.  Giving up." << G4endl;
		return 1;
	}

	// Output
	Outdata * output_data = new Outdata;
	TFile * of = new TFile(outputfn, "RECREATE");
	TTree * oT = new TTree("T","T");

	oT->Branch("totalScatteringAngle", &(output_data->totalScatteringAngle), "totalScatteringAngle/D");
	oT->Branch("scatteredAngle", &(output_data->scatteredAngle) );
	oT->Branch("scatteredPhi", &(output_data->scatteredPhi) );
	oT->Branch("edep", &(output_data->edep) );
	oT->Branch("noniedep", &(output_data->noniedep) );
	oT->Branch("processName", &(output_data->processName) );

	// Output Phantom as SD
	OutPhantomData * output_phantom_data = new OutPhantomData;
	TFile * ophf = new TFile(output_ph_fn, "RECREATE");
	TTree * ophT = new TTree("TPhantom","TPhantom");
	ophT->Branch("edep", &(output_phantom_data->edep) );

	ophT->Branch("x", &(output_phantom_data->x) );
	ophT->Branch("y", &(output_phantom_data->y) );
	ophT->Branch("z", &(output_phantom_data->z) );

	// Run manager
	//
	G4RunManager * runManager = new G4RunManager;
	//G4MTRunManager * runManager = new G4MTRunManager();
	//runManager->SetNumberOfThreads(2);
	// User Initialization classes (mandatory)
	//
	ExN02DetectorConstruction* detector = new ExN02DetectorConstruction(expt, ophT, output_phantom_data);
	runManager->SetUserInitialization(detector);

	// Physics List
	G4PhysListFactory factory;
	G4VModularPhysicsList * phys = 0;
	//G4String physName = "FTFP_BERT";
	G4String physName = "QGSP_BERT_EMY";
	// Reference PhysicsList via its name
	phys = factory.GetReferencePhysList(physName);
	runManager->SetUserInitialization(phys);

	// User Action classes
	//
	G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction(detector);
	runManager->SetUserAction(gen_action);
	//
	G4UserRunAction* run_action = new ExN02RunAction;
	runManager->SetUserAction(run_action);
	//
	G4UserEventAction* event_action = new ExN02EventAction(oT, output_data);
	runManager->SetUserAction(event_action);
	//
	G4UserSteppingAction* stepping_action = new ExN02SteppingAction(output_data);
	runManager->SetUserAction(stepping_action);


	// Initialize G4 kernel
	//
	runManager->Initialize();

#ifdef G4VIS_USE
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif    

	// Get the pointer to the User Interface manager
	//
	G4UImanager * UImanager = G4UImanager::GetUIpointer();

	if (argc >= 3)   // batch mode
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[2];
		UImanager->ApplyCommand(command+fileName);
	}
	else           // interactive mode : define UI session
	{
#ifdef G4UI_USE
		//G4UIExecutive * ui = new G4UIExecutive(argc,argv);
		G4UIterminal * ui = new G4UIterminal();
#ifdef G4VIS_USE
		UImanager->ApplyCommand("/control/execute vis.mac");
#endif
		ui->SessionStart();
		delete ui;
#endif

#ifdef G4VIS_USE
		delete visManager;
#endif     
	}

	// Free the store: user actions, physics_list and detector_description are
	//                 owned and deleted by the run manager, so they should not
	//                 be deleted in the main() program !

	// First file
	of->cd();
	oT->Write();
	of->Close();

	// Second file
	ophf->cd();
	ophT->Write();
	ophf->Close();

	delete output_data;
	delete output_phantom_data;

#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;

	G4cout << "[DONE]" << G4endl;

	return 0;
}

void SaveRandomSeed(long int seed, TString name) {

	TString ofn = "randseed_";
	ofn += name;
	ofn += ".txt";
	std::ofstream ofs (ofn.Data(), std::ofstream::out);
	ofs << seed << endl;
	ofs.close();

}

bool CheckPars(int argc, char ** argv) {

	// Input parameters
	if ( argc < 2 ) {
		G4cout << "use: " << G4endl;
		G4cout << "		" << argv[0] << "experiment(1=strips,2=silicon)  macro" << G4endl;
		return false;
	}

	return true; // OK
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

