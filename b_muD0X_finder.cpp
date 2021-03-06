/*
 * b_muD0X_finder.cpp
 *
 *  Created on: Oct 22, 2010
 *      Author: magania
 */

#include <AA.hpp>

#include <TTree.h>
#include <TFile.h>

#include "EvtSaver.h"
#include "D0Finder.h"
#include "PtlFinder.h"
#include "BMuD0XFinder.h"
#include "TagSaver.h"

void params() {
	//Parameters for Real Data and MC
	AA::debugLevel = 0;
	//  AA::firstEventIn = 6539424; // Optional first event to be processed
	AA::monitor.set("AA_ElmBox_Maps", 2, 2);
	AA::monitor.set("AA_State5_Chi2", 10, 2);
	//  AA::monitor.set("AA_DST_Field",100,1);

	AA::xBeamDefault = -0.030;
	AA::yBeamDefault = +0.030;

	AA::nTrkMinVrt = 5; //It defines the minimal number of tracks in the primary vertex.
	AA::correctMomentum = false;

	//parameters for jet finding (see PYTHIA description)
	AA::mstu[46] = 5;
	AA::paru[44] = 15.;

	//definition of input blocks
	AA::blocksIn.clear();
	//definition of output blocks
	AA::blocksOut.clear();

	AA::tuneImpact = false;

	// :)
	AA::findPVFlag = true;
	AA::findV0Flag = true;
	AA::findCascadeFlag = false;
	AA::findJetsFlag = true;
	AA::findJPsiFlag = true;
	AA::findSVFlag = true;
}

void start_ana(int argc, char** argv) {
	if (argc != 3 || !(argv[1][0] == 'f' || argv[1][0] == 'e')) {
		std::cout << "USAGE: " << argv[0] << " f file" << std::endl;
		std::cout << "       " << argv[0] << " e elist" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string file = argv[2];

	AA::newJob();
	AA::newAnaJob();
	AA::newSelectJob();
	params();
	AA::displayParameters();

	bool ok = false;
	switch (argv[1][0]) {
	case 'e':
		std::cout << argv[0] << " II: Doing event lists " << std::endl;
		ok = AA::dst.eventLst(file);
		break;
	case 'f':
		std::cout << argv[0] << " II: Doing sequential mode" << std::endl;
		ok = AA::dst.fileLst(file);
		break;
	default:
		std::cout << argv[0] << " EE: Unknown switch " << argv[1][0]
				<< std::endl;
		exit(EXIT_FAILURE);
	}
	if (!ok) {
		std::cout << argv[0] << " EE: Input stream is not valid or is empty."
				<< std::endl;
		exit(EXIT_FAILURE);
	}
}

bool nextEvent() {
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();
	return AA::dst.getEvent();
}

int main(int argc, char** argv) {
	start_ana(argc, argv);
	/* -- We will save all in this tree --*/
	TFile root_file("bMuD0X.root", "RECREATE");
	TTree tree("tree", "all info.");

	EvtSaver evt_saver(tree);
        
	D0Finder d0_finder(tree);
	PtlFinder pi_finder("dpion", tree, false);
	PtlFinder mu_finder("muon", tree, true);
        BMuD0XFinder b_finder(tree, d0_finder, pi_finder, mu_finder);

	/* -- Initilization of geometry, field and beam spot.
	 * Should be done after the fileLst or eventLst initilization -- */
	AA::det.input();
	AA::field.input();
	AA::spot.input();
#if defined(P17) || defined(MC)
//	AA:ipcalib.input();
#endif
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();

	while (nextEvent()) {
		AA::setField(); //define beam-spot and field value
		AA::spot.set();
		AA::analyse();
		AA::select(AA::TAG);
//                std::cout << "Run:" << AA::runNumber << " Evt: "  << AA::evtNumber << std::endl;

		if (!d0_finder.find())
			continue;
		//std::cout << "JPsi" << std::endl;
		if (!pi_finder.find())
			continue;
                //std::cout << "Phi" << std::endl;
		if (!mu_finder.find())
			continue;
		if (!b_finder.find())
			continue;

		b_finder.begin();
		while (b_finder.next()){
			b_finder.fill();
			evt_saver.fill();
			tree.Fill();
		}
		AA::dst.outEventLst("bMuD0X_elist");
	}//End while next event.

	tree.Write();

	root_file.Write();
	root_file.Close();
	std::cout << argv[0] << " II: b_MuD0X_finder ended succesfully." << std::endl;
	std::exit(EXIT_SUCCESS);
}
