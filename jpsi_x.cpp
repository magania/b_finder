/*
 * jpsi_finder.cpp
 *
 *  Created on: Jul 17, 2009
 *      Author: magania
 */

#include <AA.hpp>
#include <Monitor.hpp>
#include <DST.hpp>
#include <Field.hpp>
#include <Spot.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>

#include "JPsiFinder.h"

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
	AA::findV0Flag = false;
	AA::findCascadeFlag = false;
	AA::findJetsFlag = false;
	AA::findJPsiFlag = true;
	AA::findSVFlag = false;
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

        TFile root_file("jpsi_x.root", "RECREATE");
	TTree tree("tree", "all info.");	
	JPsiFinder jpsi_finder(tree);

	/* -- Initilization of geometry, field and beam spot.
	 * Should be done after the fileLst or eventLst initilization -- */
	AA::det.input();
	AA::field.input();
	AA::spot.input();
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();

        unsigned long int evc =0;

        int n_bins = 50;
        double histo_bins[2*n_bins+2];
        double min_pt = 0.5;
        histo_bins[n_bins+1] = min_pt;
        for ( int i=0; i<n_bins; i++)
           histo_bins[n_bins+2+i] = histo_bins[n_bins+1+i] + exp(i/20)*0.1;
        for ( int i=0; i<=n_bins; i++)
           histo_bins[n_bins-i] = -histo_bins[n_bins+1+i];

        for (int i=0; i<(2*n_bins+2); i++) 
         std::cout << histo_bins[i] << " ";
        std::cout << std::endl;

        
        TH2D *jpsi_histo = new TH2D("jpsi_histo", "jpsi_histo", 2*n_bins+1, histo_bins, 100, 0.0, 6.28);
	while (nextEvent()) {
                if (evc++ > 10000) break;
		AA::setField(); //define beam-spot and field value
		AA::spot.set();
		AA::analyse();
		AA::select(AA::TAG);
		if (!jpsi_finder.find())
			continue;
                while (jpsi_finder.next()){
                  jpsi_histo->Fill(jpsi_finder.getMuPlus().q()*jpsi_finder.getMuPlus().pt(), jpsi_finder.getMuPlus().phi());
                  jpsi_histo->Fill(jpsi_finder.getMuMinus().q()*jpsi_finder.getMuMinus().pt(), jpsi_finder.getMuMinus().phi());
		}
/*
		TString name = "histo_";
                name+= evc;   	        
                TH2D *histo = new TH2D(name, name, 2*n_bins+1, histo_bins, 100, 0.0, 6.28);
                const AA::PtlLst ptl_lst = *AA::ptlBox.particles();
                //if (ptl_lst.size() < 250 ) continue;
 	        for (AA::PtlLstCIt p = ptl_lst.begin(); p != ptl_lst.end(); ++p){
                  histo->Fill((*p)->q()*(*p)->pt(), (*p)->phi());
                  if ( abs((*p)->q()) != 1 ) { std::cout << "q: " <<  (*p)->q() << std::endl; exit(1); }
                  //std::cout << "q:" << (*p)->q() << " pt:" << (*p)->pt() << " phi:" << (*p)->phi() << std::endl;
                }
                histo->Print();
*/
	}//End while next event.

        root_file.Write();
        root_file.Close();
	std::cout << argv[0] << " II: jpsi_finder ended succesfully." << std::endl;
	std::exit(EXIT_SUCCESS);
}
