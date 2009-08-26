/*
 * up_finder.cpp
 *
 *  Created on: Aug 26, 2009
 *      Author: magania
 */

#include <getopt.h>

#include <AA.hpp>
#include <Ana.hpp>
#include <AA.hpp>
#include <Trg.hpp>
#include <Muo.hpp>
#include <Ptl.hpp>
#include <DST.hpp>
#include <Field.hpp>
#include <Spot.hpp>
#include <select.hpp>
#include <Monitor.hpp>
#include <V0Finder.hpp>

#include <TTree.h>
#include <TFile.h>

#include "EvtSaver.h"
#include "UpsilonFinder.h"
#include "GammaFinder.h"
#include "XYGammaFinder.h"

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
	AA::findSVFlag = true;
	AA::findV0Flag = false;
	AA::findCascadeFlag = false;
	AA::findJetsFlag = false;
	AA::findJPsiFlag = false;
}

bool nextEvent() {
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();
	return AA::dst.getEvent();
}


void usage(void) {
    printf("\n");
    printf("\033[1mUsage: up_finder [options].\033[m \n");
    printf("options:\n");
    printf("\t \033[1m-i\033[m input_file \t Set the file to process.\n\n");
    printf("\t \033[1m--mc\033[m        \t Process MC events.\n");
    printf("\t \033[1m--elist\033[m     \t Understand input_file as an elist.\n\n");
    printf("Report bugs to <magania@fnal.gov>.\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
    static int mc    = false;
    static int elist = false;

    const char *input_file = 0;

    int c;
    while (1) {
        static struct option long_options[] ={
            /* These options set a flag. */
            {"mc" , no_argument, &mc , 1},
            {"elist" , no_argument, &elist , 1},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "i:",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
            case 'i':
                input_file = optarg;
                break;
            case '?':
                usage();
            default:
                usage();
        }
    }

    if (!input_file){
    	std::cout << argv[0] << " II: Missing input_file " << std::endl;
    	usage();
    }

	AA::newJob();
	AA::newAnaJob();
	AA::newSelectJob();
	params();
	AA::displayParameters();

	bool ok = false;
	if (elist){
		std::cout << argv[0] << " II: Doing event lists " << std::endl;
		ok = AA::dst.eventLst(input_file);
	} else {
		std::cout << argv[0] << " II: Doing sequential mode" << std::endl;
		ok = AA::dst.fileLst(input_file);
	}

	if (!ok) {
		std::cout << argv[0] << " EE: Input stream is not valid or is empty." << std::endl;
		exit(EXIT_FAILURE);
	}

/* ================================        MAIN       ===================================*/
	/* -- We will save all in this tree --*/
	TFile root_file("up.root", "RECREATE");
	TTree tree("tree", "all info.");
	TTree *treeMC=0;

if (mc)
	treeMC = new TTree("treeMC", "all mc info.");
	EvtSaver evt_saver(tree);

	//BhhFinder hh_finder(tree);
	UpsilonFinder upsilon_finder(tree);
	GammaFinder gamma_finder(tree);
	XYGammaFinder xyg_finder(tree, upsilon_finder, gamma_finder);

	/* -- Initilization of geometry, field and beam spot.
	 * Should be done after the fileLst or eventLst initilization -- */
	AA::det.input();
	AA::field.input();
	AA::spot.input();
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();

	while (nextEvent()) {
		AA::setField(); //define beam-spot and field value
		AA::spot.set();
		AA::analyse();
		AA::select(AA::TAG);
        //std::cout << "Run:" << AA::runNumber << " Evt: "  << AA::evtNumber << std::endl;

		if (!upsilon_finder.find())
			continue;
		if (!gamma_finder.find())
			continue;
		if (!xyg_finder.find())
			continue;
		while (xyg_finder.next()){
//			std::cout << xyg_finder.getMass() << ' ' << std::endl;
			xyg_finder.fill();
			evt_saver.fill();
			tree.Fill();
		}
		AA::dst.outEventLst("up_elist");
	}//End while next event.

	tree.Write();
if (mc)
	treeMC->Write();

	root_file.Write();
	root_file.Close();
	std::cout << argv[0] << " II: up_finder ended succesfully." << std::endl;
	std::exit(EXIT_SUCCESS);
}
