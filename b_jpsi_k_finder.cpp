/*
 * b_jpsi_k_finder.cpp
 *
 * Created on: Mar 16, 2009
 * Author: magania
 */

#include <EvtSaver.h>
#include <PhiFinder.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <PtlFinder.h>
#include <JPsiFinder.h>
#include <V0Finder.h>
#include <TagSaver.h>

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>

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
//#if defined(P17) || defined(P20)
//#include "IPcal.hpp"
//#endif

#include <TTree.h>
#include <TFile.h>

using namespace std;
using namespace AA;

double MASS_B_MIN = 4.9;
double MASS_B_MAX = 5.7;

double PDG_B_MASS = 5.27917;

void params() {
	//Parameters for Real Data and MC
	AA::debugLevel = 0;
	// AA::firstEventIn = 6539424; // Optional first event to be processed
	AA::monitor.set("AA_ElmBox_Maps", 2, 2);
	AA::monitor.set("AA_State5_Chi2", 10, 2);
	// AA::monitor.set("AA_DST_Field",100,1);

	AA::xBeamDefault = -0.030;
	AA::yBeamDefault = +0.030;

	AA::nTrkMinVrt = 3;
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
		std::cout << " " << argv[0] << " e file" << std::endl;
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

	//Initilization of geometry, field and beam spot. Should be done
	//after the fileLst or eventLst initilization
	AA::det.input();
	AA::field.input();
	AA::spot.input();
//#ifdef P17
//	AA:ipcalib.input();
//#endif
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();
}

bool nextEvent() {
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();
	return AA::dst.getEvent();
}

int main(int argc, char** argv) {
	TFile root_file("b.root", "RECREATE");
	TTree tree("tree", "all info.");

	double dlxy, dvlxy;
	double mb, vmb, emb;

	tree.Branch("b_dlxy", &dlxy, "b_dlxy/D");
	tree.Branch("b_dvlxy", &dvlxy, "b_dvlxy/D");
	tree.Branch("b_mass", &mb, "b_mass/D");
	tree.Branch("b_mass_error", &emb, "b_mass_error/D");

	EvtSaver evt_saver(tree);
	JPsiFinder jpsi_finder(tree);
	PtlFinder k_finder("k", tree, true);
	PtlSaver b_saver("b", tree, false, false);
	VrtSaver vb_saver("b_vrt", tree);
	VrtSaver vp_saver("p_vrt", tree);
	TagSaver tag_saver(tree);

	start_ana(argc, argv);

	while (nextEvent()) {
		AA::setField(); //define beam-spot and field value
		AA::spot.set();
		AA::analyse();
		AA::select(AA::TAG);

		if (!jpsi_finder.find())
			continue;
		if (!k_finder.find())
			continue;

		jpsi_finder.begin();
		k_finder.begin();
		while (jpsi_finder.next()){
			Ptl* jpsi = &jpsi_finder.getJPsi();
			Ptl* muplus = &jpsi_finder.getMuPlus();
			Ptl* muminus = &jpsi_finder.getMuMinus();

			if ( muplus->pt() < 1.5 || muminus->pt() < 1.5) //CUT
				continue;

			double jpsi_mass = jpsi_finder.getMass();

			if ( jpsi_mass < 2.9 || jpsi_mass > 3.3 )  //CUT
				continue;

			if ( muplus->nSMT() < 1 || muminus->nSMT() < 1) //CUT
				continue;


			while (k_finder.next()) {
				Ptl* k = &k_finder.getPtl();

				if (muplus == k || muplus == k)
					continue;

				if( k->pt() < 1.5 || k->nSMT() < 1)    //CUT
					continue;

				//Make vertex of these 5 particles
				TrkLst track_list;
				track_list.push_back(muplus);
				track_list.push_back(muminus);
				track_list.push_back(k);

				Vrt b_vrt;
				if (!b_vrt.fill(jpsi->decayVertex()->x(), &track_list))
					continue;
				if (!b_vrt.filter())
					continue;
				if (b_vrt.size() != 3)
					continue;

				// Construct new particle (B) from these 5 tracks
				Ptl b;
				if (!b.combine(&b_vrt, 0))
					continue;
				if (!b.associateToVrt(&AA::vrtPBox))
					continue;

				if (b_vrt.chi2() > 15) //CUT
					continue;

				// Find primary vertex with minimal distance to B vertex
				Vrt* b_pv = 0;
				if (!b_vrt.associateToVrt(&AA::vrtPBox, b_pv))
					continue;

				vector<AA::PType> types(3);
				types[0] = AA::MU_PLUS;
				types[1] = AA::MU_MINUS;
				types[2] = AA::K_PLUS;

				//Origin and decays vertex of B
				const Vrt* pvB = b.decayVertex();
				const Vrt* pvp = b.primaryVertex();

				// Guennadi Vertex Constraint gvc
				// Bool_t gvc = false;
				{
					//Apply vertex constraint
					HepVector v = pvB->x();
					HepSymMatrix vv = pvB->v();

					//Define constraints
					CnsLst clst;
					Constraint cns;

					//Constraint on J/psi mass
					cns.type = AA::CNS_MASS;
					cns.lst.clear();
					cns.lst.push_back(muminus);
					cns.lst.push_back(muplus);
					cns.value1 = AA::mass[AA::JPSI];
					cns.masses.clear();
					cns.masses.push_back(AA::mass[AA::MU_PLUS]);
					cns.masses.push_back(AA::mass[AA::MU_MINUS]);
					clst.push_back(cns);

					//Constraint on direction of B
					cns.type = AA::CNS_ORIGIN;
					cns.lst.clear();
					cns.lst.push_back(muplus);
					cns.lst.push_back(muminus);
					cns.lst.push_back(k);
					cns.pvp = pvp;
					clst.push_back(cns);

					//Apply constraint
					if (!b.applyConstraints(&clst))
						continue;
				}

				//Decay length and uncertainty (exclude decay products from primary vertex)
				PtlLst lstDecay;
				lstDecay.push_back(muplus);
				lstDecay.push_back(muminus);
				lstDecay.push_back(k);
				b.decayLengthProper(PDG_B_MASS, dlxy, dvlxy, &lstDecay);
				//double exy = sqrt(fabs(vlxy));
				//Momentum of decay products and uncertainty
				//(taking into account the applied constraints)
				//The order of particles is:
				//moms[0] - mu+
				//moms[1] - mu-
				//moms[2] - K+
				//moms[3] - K-
				vector<HepVector> moms;
				vector<HepSymMatrix> vmoms;
				if (!b.momChildren(moms, vmoms))
					continue;

//				if (dlxy < 0.01)    //CUT
//					continue;

				//mass of Bs and uncertainty
				b.mass(types, mb, vmb);
				emb = sqrt(fabs(vmb));

				if (mb < MASS_B_MIN || MASS_B_MAX < mb)
					continue;

				/* --------------------- savers -------------------------*/
				b_saver.fill(b);
				vb_saver.fill(*pvB);
				vp_saver.fill(*pvp);
				evt_saver.fill();
				jpsi_finder.fill();
				k_finder.fill();
				tag_saver.fill(b.mom(), pvp, &lstDecay);  // Tagging info.
				/* --------------------- savers ------------------------*/

				tree.Fill();
				AA::dst.outEventLst("b_elist");
			}
		}
	}//End while next event.

	tree.Write();
	root_file.Write();
	root_file.Close();
	std::cout << argv[0] << " II: jpsi_k ended succesfully." << std::endl;
	std::exit(EXIT_SUCCESS);
}

