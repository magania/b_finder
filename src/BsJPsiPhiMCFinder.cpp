/*
 * BsJPsiPhiMCFinder.cpp
 *
 *  Created on: Jun 15, 2009
 *      Author: magania
 */

#include <BsJPsiPhiMCFinder.h>

BsJPsiPhiMCFinder::BsJPsiPhiMCFinder(TTree &tree, TTree &treeMC) {
	treeMC.Branch("mc_nbs", &mc_nbs, "mc_nbs/I");
	treeMC.Branch("mc_nbsbar", &mc_nbsbar, "mc_nbsbar/I");
	treeMC.Branch("mc_run", &mc_run, "mc_run/I");
	treeMC.Branch("mc_evt", &mc_evt, "mc_evt/I");
	treeMC.Branch("mc_match", &mc_match, "mc_match/I");

        tree.Branch("match_muplus_chi2",  &match_muplus_chi2,  "match_muplus_chi2/D");
        tree.Branch("match_muminus_chi2", &match_muminus_chi2, "match_muminus_chi2/D");
        tree.Branch("match_kplus_chi2",   &match_kplus_chi2,   "match_kplus_chi2/D");
        tree.Branch("match_kminus_chi2",  &match_kminus_chi2,   "match_kminus_chi2/D");

	tree.Branch("mc_match", &mc_match, "mc_match/I");
	tree.Branch("mc_match_sw", &mc_match_sw, "mc_match_sw/I");
}

bool BsJPsiPhiMCFinder::find(){
	mc_nbs = DecayMC::findAndAsignDecay(531, &bs_gen,
			443, &jpsi_gen, 333, &phi_gen,
			-13, &muplus_gen, 13, &muminus_gen,
			321, &kplus_gen, -321,&kminus_gen);
	mc_nbsbar = DecayMC::findAndAsignDecay(-531, &bs_gen,
			443, &jpsi_gen, 333, &phi_gen,
			-13, &muplus_gen, 13, &muminus_gen,
			-321, &kplus_gen,321, &kminus_gen);
	if ( mc_nbs == 0 && mc_nbsbar ==0 ){
	        muplus_gen=0;
	        muminus_gen=0;
	        kplus_gen=0;
	        kminus_gen=0;
	}

	mc_run = runNumber;
	mc_evt = evtNumber;
        mc_match = 0;
        if ( (mc_nbs==1&&mc_nbsbar==0) || (mc_nbs==0&&mc_nbsbar==1) )
	  mc_match = (int) DecayMC::match(muplus_gen, &muplus_reco,
	         		muminus_gen, &muminus_reco,
	         		kplus_gen, &kplus_reco,
		        	kminus_gen, &kminus_reco, 99);

	return mc_match;
}

void BsJPsiPhiMCFinder::fill(Ptl* muplus, Ptl* muminus, Ptl* kplus, Ptl* kminus){
	PtlMC*  muplus_inv, *muminus_inv, *kplus_inv, *kminus_inv;
	muplus_inv = 0; muminus_inv = 0; kplus_inv = 0; kminus_inv = 0;
	muplus_inv  = DecayMC::imatch(muplus, match_muplus_chi2);
	muminus_inv = DecayMC::imatch(muminus, match_muminus_chi2);
        kplus_inv   = DecayMC::imatch(kplus, match_kplus_chi2);
        kminus_inv  = DecayMC::imatch(kminus, match_kminus_chi2);

       if (muplus_inv && muminus_inv && kplus_inv && kminus_inv) {
	    mc_match    =(int)(muplus_inv == muplus_gen && muminus_inv == muminus_gen
    	  	               && kplus_inv == kplus_gen && kminus_inv == kminus_gen);
	    mc_match_sw =(int)(muplus_inv == muplus_gen && muminus_inv == muminus_gen
    	           	       && kplus_inv == kminus_gen && kminus_inv == kplus_gen);
	} else {
	 	mc_match=0;
		mc_match_sw=0;
	}
//    mc_match = (int)(muplus == muplus_reco && muminus == muminus_reco
//                         && kplus == kplus_reco && kminus == kminus_reco);
//    mc_match_sw = (int)(muplus == muplus_reco && muminus == muminus_reco
//                         && kplus == kminus_reco && kminus == kplus_reco);
}
