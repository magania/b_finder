/*
 * BdJPsiKstarMCFinder.cpp
 *
 *  Created on: Jul 17, 2009
 *      Author: magania
 */

#include <BdJPsiKstarMCFinder.h>

BdJPsiKstarMCFinder::BdJPsiKstarMCFinder(TTree &tree, TTree &treeMC) {
	treeMC.Branch("mc_nbd", &mc_nbd, "mc_nbd/I");
	treeMC.Branch("mc_nbdbar", &mc_nbdbar, "mc_nbdbar/I");
	treeMC.Branch("mc_run", &mc_run, "mc_run/I");
	treeMC.Branch("mc_evt", &mc_evt, "mc_evt/I");
	treeMC.Branch("mc_match", &mc_match, "mc_match/I");

//	tree.Branch("mc_inverse", &mc_inverse, "mc_inverse/I");
	tree.Branch("mc_match", &mc_match, "mc_match/I");
	tree.Branch("mc_match_sw", &mc_match_sw, "mc_match_sw/I");
//        tree.Branch("mc_match_fw", &mc_match_fw, "mc_match_fw/I");
}

bool BdJPsiKstarMCFinder::find(){
	mc_nbd = DecayMC::findAndAsignDecay(511, &bd_gen,
			443, &jpsi_gen, 313, &kstar_gen,
			-13, &muplus_gen, 13, &muminus_gen,
			321, &kaon_gen, -211,&pion_gen);
	mc_nbdbar = DecayMC::findAndAsignDecay(-511, &bd_gen,
			443, &jpsi_gen, 313, &kstar_gen,
			-13, &muplus_gen, 13, &muminus_gen,
			-321, &kaon_gen,211, &pion_gen);
	mc_run = runNumber;
	mc_evt = evtNumber;
        mc_match = 0;
        if ( (mc_nbd==1&&mc_nbdbar==0) || (mc_nbd==0&&mc_nbdbar==1) )
	  mc_match = (int) DecayMC::match(muplus_gen, &muplus_reco,
	         		muminus_gen, &muminus_reco,
	         		kaon_gen, &kaon_reco,
		        	pion_gen, &pion_reco, 99);

	return mc_match;
}

void BdJPsiKstarMCFinder::fill(Ptl* muplus, Ptl* muminus, Ptl* kaon, Ptl* pion){
//	PtlMC*  muplus_inv, *muminus_inv, *kplus_inv, *kminus_inv;
//	muplus_inv = 0; muminus_inv = 0; kplus_inv = 0; kminus_inv = 0;
//	muplus_inv  = DecayMC::imatch(muplus);
//	muminus_inv = DecayMC::imatch(muminus);
//    kplus_inv   = DecayMC::imatch(kplus);
//    kminus_inv  = DecayMC::imatch(kminus);

//    mc_inverse = (int)(muplus_inv && muminus_inv && kplus_inv && kminus_inv);
//    mc_match =   (int)(muplus_inv == muplus_gen && muminus_inv == muminus_gen
//    		        && kplus_inv == kplus_gen && kminus_inv == kminus_inv);
//    mc_match_sw = (int)(muplus_inv == muplus_gen && muminus_inv == muminus_gen
//    	           	 && kplus_inv == kminus_gen && kminus_inv == kplus_inv);
    mc_match = (int)(muplus == muplus_reco && muminus == muminus_reco
                         && kaon == kaon_reco && pion == pion_reco);
    mc_match_sw = (int)(muplus == muplus_reco && muminus == muminus_reco
                         && kaon == pion_reco && pion == kaon_reco);
}
