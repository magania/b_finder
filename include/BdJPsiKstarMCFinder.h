/*
 * BdJPsiKstarMCFinder.h
 *
 *  Created on: Jul 17, 2009
 *      Author: magania
 */
#include <DecayMC.h>
#include <TTree.h>

#ifndef BDJPSIKSTARMCFINDER_H_
#define BDJPSIKSTARMCFINDER_H_

class BdJPsiKstarMCFinder {
public:
	BdJPsiKstarMCFinder(TTree &tree, TTree &treeMC);
	bool find();
	void fill(Ptl* muplus, Ptl* muminus, Ptl* kaon, Ptl* pion);

private:
	Int_t mc_nbd, mc_nbdbar, mc_run, mc_evt, mc_match, mc_match_sw;
	PtlMC *bd_gen, *jpsi_gen, *muplus_gen, *muminus_gen, *kstar_gen, *kaon_gen, *pion_gen;
	Ptl *muplus_reco, *muminus_reco, *kaon_reco, *pion_reco;

};

#endif /* BDJPSIKSTARMCFINDER_H_ */
