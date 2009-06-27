/*
 * BsJPsiPhiMCFinder.h
 *
 *  Created on: Jun 15, 2009
 *      Author: magania
 */
#include <DecayMC.h>
#include <TTree.h>

#ifndef BSJPSIPHIMCFINDER_H_
#define BSJPSIPHIMCFINDER_H_

class BsJPsiPhiMCFinder {
public:
	BsJPsiPhiMCFinder(TTree &tree, TTree &treeMC);
	bool find();
	void fill(Ptl* muplus, Ptl* muminus, Ptl* kplus, Ptl* kminus);

private:
	Int_t mc_nbs, mc_nbsbar, mc_run, mc_evt, mc_match, mc_match_sw;
	PtlMC *bs_gen, *jpsi_gen, *muplus_gen, *muminus_gen, *phi_gen, *kplus_gen, *kminus_gen;
	Ptl *muplus_reco, *muminus_reco, *kplus_reco, *kminus_reco;

};

#endif /* BSJPSIPHIMCFINDER_H_ */
