/*
 * TagSaver.h
 *
 *  Created on: Mar 19, 2009
 *      Author: magania & yorch
 */

#ifndef TAGSAVER_H_
#define TAGSAVER_H_

#include <vector>
#include <string>
#include <list>

#include <AA.hpp>
#include <Trg.hpp>
#include <FlavorTag.hpp>
#include <FlavorSSTag.hpp>

#include <TTree.h>

class TagSaver {
public:
	TagSaver(TTree &tree);
	~TagSaver();

	void fill(const HepVector &v, const AA::Vrt *pv, const AA::PtlLst *plst, double ptcut = 0.5);

private:

	//New Tagging
	AA::TagCombinedAll tagComb;
	
	//Old Tagging
	Double_t dilution;
	Int_t dilution_defined;
	
	//tagOpposite
	Double_t newtag_ost;
	Int_t newtag_ost_defined;
	Double_t newtag_muon;//
 	Int_t    newtag_muon_defined;
	Double_t newtag_electron;//
	Int_t    newtag_electron_defined;
	Double_t newtag_electronJet;//
	Int_t    newtag_electronJet_defined;
	Double_t newtag_electronSVCharge;//
        Int_t    newtag_electronSVCharge_defined;
	Double_t newtag_muonJet;//Muon's jet charge
	Int_t    newtag_muonJet_defined;
	Double_t newtag_muonPtRel;//
	Int_t    newtag_muonPtRel_defined;
	Double_t newtag_muonSVCharge;//
	Int_t    newtag_muonSVCharge_defined;
	Double_t newtag_svCharge;//
	Int_t    newtag_svCharge_defined;
	Double_t newtag_sv_Pt;//
	Int_t    newtag_sv_Pt_defined;
	Double_t newtag_Opp_eventCharge;//
	Int_t    newtag_Opp_eventCharge_defined;
	//tagSame
	Double_t newtag_sst;
	Int_t newtag_sst_defined;
	Double_t newtag_minDeltaR;//
	Int_t    newtag_minDeltaR_defined;
	Double_t newtag_QjetPt;//Jet charge weighted with Pt
	Int_t    newtag_QjetPt_defined;
	//tagEventCharge
	Double_t newtag_SST_eventCharge;//
	Int_t    newtag_SST_eventCharge_defined;

};

#endif /* TAGSAVER_H_ */
