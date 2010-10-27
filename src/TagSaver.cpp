/*
 * EvtSaver.cpp
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#include <TagSaver.h>

TagSaver::TagSaver(TTree &tree){
  //Old Tagging
  tree.Branch("dilution"              , &dilution              , "dilution/D"              );//
  tree.Branch("dilution_defined"              , &dilution_defined              , "dilution_defined/I"              );//

 //New Tagging
  tree.Branch("newtag_ost_defined"             , &newtag_ost_defined             , "newtag_ost_defined/I"             );//
  tree.Branch("newtag_muon_defined"            , &newtag_muon_defined            , "newtag_muon_defined/I"            );//
  tree.Branch("newtag_electron_defined"        , &newtag_electron_defined        , "newtag_electron_defined/I"        );//
  tree.Branch("newtag_electronJet_defined"     , &newtag_electronJet_defined     , "newtag_electronJet_defined/I"     );//
  tree.Branch("newtag_electronSVCharge_defined", &newtag_electronSVCharge_defined, "newtag_electronSVCharge_defined/I");//
  tree.Branch("newtag_muonJet_defined"         , &newtag_muonJet_defined         , "newtag_muonJet_defined/I"         );//Muon's jet charge
  tree.Branch("newtag_muonPtRel_defined"       , &newtag_muonPtRel_defined       , "newtag_muonPtRel_defined/I"       );//
  tree.Branch("newtag_muonSVCharge_defined"    , &newtag_muonSVCharge_defined    , "newtag_muonSVCharge_defined/I"    );//
  tree.Branch("newtag_svCharge_defined"        , &newtag_svCharge_defined        , "newtag_svCharge_defined/I"        );//
  tree.Branch("newtag_sv_Pt_defined"           , &newtag_sv_Pt_defined           , "newtag_sv_Pt_defined/I"           );//
  tree.Branch("newtag_Opp_eventCharge_defined" , &newtag_Opp_eventCharge_defined , "newtag_Opp_eventCharge_defined/I" );//
  //tagSame
  tree.Branch("newtag_sst_defined"             , &newtag_sst_defined             , "newtag_sst_defined/I"              );//
  tree.Branch("newtag_minDeltaR_defined"       , &newtag_minDeltaR_defined       , "newtag_minDeltaR_defined/I"       );//
  tree.Branch("newtag_QjetPt_defined"          , &newtag_QjetPt_defined          , "newtag_QjetPt_defined/I"          );//Jet charge weighted with Pt
  //tagEventCharge
  tree.Branch("newtag_SST_eventCharge_defined" , &newtag_SST_eventCharge_defined , "newtag_SST_eventCharge_defined/I" );//

  //New Tagging
  tree.Branch("newtag_ost"             , &newtag_ost             , "newtag_ost/D"              );//
  tree.Branch("newtag_muon"            , &newtag_muon            , "newtag_muon/D"            );//
  tree.Branch("newtag_electron"        , &newtag_electron        , "newtag_electron/D"        );//
  tree.Branch("newtag_electronJet"     , &newtag_electronJet     , "newtag_electronJet/D"     );//
  tree.Branch("newtag_electronSVCharge", &newtag_electronSVCharge, "newtag_electronSVCharge/D");//
  tree.Branch("newtag_muonJet"         , &newtag_muonJet         , "newtag_muonJet/D"         );//Muon's jet charge
  tree.Branch("newtag_muonPtRel"       , &newtag_muonPtRel       , "newtag_muonPtRel/D"       );//
  tree.Branch("newtag_muonSVCharge"    , &newtag_muonSVCharge    , "newtag_muonSVCharge/D"    );//
  tree.Branch("newtag_svCharge"        , &newtag_svCharge        , "newtag_svCharge/D"        );//
  tree.Branch("newtag_sv_Pt"           , &newtag_sv_Pt           , "newtag_sv_Pt/D"           );//
  tree.Branch("newtag_Opp_eventCharge" , &newtag_Opp_eventCharge , "newtag_Opp_eventCharge/D" );//
  //tagSame
  tree.Branch("newtag_sst"             , &newtag_sst             , "newtag_sst/D"            );//
  tree.Branch("newtag_minDeltaR"       , &newtag_minDeltaR       , "newtag_minDeltaR/D"       );//
  tree.Branch("newtag_QjetPt"          , &newtag_QjetPt          , "newtag_QjetPt/D"          );//Jet charge weighted with Pt
  //tagEventCharge
  tree.Branch("newtag_SST_eventCharge" , &newtag_SST_eventCharge , "newtag_SST_eventCharge/D" );//

}

TagSaver::~TagSaver(){
	return;
}

void TagSaver::fill(const HepVector &v, const AA::Vrt *pv, const AA::PtlLst *plst, double ptcut){
  tagComb.fill (v,pv,plst,ptcut);
  //Old Tagging
  dilution = tagComb.dilution();
  dilution_defined = tagComb.defined();
  //New Tagging
  newtag_ost_defined              = tagComb.tagOpposite()->defined();
  newtag_muon_defined             = tagComb.tagOpposite()->muon()->defined();//
  newtag_electron_defined         = tagComb.tagOpposite()->electron()->defined();//
  newtag_electronJet_defined      = tagComb.tagOpposite()->electronJet()->defined();//
  newtag_electronSVCharge_defined = tagComb.tagOpposite()->electronSVCharge()->defined();//
  newtag_muonJet_defined          = tagComb.tagOpposite()->muonJet()->defined();//Muon's jet charge
  newtag_muonPtRel_defined        = tagComb.tagOpposite()->muonPtRel()->defined();//
  newtag_muonSVCharge_defined     = tagComb.tagOpposite()->muonSVCharge()->defined();//
  newtag_svCharge_defined         = tagComb.tagOpposite()->svCharge()->defined();//
  newtag_sv_Pt_defined            = tagComb.tagOpposite()->sv_Pt()->defined();//
  newtag_Opp_eventCharge_defined  = tagComb.tagOpposite()->eventCharge()->defined();//
  //tagSame
  newtag_sst_defined              = tagComb.tagSame()->defined();
  newtag_minDeltaR_defined        = tagComb.tagSame()->minDeltaR()->defined();//
  newtag_QjetPt_defined           = tagComb.tagSame()->QjetPt()->defined();//Jet charge weighted with Pt
  //tagEventCharge
  newtag_SST_eventCharge_defined  = tagComb.tagEventCharge()->defined();//

  newtag_ost              = tagComb.tagOpposite()->dilution();

  newtag_muon             = tagComb.tagOpposite()->muon()->value();//
  newtag_electron         = tagComb.tagOpposite()->electron()->value();//
  newtag_electronJet      = tagComb.tagOpposite()->electronJet()->value();//
  newtag_electronSVCharge = tagComb.tagOpposite()->electronSVCharge()->value();//
  newtag_muonJet          = tagComb.tagOpposite()->muonJet()->value();//Muon's jet charge
  newtag_muonPtRel        = tagComb.tagOpposite()->muonPtRel()->value();//
  newtag_muonSVCharge     = tagComb.tagOpposite()->muonSVCharge()->value();//
  newtag_svCharge         = tagComb.tagOpposite()->svCharge()->value();//
  newtag_sv_Pt            = tagComb.tagOpposite()->sv_Pt()->value();//
  newtag_Opp_eventCharge  = tagComb.tagOpposite()->eventCharge()->value();//
  //tagSame
  newtag_sst               = tagComb.tagSame()->dilution();
  newtag_minDeltaR        = tagComb.tagSame()->minDeltaR()->value();//
  newtag_QjetPt           = tagComb.tagSame()->QjetPt()->value();//Jet charge weighted with Pt
  //tagEventCharge
  newtag_SST_eventCharge  = tagComb.tagEventCharge()->dilution();//

}
