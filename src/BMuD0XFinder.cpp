/*
 * BMuD0XFinder.cpp
 *
 *  Created on: Oct 22, 2010
 *      Author: magania
 */

#include "BMuD0XFinder.h"

BMuD0XFinder::BMuD0XFinder(TTree &tree, D0Finder &D0, PtlFinder &pion, PtlFinder &muon) :
        tag_saver(tree),
	b_saver("b",tree, false, false),
	vb_saver("b_vrt", tree),
	vp_saver("b_pv_vrt", tree)
{
	index = -1;

	_d0_finder = &D0;
	_pi_finder = &pion;
	_mu_finder = &muon;

	tree.Branch("dstar_mass", &dstar_mass, "dstar_mass/D");
	tree.Branch("b_dl",  &b_dl,  "b_dl/D");
	tree.Branch("b_vdl",  &b_vdl,  "b_vdl/D");
	tree.Branch("b_pdl",  &b_pdl,  "b_pdl/D");
	tree.Branch("b_epdl",  &b_epdl,  "b_epdl/D");
	tree.Branch("d0_dl", &d0_dl, "d0_dl/D");
	tree.Branch("d0_vdl", &d0_vdl, "d0_vdl/D");
}

BMuD0XFinder::~BMuD0XFinder() {
}

void BMuD0XFinder::clean(){
	_boxp.clear();
	_boxv.clear();

        v_d0_index.clear();
	v_mu_index.clear();
        v_b.clear();
	v_mdstar.clear();
}

int BMuD0XFinder::find(){
    //std::cout << "Finding B .. " << std::endl;
	clean();
	_d0_finder->begin();
	_mu_finder->begin();
	while (_mu_finder->next()){
		AA::Ptl* muon = &_mu_finder->getPtl();
                AA::MuoGlb *glbMuon = muon->muon();
                if ( !glbMuon) continue;

                if( muon->ptot() < PTOT_MU_MIN ) continue;                       
                if( muon->pt() < PT_MU_MIN ) continue;                       
                if( muon->nSMT() < MUON_NSMT_MIN || muon->nCFT()<MUON_NCFT_MIN ) continue;                       
                if( !muon->primaryVertex() ) continue;
                if (glbMuon->nseg() < MUON_NSEG_MIN) continue;
                if (glbMuon->chisqloc() < MUON_CHISQLOC_MIN ) continue;

//                std::cout << "Muon" << std::endl;

		while (_d0_finder->next()) {
			AA::Ptl* d0 = &_d0_finder->getD0();
			AA::Ptl* kaon = &_d0_finder->getKaon();
			AA::Ptl* pion = &_d0_finder->getPion();
		
			/* ---------------- Selection ------------- */
                        if ( kaon->chi2Vrt()<PI_K_CHI2VRT_MIN || pion->chi2Vrt()<PI_K_CHI2VRT_MIN ) continue;

                        if (muon->q()*kaon->q() < 0 )
                                continue;

			if (muon == kaon || muon == pion)
				continue;

			if (muon->primaryVertex() != kaon->primaryVertex() || muon->primaryVertex() != pion->primaryVertex())
				continue;

                        if (!muon->jet() || muon->jet() != kaon->jet() || muon->jet() != pion->jet()) continue;

   		        if (!d0->associateToVrt(&AA::vrtPBox) || d0->primaryVertex() != muon->primaryVertex())
  		            continue;

                        if (d0->cosXY() < COSXY_D0_MIN)
                            continue;

//                        std::cout << "First Cuts" << std::endl;

                        double d0lxy,d0vlxy;
                        d0->decayVertex()->distanceXY(muon->primaryVertex(),d0lxy,d0vlxy);
                        if ( d0lxy*d0lxy < SIGMA_B_LXY*d0vlxy ) continue;

                        /* -- B vertex -- */
                        AA::TrkLst track_list;
                        track_list.push_back(d0);
                        track_list.push_back(muon);
                        AA::Vrt *b_vrt = _boxv.newVrt();
                        if (!b_vrt->fill(muon->primaryVertex()->x(), &track_list))
                           continue;
                        if (!b_vrt->filter())
                           continue;
                        if (b_vrt->size() != 2)
                           continue;

                       if (b_vrt->chi2() > CHI2_BVRT_MIN ) continue;
//                      std::cout << "First Cuts 1" << std::endl;

                        double blxy,bvlxy;
                        b_vrt->distanceXY(muon->primaryVertex(),blxy,bvlxy);
                        if (blxy*blxy > d0lxy*d0lxy) continue;
                        
                        double bdlxy,bdvlxy;
                        d0->decayVertex()->distanceXY(b_vrt,bdlxy,bdvlxy);
                        if ( bdlxy*bdlxy < SIGMA_D0_LXY*bdvlxy ) continue;

                      /* -- Construct new particle B from these 2 tracks -- */
                        AA::Ptl* b = _boxp.newPtl();
                        if (!b->combine(b_vrt, muon->q()))
                           continue;

//                      std::cout << "First Cuts 2" << std::endl;
                      /* Compute mass of l + D0 */
                      std::list<AA::Chain*> chb;
           	      b_vrt->chains(&chb);
                      std::list<AA::Chain*>::iterator pc = chb.begin();
                      AA::Chain* pcD = *pc;
 	              AA::Chain* pcl = *(++pc);

                      double xmb = AA::xMass(b_vrt,pcl,pcD,AA::mass[AA::MU_PLUS],_d0_finder->getMass());
//                      std::cout << "First Cuts 3" << std::endl;

                      /* -- cut on B mass -- */
                      if ( xmb < MASS_B_MIN || MASS_B_MAX < xmb)
                          continue;

   		    if (!b->associateToVrt(&AA::vrtPBox) || muon->primaryVertex() != b->primaryVertex())
  		          continue;

                    if (b->cosXY() < COSXY_B_MIN)
                         continue;

                    bool dstarFound = false;
                    double dstarMass = -1.0;
		    _pi_finder->begin();
                    while(_pi_finder->next()){
			AA::Ptl* dpion = &_pi_finder->getPtl();
		
           	        if (dpion == muon || dpion == kaon || dpion == pion)
                           continue;
			
			if (dpion->q() == 0)
				continue;

			if (dpion->primaryVertex() != muon->primaryVertex())
				continue;

			//Find the B0 vertex
      			AA::TrkLst lstT;
     			lstT.push_back(muon);
      			lstT.push_back(d0);
      			lstT.push_back(dpion);
                        AA::Vrt *b0_vrt = _boxv.newVrt();
      			if(!b0_vrt->fill(b->decayVertex()->x(),&lstT)) continue;
      			if(!b0_vrt->filter()) continue;
      			if(b0_vrt->size() != 3) continue;
//      			if(b0_vrt->chi2() > 20.) continue;
 
                        /* -- D* vertex -- */
                        AA::PtlLst dtrack_list;
                        dtrack_list.push_back(d0);
                        dtrack_list.push_back(dpion);

                        /* -- Construct new particle D* from these 2 tracks -- */
                        AA::Ptl* dstar = _boxp.newPtl();
                        if (!dstar->combine(b0_vrt, &dtrack_list, dpion->q()))
                           continue;

       			//Mass of Dstar
      			double mDst = AA::xMass(b0_vrt,dstar->chain(),dpion->chain(), _d0_finder->getMass() ,AA::mass[AA::PI_MINUS]);

                        /* -- cut on B mass -- */
                        double mdiff = mDst - _d0_finder->getMass();
                        if ( mdiff < MASS_DSmD0_MIN || MASS_DSmD0_MAX < mdiff)
                           continue;

                        dstarFound = true;
                        dstarMass = mDst;
                    }

//                    std::cout << "dstarFound: " << dstarMass << std::endl;
            //          << " Chi: " << jpsi_finder->getIndex()
            //          << " Kstar:" << kstar_finder->getIndex() << std::endl;

		    v_d0_index.push_back(_d0_finder->getIndex());
		    v_mu_index.push_back(_mu_finder->getIndex());

		    v_b.push_back(b);
		    v_mdstar.push_back(dstarMass);
		}
         }
//  std::cout << "Found: " << v_b.size() << std::endl;
  return v_b.size();
}

void BMuD0XFinder::begin(){
	index = -1;
}

bool BMuD0XFinder::next(){
	if (++index < v_b.size())
		return true;
	else
		index=-1;
	return false;
}

double BMuD0XFinder::getDstarMass(){
	return v_mdstar[index];
}

void BMuD0XFinder::fill(){
//    std::cout << "Filling: " << index << " D0:" << v_d0_index[index] << " mu:" << v_mu_index[index] << std::endl;
//    std::cout << "Filling: " << index << " B:" << v_b[index] << " mu:" << v_mu_index[index] << std::endl;
	_d0_finder->setIndex(v_d0_index[index]);
	_mu_finder->setIndex(v_mu_index[index]);

	_d0_finder->fill();
	_mu_finder->fill();
	b_saver.fill(*v_b[index]);
	vb_saver.fill(*v_b[index]->decayVertex());  // B vertex info.
	vp_saver.fill(*v_b[index]->primaryVertex());    // Primary vertex info.

	/* -- Info missing by the savers -- */
	dstar_mass = v_mdstar[index];
	/* -- lifetime & lifetime_error -- */
	v_b[index]->decayLengthXY(v_b[index]->primaryVertex(),b_dl,b_vdl);
	_d0_finder->getD0().decayLengthXY(v_b[index]->decayVertex(),d0_dl,d0_vdl);

        AA::PtlLst particle_list;
        particle_list.push_back(&_d0_finder->getKaon());
        particle_list.push_back(&_d0_finder->getPion());
        particle_list.push_back(&_mu_finder->getPtl());
 
        double ctau, vctau;
 	v_b[index]->decayLengthProper(PDG_B_MASS, ctau, vctau, &particle_list);
	b_pdl = ctau; b_epdl=sqrt(fabs(vctau));

        tag_saver.fill(v_b[index]->mom(), v_b[index]->primaryVertex(), &particle_list);  // Tagging info.
}

