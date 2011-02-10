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

        tree.Branch("mc_K", &mc_K, "mc_K/D");
        tree.Branch("mc_L", &mc_L, "mc_L/D");
        tree.Branch("mc_ct", &mc_ct, "mc_ct/D");
        tree.Branch("mc_b", &mc_b, "mc_b/I");
        tree.Branch("mc_d", &mc_d, "mc_d/I");
        tree.Branch("mc_dc", &mc_dc, "mc_dc/I");
        tree.Branch("mc_dmatch", &mc_dmatch, "mc_dmatch/I");
        tree.Branch("mc_smatch", &mc_smatch, "mc_smatch/I");
        tree.Branch("mc_mmatch", &mc_mmatch, "mc_mmatch/I");
}

BMuD0XFinder::~BMuD0XFinder() {
}

void BMuD0XFinder::clean(){
	_boxp.clear();
	_boxv.clear();

        v_d0_index.clear();
	v_mu_index.clear();
	v_pi_index.clear();
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

                        if (fabs(d0->cosXY()) < COSXY_D0_MIN)
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

//                        double blxy,bvlxy;
//                        b_vrt->distanceXY(muon->primaryVertex(),blxy,bvlxy);
//                        if (blxy*blxy > d0lxy*d0lxy) continue;

                        double bdlxy,bdvlxy;
                        d0->decayVertex()->distanceXY(b_vrt,bdlxy,bdvlxy);
                        if ( bdlxy*bdlxy > SIGMA_D0_LXY*bdvlxy ) continue;

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

//                    if (fabs(b->cosXY()) < COSXY_B_MIN)
//                         continue;

                    double dstarMass = -1.0;
                    int piFinderIndex = 0;
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

			if (dstarMass != -1.0 && TMath::Abs(mDst - 2010.25) > TMath::Abs(dstarMass - 2010.25))
    			   continue;
	
			piFinderIndex = _pi_finder->getIndex();
                        dstarMass = mDst;
                    }

//                    std::cout << "dstarFound: " << dstarMass << std::endl;
            //          << " Chi: " << jpsi_finder->getIndex()
            //          << " Kstar:" << kstar_finder->getIndex() << std::endl;

		    v_d0_index.push_back(_d0_finder->getIndex());
		    v_mu_index.push_back(_mu_finder->getIndex());
		    v_pi_index.push_back(piFinderIndex);

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
	_pi_finder->setIndex(v_pi_index[index]);

	_d0_finder->fill();
	_mu_finder->fill();
	_pi_finder->fill();
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
 	v_b[index]->decayLengthProper(PDG_BP_MASS, ctau, vctau, &particle_list);
	b_pdl = ctau; b_epdl=sqrt(fabs(vctau));
//	cout << "L = " << b_pdl << endl;

        tag_saver.fill(v_b[index]->mom(), v_b[index]->primaryVertex(), &particle_list);  // Tagging info.

//	cout << "MC Fill" << endl;
	mc_dc = -1;
	mc_b = -1;
	mc_d = -1;
	mc_L = -10;
	mc_ct = -10;

#ifdef MC
	cout << "MC Fill" << endl;
	Ptl* muon = &_mu_finder->getPtl();
	Ptl* kaon = &_d0_finder->getKaon();
	Ptl* pion = &_d0_finder->getPion();
	Ptl* dpion = &_pi_finder->getPtl();

//	cout << muon << " " << kaon << " " << pion << endl;

        PtlMC* muonMC, *kaonMC, *pionMC, *dpionMC;
	double chi;

	muonMC = DecayMC::imatch(muon,chi);
	kaonMC = DecayMC::imatch(kaon,chi);
	pionMC = DecayMC::imatch(pion,chi);
	dpionMC = DecayMC::imatch(dpion,chi);

//	cout << muonMC << " " << kaonMC << " " << pionMC << endl;
	if ( !muonMC )
		return;	

        PtlMC* bMC = DecayMC::getParent(muonMC);
	mc_b = bMC->idPDG();
	cout << "** ";
	printChain(bMC);
	cout << endl;

	if (bMC && kaonMC && pionMC && muonMC)
	{
	        double px = kaonMC->mom4(1) + pionMC->mom4(1) + muonMC->mom4(1);
	        double py = kaonMC->mom4(2) + pionMC->mom4(2) + muonMC->mom4(2);
	        double ptMuD0 = TMath::Sqrt(px*px+py*py);

	        mc_K = ptMuD0/bMC->pt();
//		cout << "ptB = " << bMC->pt() << "   ptMuD0 = " << ptMuD0 << "   K = " << mc_K << endl;
		VrtMC *vStart = bMC->vrtStart();
		VrtMC *vEnd = bMC->vrtEnd();

//		cout << vStart->x() << endl << vEnd->x()<< endl;
		double lx = vEnd->x(1) - vStart->x(1);
		double ly = vEnd->x(2) - vStart->x(2);
		
		if ( bMC->idPDG() == 521){
			mc_L  = PDG_BP_MASS*(lx*px+ly*py)/ptMuD0/ptMuD0;
			mc_ct = PDG_BP_MASS*(lx*bMC->mom4(1)+ly*bMC->mom4(2))/bMC->pt()/bMC->pt();
		}
		if ( bMC->idPDG() == 511){
			mc_L  = PDG_B0_MASS*(lx*px+ly*py)/ptMuD0/ptMuD0;
			mc_ct = PDG_B0_MASS*(lx*bMC->mom4(1)+ly*bMC->mom4(2))/bMC->pt()/bMC->pt();
		}
		if ( bMC->idPDG() == 531){
			mc_L  = PDG_BS_MASS*(lx*px+ly*py)/ptMuD0/ptMuD0;
			mc_ct = PDG_BS_MASS*(lx*bMC->mom4(1)+ly*bMC->mom4(2))/bMC->pt()/bMC->pt();
		}

//		cout << "L = (" << lx <<","<< lx << ")  mc_L = " << mc_L << endl;
	} else {
		mc_K = -1.0;
//		cout << "K = " << -1.0 << endl;
	}

	mc_dmatch = matchD0(bMC,kaonMC,pionMC);
	mc_smatch = matchD0Star(bMC,kaonMC,pionMC);
	mc_mmatch = matchD0Sm(bMC,kaonMC,pionMC,dpionMC);


	mc_dc = 0;
        PtlMCLst child_list = bMC->vrtEnd()->children();
        for (PtlMCLstCIt c = child_list.begin(); c != child_list.end(); ++c) 
	{
		mc_dc++;
		if ( (*c)->idPDG() == -13 )
			continue;
		if ( (*c)->idPDG() == 14 )
			continue;
		if ( matchD0((*c),kaonMC,pionMC) ){
			mc_d = (*c)->idPDG();
			continue;
		}
	}

	cout << "dmatch:" << mc_dmatch << "  smatch:" << mc_smatch << "  mmatch:" << mc_mmatch << "  b:" << mc_b << "  d:" << mc_d << "  dc:" << mc_dc << endl;
#endif
}

bool BMuD0XFinder::matchD0Sm(PtlMC* p, PtlMC* kaonMC, PtlMC* pionMC, PtlMC* dpionMC)
{
        PtlMC *chD0, *chdPion, *chKaon, *chPion;
	if ( p->idPDG() == -413 &&
		 DecayMC::asignChildren(p, -421, &chD0, -211, &chdPion) &&
                 DecayMC::asignChildren(chD0, 321, &chKaon, -211, &chPion) &&
                 chKaon == kaonMC &&
                 chPion == pionMC &&
                 chdPion == dpionMC )
		return true;
        if ( p->vrtEnd() ){
                PtlMCLst child = p->vrtEnd()->children();
                for (PtlMCLstCIt c = child.begin(); c != child.end(); ++c)
                        if ( matchD0Sm((*c),kaonMC,pionMC,dpionMC) )
                                return true;
        }
        return false;
}

bool BMuD0XFinder::matchD0Star(PtlMC* p, PtlMC* kaonMC, PtlMC* pionMC)
{
        PtlMC *chD0, *chPiOrGamma, *chKaon, *chPion;
	if ( p->idPDG() == -423 && 
		 (DecayMC::asignChildren(p, -421, &chD0, 111, &chPiOrGamma) || DecayMC::asignChildren(p, -421, &chD0, 22, &chPiOrGamma)) &&
	         DecayMC::asignChildren(chD0, 321, &chKaon, -211, &chPion) &&
                 chKaon == kaonMC &&
                 chPion == pionMC )
		return true;
        if ( p->vrtEnd() ){
                PtlMCLst child = p->vrtEnd()->children();
                for (PtlMCLstCIt c = child.begin(); c != child.end(); ++c)
                        if ( matchD0Star((*c),kaonMC,pionMC) )
                                return true;
        }
        return false;
}

bool BMuD0XFinder::matchD0(PtlMC* p, PtlMC* kaonMC, PtlMC* pionMC)
{
        PtlMC *chKaon, *chPion;
	if ( p->idPDG() == -421  &&
		 DecayMC::asignChildren(p, 321, &chKaon, -211, &chPion)  &&
		 chKaon == kaonMC &&
		 chPion == pionMC )
		return true;
	if ( p->vrtEnd() ){
		PtlMCLst child = p->vrtEnd()->children();
	        for (PtlMCLstCIt c = child.begin(); c != child.end(); ++c)
			if ( matchD0((*c),kaonMC,pionMC) )
				return true;
	}
	return false;
}

void BMuD0XFinder::printChain(PtlMC* p)
{
	cout << DecayMC::getIdPdg(p);
	if ( p->vrtEnd() )
	{
		cout << "[";
        	PtlMCLst child = p->vrtEnd()->children();
                for (PtlMCLstCIt c = child.begin(); c != child.end(); ++c)
			printChain(*c);
		cout << "] ";
	} else {
		cout << " ";
	}
	
}

