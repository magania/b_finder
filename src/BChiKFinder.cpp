/*
 * BChiKFinder.cpp
 *
 *  Created on: Sep 17, 2009
 *      Author: magania
 */

#include "BChiKFinder.h"

BChiKFinder::BChiKFinder(TTree &tree, ChiJPsiGFinder &chi, PtlFinder &kaon) :
	b_saver("b",tree, false, false),
	vb_saver("b_vrt", tree),
	vp_saver("b_pv_vrt", tree)
{
	index = -1;

	_chi_finder = &chi;
	_kaon_finder = &kaon;

	tree.Branch("b_mass", &b_mass, "b_mass/D");
	tree.Branch("b_mass_error", &b_mass_error, "b_mass_error/D");
	tree.Branch("b_dl", &b_dl, "b_dl/D");
}

BChiKFinder::~BChiKFinder() {
}

void BChiKFinder::clean(){
	_boxp.clear();
	_boxv.clear();

    v_chi_index.clear();
	v_kaon_index.clear();
    v_b.clear();
	v_mb.clear();
	v_emb.clear();
}

int BChiKFinder::find(){
    //std::cout << "Finding B .. " << std::endl;
	clean();
	_chi_finder->begin();
	_kaon_finder->begin();
	while (_chi_finder->next())
		while (_kaon_finder->next()) {
			Ptl* chi = &_chi_finder->getChi();
			Ptl* jpsi = &_chi_finder->getJPsi();
			Ptl* mu_plus = &_chi_finder->getMuPlus();
			Ptl* mu_minus = &_chi_finder->getJPsi();
			Ptl* gamma = &_chi_finder->getGamma();
			Ptl* e_plus = &_chi_finder->getEPlus();
			Ptl* e_minus = &_chi_finder->getEMinus();
			Ptl* kaon = &_kaon_finder->getPtl();

			/* ---------------- Selection ------------- */
			if (kaon == mu_plus || kaon == mu_minus || kaon == e_plus || kaon == e_minus)
				continue;

			if (kaon->pt() < PT_KAON_MIN)
				continue;

		    /* -- b vertex -- */
		    TrkLst track_list;
		    track_list.push_back(mu_plus);
		    track_list.push_back(mu_minus);
		    track_list.push_back(gamma);
		    track_list.push_back(kaon);
		    Vrt *b_vrt = _boxv.newVrt();
		    if (!b_vrt->fill(jpsi->decayVertex()->x(), &track_list))
		      continue;
		    if (!b_vrt->filter())
		      continue;
		    if (b_vrt->size() != 4)
		      continue;

		    /* -- Construct new particle (B) from these 4 tracks -- */
		    Ptl* b = _boxp.newPtl();
		    if (!b->combine(b_vrt, 0))
		      continue;
		    if (!b->associateToVrt(&AA::vrtPBox))
		      continue;

		    /* -- mass of B and uncertainty --*/
		    vector<PType> types(4);
		    types[0] = AA::MU_PLUS;
		    types[1] = AA::MU_MINUS;
		    types[2] = AA::GAMMA_STAR;
		    types[3] = AA::K_PLUS;
		    double mb, vmb, emb;
		    b->mass(types, mb, vmb);
		    emb = sqrt(fabs(vmb));

		    /* -- cut on B mass -- */
		    if ( mb < MASS_B_MIN || MASS_B_MAX < mb)
		  	  continue;

            //std::cout << "Bs: " << v_bd.size()
            //          << " Chi: " << jpsi_finder->getIndex()
            //          << " Kstar:" << kstar_finder->getIndex() << std::endl;
		    v_chi_index.push_back(_chi_finder->getIndex());
		    v_kaon_index.push_back(_kaon_finder->getIndex());

		    v_b.push_back(b);
		    v_mb.push_back(mb);
		    v_emb.push_back(emb);
		}
  //std::cout << "Found: " << v_bd.size() << std::endl;
  return v_b.size();
}

void BChiKFinder::begin(){
	index = -1;
}

bool BChiKFinder::next(){
	if (++index < v_b.size())
		return true;
	else
		index=-1;
	return false;
}

double BChiKFinder::getMass(){
	return v_mb[index];
}

void BChiKFinder::fill(){
    //std::cout << "Filling: " << index << " Chi:" << v_jpsi_index[index] << " Kstar:" << v_kstar_index[index] << std::endl;
	_chi_finder->setIndex(v_chi_index[index]);
	_kaon_finder->setIndex(v_kaon_index[index]);

	_chi_finder->fill();
	_kaon_finder->fill();
	b_saver.fill(*v_b[index]);  // Bs info.
	vb_saver.fill(*v_b[index]->decayVertex());  // B vertex info.
	vp_saver.fill(*v_b[index]->primaryVertex());    // Primary vertex info.

	/* -- Info missing by the savers -- */
	b_mass = v_mb[index];
	b_mass_error = v_emb[index];

	/* -- lifetime & lifetime_error -- */
    double lxy, vlxy;
	v_b[index]->decayLengthXY(v_b[index]->primaryVertex(),lxy,vlxy);
	b_dl = lxy;
}

