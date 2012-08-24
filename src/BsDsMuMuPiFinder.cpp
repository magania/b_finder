/*
 * BsDsMuMuPiFinder.cpp
 *
 *  Created on: Nov 12, 2011
 *      Author: magania
 */

#include "BsDsMuMuPiFinder.h"

BsDsMuMuPiFinder::BsDsMuMuPiFinder(TTree &tree, DsFinder &ds, SSMuFinder &ssmu, PtlFinder &pion) :
	b_saver("b",tree, false, false),
	vb_saver("b_vrt", tree),
	vp_saver("b_pv_vrt", tree)
{
	index = -1;

	_ds_finder = &ds;
	_ssmu_finder = &ssmu;
	_pion_finder = &pion;

	tree.Branch("b_mass", &b_mass, "b_mass/D");
	tree.Branch("b_mass_error", &b_mass_error, "b_mass_error/D");
	tree.Branch("b_dl", &b_dl, "b_dl/D");
}

BsDsMuMuPiFinder::~BsDsMuMuPiFinder() {
}

void BsDsMuMuPiFinder::clean(){
	_boxp.clear();
	_boxv.clear();

        v_ds_index.clear();
	v_ssmu_index.clear();
	v_pion_index.clear();
        v_b.clear();
	v_mb.clear();
	v_emb.clear();
}

int BsDsMuMuPiFinder::find(){
    //std::cout << "Finding B .. " << std::endl;
	clean();
	_ds_finder->begin();
	_ssmu_finder->begin();
	_pion_finder->begin();
	while (_ds_finder->next())
		while (_ssmu_finder->next())
			while (_pion_finder->next()) {
				AA::Ptl* ds = &_ds_finder->getDs();
				AA::Ptl* ds_kplus = &_ds_finder->getKPlus();
				AA::Ptl* ds_kminus = &_ds_finder->getKMinus();
				AA::Ptl* ds_pion = &_ds_finder->getPion();
				AA::Ptl* mu1 = &_ssmu_finder->getMu1();
				AA::Ptl* mu2 = &_ssmu_finder->getMu2();
				AA::Ptl* pion = &_pion_finder->getPtl();
			
				/* ---------------- Selection ------------- */
				if (ds_kplus == mu1 || ds_kplus == mu2 || ds_kplus == pion ||
				    ds_kminus == mu1 || ds_kminus == mu2 || ds_kminus == pion ||
				    ds_pion == mu1 || ds_pion == mu2 || ds_pion == pion ||
				    mu1 == pion || mu2 == pion)
					continue;

				if ( ((ds->q()) * (pion->q())) < 0 )
					continue;
				if ( ((pion->q()) * (mu1->q())) > 0 )
					continue;

				if (pion->pt() < PT_PION_MIN)
					continue;

	    			/* -- b vertex 	-- */
	    			AA::TrkLst track_list;
	    			track_list.push_back(ds);
	    			track_list.push_back(mu1);
	    			track_list.push_back(mu2);
	    			track_list.push_back(pion);
	    			AA::Vrt *b_vrt = _boxv.newVrt();
	    			if (!b_vrt->fill(_ssmu_finder->getVrt().x(), &track_list))
				{
					_boxv.erase(b_vrt);
	      				continue;
				}
	    			if (!b_vrt->filter())
				{
					_boxv.erase(b_vrt);
	      				continue;
				}
	    			if (b_vrt->size() != 4)
				{
					_boxv.erase(b_vrt);
	      				continue;
				}

//				std::cout << "V chi2 = " << b_vrt->chi2() << std::endl;
//				if (b_vrt->chi2() > 10 )
//					continue;

	   			/* -- Construct new particle (B) from these 4 tracks -- */
	    			AA::Ptl* b = _boxp.newPtl();
	    			if (!b->combine(b_vrt, 0))
				{
					_boxp.erase(b);
					_boxv.erase(b_vrt);
	      				continue;
				}
	    			if (!b->associateToVrt(&AA::vrtPBox))
				{
					_boxp.erase(b);
					_boxv.erase(b_vrt);
	      				continue;
				}
				
	    			/* -- mass of B and uncertainty --*/
	    			std::vector<AA::PType> types(4);
	    			types[0] = AA::DS_PLUS;	
	    			types[1] = AA::MU_MINUS;
	    			types[2] = AA::MU_MINUS;
	    			types[3] = AA::PI_PLUS;
	    			double mb, vmb, emb;
	    			b->mass(types, mb, vmb);
	    			emb = sqrt(fabs(vmb));

	    			/* -- cut on B mass -- */
	    			if ( mb < MASS_B_MIN || MASS_B_MAX < mb)
				{
					_boxp.erase(b);
					_boxv.erase(b_vrt);
	      				continue;
				}
//				std::cout << "Here 3 " << std::endl;

				if ( blind ) {
		    			if ( mb < MASS_B_RIGHT && MASS_B_LEFT < mb)
					{
						_boxp.erase(b);
						_boxv.erase(b_vrt);
	      					continue;
					}
				}

	    			v_ds_index.push_back(_ds_finder->getIndex());
	    			v_ssmu_index.push_back(_ssmu_finder->getIndex());
	    			v_pion_index.push_back(_pion_finder->getIndex());

	    			v_b.push_back(b);
	    			v_mb.push_back(mb);
	    			v_emb.push_back(emb);
				}
//  std::cout << "Found: " << v_b.size() << std::endl;
  return v_b.size();
}

void BsDsMuMuPiFinder::begin(){
	index = -1;
}

bool BsDsMuMuPiFinder::next(){
	if (++index < v_b.size())
		return true;
	else
		index=-1;
	return false;
}

double BsDsMuMuPiFinder::getMass(){
	return v_mb[index];
}

void BsDsMuMuPiFinder::fill(){
    //std::cout << "Filling: " << index << " Chi:" << v_jpsi_index[index] << " Kstar:" << v_kstar_index[index] << std::endl;
	_ds_finder->setIndex(v_ds_index[index]);
	_ssmu_finder->setIndex(v_ssmu_index[index]);
	_pion_finder->setIndex(v_pion_index[index]);

	_ds_finder->fill();
	_ssmu_finder->fill();
	_pion_finder->fill();
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
