/*
 * DsFinder.cpp
 *
 *  Created on: Nov 12, 2011
 *      Author: magania
 */

#include "DsFinder.h"

DsFinder::DsFinder(TTree &tree, PhiFinder &phi, PtlFinder &pion) :
	ds_saver("ds",tree, false, false)
{
	index = -1;

	_phi_finder = &phi;
	_pion_finder = &pion;

	tree.Branch("ds_mass", &ds_mass, "ds_mass/D");
	tree.Branch("ds_mass_error", &ds_mass_error, "ds_mass_error/D");
}

DsFinder::~DsFinder() {
}

void DsFinder::clean(){
	_boxp.clear();
	_boxv.clear();

	v_phi_index.clear();
	v_pion_index.clear();
        v_ds.clear();
	v_mds.clear();
	v_emds.clear();
}

int DsFinder::find(){
        //std::cout << "Finding Ds .. " << std::endl;
	clean();
	_phi_finder->begin();
	_pion_finder->begin();
	while (_phi_finder->next())
		while (_pion_finder->next()) {
			AA::Ptl* phi = &_phi_finder->getPhi();
			AA::Ptl* pion = &_pion_finder->getPtl();
			AA::Ptl* k_plus = &_phi_finder->getKPlus();
			AA::Ptl* k_minus = &_phi_finder->getKMinus();

			/* ---------------- Selection ------------- */
			if (pion == k_plus || pion == k_minus)
				continue;

			if (pion->pt() < PT_PION_MIN)
				continue;

		    /* -- Ds vertex -- */
		    AA::TrkLst track_list;
		    track_list.push_back(k_plus);
		    track_list.push_back(k_minus);
		    track_list.push_back(pion);
		    AA::Vrt *ds_vrt = _boxv.newVrt();
		    if (!ds_vrt->fill(phi->decayVertex()->x(), &track_list))
		      continue;
		    if (!ds_vrt->filter())
		      continue;
		    if (ds_vrt->size() != 3)
		      continue;

		    /* -- Construct new particle (Ds) from these 3 tracks -- */
		    AA::Ptl* ds = _boxp.newPtl();
		    if (!ds->combine(ds_vrt, pion->q()))
		      continue;

		    /* -- mass of Ds and uncertainty --*/
		    std::vector<AA::PType> types(3);
		    types[0] = AA::K_PLUS;
		    types[1] = AA::K_MINUS;
		    types[2] = AA::PI_PLUS;
		    double mds, vmds, emds;
		    ds->mass(types, mds, vmds);
		    emds = sqrt(fabs(vmds));

		    /* -- cut on B mass -- */
		    if ( mds < MASS_DS_MIN || MASS_DS_MAX < mds)
		  	  continue;

            //std::cout << "Bs: " << v_bd.size()
            //          << " Chi: " << jpsi_finder->getIndex()
            //          << " Kstar:" << kstar_finder->getIndex() << std::endl;
		    v_phi_index.push_back(_phi_finder->getIndex());
		    v_pion_index.push_back(_pion_finder->getIndex());

		    v_ds.push_back(ds);
		    v_mds.push_back(mds);
		    v_emds.push_back(emds);
		}
  //std::cout << "Found: " << v_bd.size() << std::endl;
  return v_ds.size();
}

void DsFinder::begin(){
	index = -1;
}

int DsFinder::getIndex(){
        return index;
}

void DsFinder::setIndex(int i){
	if ( i > v_ds.size()  || i < 0){
		std::cout << "DsFinder: Error seting index " << i << " max: " << v_ds.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

bool DsFinder::next(){
	if (++index < v_ds.size())
		return true;
	else
		index=-1;
	return false;
}

double DsFinder::getMass(){
	return v_mds[index];
}

AA::Ptl& DsFinder::getDs(){
  	return *v_ds[index];
}

AA::Ptl& DsFinder::getKPlus(){
	_phi_finder->setIndex(v_phi_index[index]);
	return _phi_finder->getKPlus();;
}

AA::Ptl& DsFinder::getKMinus(){
	_phi_finder->setIndex(v_phi_index[index]);
	return _phi_finder->getKMinus();;
}

AA::Ptl& DsFinder::getPion(){
	_pion_finder->setIndex(v_pion_index[index]);
	return _pion_finder->getPtl();;
}

void DsFinder::fill(){
    //std::cout << "Filling: " << index << " Chi:" << v_jpsi_index[index] << " Kstar:" << v_kstar_index[index] << std::endl;
	_phi_finder->setIndex(v_phi_index[index]);
	_pion_finder->setIndex(v_pion_index[index]);

	_phi_finder->fill();
	_pion_finder->fill();
	ds_saver.fill(*v_ds[index]);  // Bs info.

	/* -- Info missing by the savers -- */
	ds_mass = v_mds[index];
	ds_mass_error = v_emds[index];
}

