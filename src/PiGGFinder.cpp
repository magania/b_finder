/*
 * PiGGFinder.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: magania
 */

#include "PiGGFinder.h"

PiGGFinder::PiGGFinder(TTree &tree, GammaFinder& gamma_finder) :
	pi_saver("pi",tree, false, false),
	vrt_saver("pi_vrt",tree)
{
	_gamma_finder = &gamma_finder;
	index = -1;
	tree.Branch("pi_index", &index, "pi_index/D");
	tree.Branch("pi_mass", &masa, "pi_mass/D");
}

PiGGFinder::~PiGGFinder() {
}

void PiGGFinder::clean(){
	Pi.clear();
	vtx.clear();
	gamma1.clear();
	gamma2.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int PiGGFinder::find(){
	clean();
	_gamma_finder->begin();

	while (_gamma_finder->next()){
		Ptl *gamma1_ptl = &_gamma_finder->getGamma();
		Ptl *e1_plus_ptl = &_gamma_finder->getEPlus();
		Ptl *e1_minus_ptl = &_gamma_finder->getEMinus();
		int gamma1_index = _gamma_finder->getIndex();
		while (_gamma_finder->next()){
			Ptl *gamma2_ptl = &_gamma_finder->getGamma();
			Ptl *e2_plus_ptl = &_gamma_finder->getEPlus();
			Ptl *e2_minus_ptl = &_gamma_finder->getEMinus();
			int gamma2_index = _gamma_finder->getIndex();

		    if (!gamma1_ptl->associateToVrt(&AA::vrtPBox))
		    	continue;
		    if (!gamma2_ptl->associateToVrt(&AA::vrtPBox))
		    	continue;
		    if (gamma1_ptl->primaryVertex() != gamma2_ptl->primaryVertex())
		    	continue;

		    /* -- Pi vertex -- */
		    TrkLst track_list;
		    track_list.push_back(gamma1_ptl);
		    track_list.push_back(gamma2_ptl);
		    Vrt *pi_vrt = _boxv.newVrt();
		    if (!pi_vrt->fill(gamma1_ptl->primaryVertex()->x(), &track_list))
		      continue;
		    if (!pi_vrt->filter())
		      continue;
		    if (pi_vrt->size() != 2)
		      continue;

		    /* -- Construct new particle (X) from these 4 tracks -- */
		    Ptl* pi_ptl = _boxp.newPtl();
		    if (!pi_ptl->combine(pi_vrt, 0))
		      continue;

		    /* -- mass of X and uncertainty --*/
		    vector<PType> types(2);
		    types[0] = AA::GAMMA_STAR;
		    types[2] = AA::GAMMA_STAR;
		    double mpi, vmpi, empi;
		    pi_ptl->mass(types, mpi, vmpi);
		    empi = sqrt(fabs(vmpi));

			Pi.push_back(pi_ptl);
			vtx.push_back(pi_vrt);
			mass.push_back(mpi);
			gamma1.push_back(gamma1_index);
			gamma2.push_back(gamma2_index);
		}
		_gamma_finder->setIndex(gamma1_index);
	}
	return Pi.size();
}

void PiGGFinder::begin(){
	index = -1;
}

bool PiGGFinder::next(){
	if (++index < Pi.size())
		return true;
	else
		index=-1;
	return false;
}

int PiGGFinder::getIndex(){
	return index;
}

void PiGGFinder::setIndex(int i){
	if ( i > Pi.size()  || i < 0){
		std::cout << "PiGGFinder: Error seting index " << i << " max: " << Pi.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& PiGGFinder::getGamma1(){
	_gamma_finder->setIndex(gamma1[index]);
	return _gamma_finder->getGamma();
}

AA::Ptl& PiGGFinder::getGamma2(){
	_gamma_finder->setIndex(gamma2[index]);
	return _gamma_finder->getGamma();
}

AA::Ptl& PiGGFinder::getPi(){
	return *Pi[index];
}

AA::Vrt& PiGGFinder::getVrt(){
	return *vtx[index];
}

double PiGGFinder::getMass(){
	return mass[index];
}

void PiGGFinder::fill(){
	_gamma_finder->fill(gamma1[index],gamma2[index]);
	pi_saver.fill(getPi());
	vrt_saver.fill(getVrt());
	masa = getMass();
}
