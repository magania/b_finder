/*
 * ChiJPsiGFinder.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: magania
 */

#include "ChiJPsiGFinder.h"

ChiJPsiGFinder::ChiJPsiGFinder(TTree &tree, JPsiFinder& jpsi_finder, GammaFinder& gamma_finder) :
	chi_saver("chi",tree, false, false),
	vrt_saver("chi_vrt",tree)
{
	_jpsi_finder = &jpsi_finder;
	_gamma_finder = &gamma_finder;
	index = -1;
	tree.Branch("chi_index", &index, "chi_index/D");
	tree.Branch("chi_mass", &masa, "chi_mass/D");
}

ChiJPsiGFinder::~ChiJPsiGFinder() {
}

void ChiJPsiGFinder::clean(){
	chi.clear();
	vtx.clear();
	jpsi.clear();
	gamma.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int ChiJPsiGFinder::find(){
	clean();
	_jpsi_finder->begin();
	_gamma_finder->begin();

	while (_jpsi_finder->next())
		while (_gamma_finder->next()){
			Ptl *jpsi_ptl = &_jpsi_finder->getJPsi();
			Ptl *mu_plus_ptl = &_jpsi_finder->getMuPlus();
			Ptl *mu_minus_ptl = &_jpsi_finder->getMuMinus();
			Ptl *gamma_ptl = &_gamma_finder->getGamma();
			Ptl *e_plus_ptl = &_gamma_finder->getEPlus();
			Ptl *e_minus_ptl = &_gamma_finder->getEMinus();

			if(mu_plus_ptl == e_plus_ptl || mu_minus_ptl == e_minus_ptl)
				continue;

		    /* -- chi vertex -- */
		    TrkLst track_list;
		    track_list.push_back(mu_plus_ptl);
		    track_list.push_back(mu_minus_ptl);
		    track_list.push_back(gamma_ptl);
		    Vrt *chi_vrt = _boxv.newVrt();
		    if (!chi_vrt->fill(_jpsi_finder->getJPsi().decayVertex()->x(), &track_list))
		      continue;
		    if (!chi_vrt->filter())
		      continue;
		    if (chi_vrt->size() != 3)
		      continue;

		    /* -- Construct new particle (Chi) from these 4 tracks -- */
		    Ptl* chi_ptl = _boxp.newPtl();
		    if (!chi_ptl->combine(chi_vrt, 0))
		      continue;

		    /* -- mass of Chi and uncertainty --*/
		    vector<PType> types(3);
		    types[0] = AA::MU_PLUS;
		    types[1] = AA::MU_MINUS;
		    types[2] = AA::GAMMA_STAR;
		    double mchi, vmchi, emchi;
		    chi_ptl->mass(types, mchi, vmchi);
		    emchi = sqrt(fabs(vmchi));

		    /* -- cut on B mass -- */
		    if ( mchi < MASS_CHI_MIN || MASS_CHI_MAX < mchi)
		  	  continue;

			chi.push_back(chi_ptl);
			vtx.push_back(chi_vrt);
			mass.push_back(mchi);
			jpsi.push_back(_jpsi_finder->getIndex());
			gamma.push_back(_gamma_finder->getIndex());
		}
	return chi.size();
}

void ChiJPsiGFinder::begin(){
	index = -1;
}

bool ChiJPsiGFinder::next(){
	if (++index < chi.size())
		return true;
	else
		index=-1;
	return false;
}

int ChiJPsiGFinder::getIndex(){
	return index;
}

void ChiJPsiGFinder::setIndex(int i){
	if ( i > chi.size()  || i < 0){
		std::cout << "ChiJPsiGFinder: Error seting index " << i << " max: " << chi.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& ChiJPsiGFinder::getJPsi(){
	_jpsi_finder->setIndex(jpsi[index]);
	return _jpsi_finder->getJPsi();
}

AA::Ptl& ChiJPsiGFinder::getGamma(){
	_gamma_finder->setIndex(gamma[index]);
	return _gamma_finder->getGamma();
}

AA::Ptl& ChiJPsiGFinder::getChi(){
	return *chi[index];
}

AA::Vrt& ChiJPsiGFinder::getVrt(){
	return *vtx[index];
}

double ChiJPsiGFinder::getMass(){
	return mass[index];
}

void ChiJPsiGFinder::fill(){
	_jpsi_finder->setIndex(jpsi[index]);
	_gamma_finder->setIndex(gamma[index]);
	_jpsi_finder->fill();
	_gamma_finder->fill();
	chi_saver.fill(getChi());
	vrt_saver.fill(getVrt());
	masa = getMass();
}
