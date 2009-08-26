/*
 * XYGammaFinder.cpp
 *
 *  Created on: Aug 13, 2009
 *      Author: magania
 */

#include "XYGammaFinder.h"

XYGammaFinder::XYGammaFinder(TTree &tree, UpsilonFinder& upsilon_finder, GammaFinder& gamma_finder) :
	x_saver("b",tree, false, false),
	vrt_saver("b_vrt",tree),
	pv_saver("b_pv",tree)
{
	_upsilon_finder = &upsilon_finder;
	_gamma_finder = &gamma_finder;
	index = -1;
	tree.Branch("b_index", &index, "b_index/D");
	tree.Branch("b_mass", &masa, "b_mass/D");
}

XYGammaFinder::~XYGammaFinder() {
}

void XYGammaFinder::clean(){
	X.clear();
	vtx.clear();
	upsilon.clear();
	gamma.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int XYGammaFinder::find(){
	clean();
	_upsilon_finder->begin();
	_gamma_finder->begin();

	while (_upsilon_finder->next())
		while (_gamma_finder->next()){
			Ptl *upsilon_ptl = &_upsilon_finder->getUpsilon();
			Ptl *mu_plus_ptl = &_upsilon_finder->getMuPlus();
			Ptl *mu_minus_ptl = &_upsilon_finder->getMuMinus();
			Ptl *gamma_ptl = &_gamma_finder->getGamma();
			Ptl *e_plus_ptl = &_gamma_finder->getEPlus();
			Ptl *e_minus_ptl = &_gamma_finder->getEMinus();

		    /* -- x vertex -- */
		    TrkLst track_list;
		    track_list.push_back(mu_plus_ptl);
		    track_list.push_back(mu_minus_ptl);
		    track_list.push_back(gamma_ptl);
		    Vrt *x_vrt = _boxv.newVrt();
		    if (!x_vrt->fill(_upsilon_finder->getVrt().x(), &track_list))
		      continue;
		    if (!x_vrt->filter())
		      continue;
		    if (x_vrt->size() != 3)
		      continue;

		    /* -- Construct new particle (X) from these 4 tracks -- */
		    Ptl* x_ptl = _boxp.newPtl();
		    if (!x_ptl->combine(x_vrt, 0))
		      continue;

		    /* -- mass of X and uncertainty --*/
		    vector<PType> types(3);
		    types[0] = AA::MU_PLUS;
		    types[1] = AA::MU_MINUS;
		    types[2] = AA::GAMMA_STAR;
		    double mx, vmx, emx;
		    x_ptl->mass(types, mx, vmx);
		    emx = sqrt(fabs(vmx));

		    /* -- cut on B mass -- */
		    if ( mx < MASS_X_MIN || MASS_X_MAX < mx)
		  	  continue;

			X.push_back(x_ptl);
			vtx.push_back(x_vrt);
			mass.push_back(mx);
			upsilon.push_back(_upsilon_finder->getIndex());
			gamma.push_back(_gamma_finder->getIndex());
		}
	return X.size();
}

void XYGammaFinder::begin(){
	index = -1;
}

bool XYGammaFinder::next(){
	if (++index < X.size())
		return true;
	else
		index=-1;
	return false;
}

int XYGammaFinder::getIndex(){
	return index;
}

void XYGammaFinder::setIndex(int i){
	if ( i > X.size()  || i < 0){
		std::cout << "XYGammaFinder: Error seting index " << i << " max: " << X.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& XYGammaFinder::getUpsilon(){
	return _upsilon_finder->getUpsilon();
}

AA::Ptl& XYGammaFinder::getGamma(){
	return _gamma_finder->getGamma();
}

AA::Ptl& XYGammaFinder::getX(){
	return *X[index];
}

AA::Vrt& XYGammaFinder::getVrt(){
	return *vtx[index];
}

double XYGammaFinder::getMass(){
	return mass[index];
}

void XYGammaFinder::fill(){
	_upsilon_finder->setIndex(upsilon[index]);
	_gamma_finder->setIndex(gamma[index]);
	_upsilon_finder->fill();
	_gamma_finder->fill();
	x_saver.fill(getX());
	vrt_saver.fill(getVrt());
	masa = getMass();
}
