/*
 * PhiFinder.cpp
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#include "PhiFinder.h"

PhiFinder::PhiFinder(TTree &tree) :
	k_plus_saver("k_plus",tree, true, false),
	k_minus_saver("k_minus",tree, true, false),
	phi_saver("phi",tree, false, false),
	vrt_saver("phi_vrt",tree)
{
	index = -1;
	tree.Branch("phi_index", &index, "phi_index/D");
	tree.Branch("phi_mass", &masa, "phi_mass/D");

	tree.Branch("k_leading_pt", &k_leading_pt, "k_leading_pt/D");
	tree.Branch("k_trailing_pt", &k_trailing_pt, "k_trailing_pt/D");
	tree.Branch("k_leading_ptot", &k_leading_ptot, "k_leading_ptot/D");
	tree.Branch("k_trailing_ptot", &k_trailing_ptot, "k_trailing_ptot/D");
}

PhiFinder::~PhiFinder() {
	// TODO Auto-generated destructor stub
}

void PhiFinder::clean(){
	phi.clear();
	vtx.clear();
	k_plus.clear();
	k_minus.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int PhiFinder::find(){
	clean();

	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();

	for (AA::PtlLst::const_iterator k1 = ptl_lst->begin(); k1 != ptl_lst->end(); ++k1)
		for (AA::PtlLst::const_iterator k2 = k1; k2 != ptl_lst->end(); ++k2){
			if ((*k1) == (*k2))
				continue;

			int qtot = ((*k1)->q()) * ((*k2)->q());

			if (qtot > 0)
				continue;

			if ((*k1)->pt() < K_MINIMUM_PT || (*k2)->pt() < K_MINIMUM_PT)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt vrt;
			if (!AA::findVrt((*k1), (*k2), &vrt))
				continue;

			//Construct new particle from these 2 tracks
			AA::Ptl* ptl = _boxp.newPtl();
			AA::Vrt* pve = _boxv.newVrt(vrt);
			if (!ptl->combine(pve, 0)) {
				_boxp.erase(ptl);
				_boxv.erase(pve);
				continue;
			}

			double masa = ptl->mass(AA::K_PLUS, AA::K_MINUS);

			if (MASS_PHI_MIN > masa || masa > MASS_PHI_MAX)
				continue;

			phi.push_back(ptl);
			vtx.push_back(pve);
			mass.push_back(masa);
			if ((*k1)->q() > 0){
				k_plus.push_back(*k1);
				k_minus.push_back(*k2);
			} else {
				k_plus.push_back(*k2);
				k_minus.push_back(*k1);
			}
		}

	return phi.size();
}

void PhiFinder::begin(){
	index = -1;
}

bool PhiFinder::next(){
	if (++index < phi.size())
		return true;
	else
		index=-1;
	return false;
}

int PhiFinder::getIndex(){
	return index;
}

void PhiFinder::setIndex(int i){
	if ( i > phi.size()  || i < 0){
		std::cout << "phiFinder: Error seting index " << i << " max: " << phi.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = 1;
	}
}

AA::Ptl& PhiFinder::getKPlus(){
	return *k_plus[index];
}

AA::Ptl& PhiFinder::getKMinus(){
	return *k_minus[index];
}

AA::Ptl& PhiFinder::getPhi(){
	return *phi[index];
}

AA::Vrt& PhiFinder::getVrt(){
	return *vtx[index];
}

double PhiFinder::getMass(){
	return mass[index];
}

void PhiFinder::fill(){
	k_plus_saver.fill(getKPlus());
	k_minus_saver.fill(getKMinus());
	phi_saver.fill(getPhi());
	vrt_saver.fill(getVrt());
	masa = getMass();
	if (getKPlus().pt() > getKMinus().pt()){
		k_leading_pt = getKPlus().pt();
		k_trailing_pt = getKMinus().pt();
	} else {
		k_leading_pt = getKMinus().pt();
		k_trailing_pt = getKPlus().pt();
	}
	if (getKPlus().ptot() > getKMinus().ptot()){
		k_leading_ptot = getKPlus().ptot();
		k_trailing_ptot = getKMinus().ptot();
	} else {
		k_leading_ptot = getKMinus().ptot();
		k_trailing_ptot = getKPlus().ptot();
	}
}
