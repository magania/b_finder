/*
 * UpsilonFinder.cpp
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#include "UpsilonFinder.h"

UpsilonFinder::UpsilonFinder(TTree &tree) :
	mu_plus_saver("mu_plus",tree, true, true),
	mu_minus_saver("mu_minus",tree, true, true),
	upsilon_saver("upsilon",tree, false, false),
	vrt_saver("upsilon_vrt",tree)
{
	index = -1;
	tree.Branch("upsilon_index", &index, "upsilon_index/D");
	tree.Branch("upsilon_mass", &masa, "upsilon_mass/D");
}

UpsilonFinder::~UpsilonFinder() {
	// TODO Auto-generated destructor stub
}

void UpsilonFinder::clean(){
	upsilon.clear();
	vtx.clear();
	mu_plus.clear();
	mu_minus.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int UpsilonFinder::find(){
	clean();

	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();

	for (AA::PtlLst::const_iterator mu1 = ptl_lst->begin(); mu1 != ptl_lst->end(); ++mu1)
		for (AA::PtlLst::const_iterator mu2 = mu1; mu2 != ptl_lst->end(); ++mu2){
			if ((*mu1) == (*mu2))
				continue;

			if ((*mu1)->primaryVertex() != (*mu2)->primaryVertex())
				continue;

			AA::MuoGlb *muon1 = (*mu1)->muon();
			AA::MuoGlb *muon2 = (*mu2)->muon();

			if ( !muon1 || !muon2 || muon1->nseg()!=3 || muon2->nseg()!=3 )
				continue;

			int qtot = ((*mu1)->q()) * ((*mu2)->q());
			if (qtot > 0)
				continue;

			if ((*mu1)->pt() < MU_MINIMUM_PT || (*mu2)->pt() < MU_MINIMUM_PT)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt vrt;
			if (!AA::findVrt((*mu1), (*mu2), &vrt))
				continue;

			//Construct new particle from these 2 tracmus
			AA::Ptl* ptl = _boxp.newPtl();
			AA::Vrt* pve = _boxv.newVrt(vrt);
			if (!ptl->combine(pve, 0)) {
				_boxp.erase(ptl);
				_boxv.erase(pve);
				continue;
			}

			double masa = ptl->mass(AA::MU_PLUS, AA::MU_MINUS);

			if (MASS_UPSILON_MIN > masa || masa > MASS_UPSILON_MAX)
				continue;

			upsilon.push_back(ptl);
			vtx.push_back(pve);
			mass.push_back(masa);
			if ((*mu1)->q() > 0){
				mu_plus.push_back(*mu1);
				mu_minus.push_back(*mu2);
			} else {
				mu_plus.push_back(*mu2);
				mu_minus.push_back(*mu1);
			}
		}

	return upsilon.size();
}

void UpsilonFinder::begin(){
	index = -1;
}

bool UpsilonFinder::next(){
	if (++index < upsilon.size())
		return true;
	else
		index=-1;
	return false;
}

int UpsilonFinder::getIndex(){
	return index;
}

void UpsilonFinder::setIndex(int i){
	if ( i > upsilon.size()  || i < 0){
		std::cout << "upsilonFinder: Error seting index " << i << " max: " << upsilon.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& UpsilonFinder::getMuPlus(){
	return *mu_plus[index];
}

AA::Ptl& UpsilonFinder::getMuMinus(){
	return *mu_minus[index];
}

AA::Ptl& UpsilonFinder::getUpsilon(){
	return *upsilon[index];
}

AA::Vrt& UpsilonFinder::getVrt(){
	return *vtx[index];
}

double UpsilonFinder::getMass(){
	return mass[index];
}

void UpsilonFinder::fill(){
	mu_plus_saver.fill(getMuPlus());
	mu_minus_saver.fill(getMuMinus());
	upsilon_saver.fill(getUpsilon());
	vrt_saver.fill(getVrt());
	masa = getMass();

}
