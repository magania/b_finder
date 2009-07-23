/*
 * KstarFinder.cpp
 *
 *  Created on: Jul 16, 2009
 *      Author: magania
 */

#include "KstarFinder.h"

KstarFinder::KstarFinder(TTree &tree) :
	kaon_saver("kaon",tree, true, false),
	pion_saver("pion",tree, true, false),
	kstar_saver("kstar",tree, false, false),
	vrt_saver("kstar_vrt",tree)
{
	index = -1;
	tree.Branch("kstar_index", &index, "kstar_index/D");
	tree.Branch("kstar_mass", &masa, "kstar_mass/D");

	MINIMUM_PT = KAON_MINIMUM_PT;
	if (PION_MINIMUM_PT < MINIMUM_PT)
		MINIMUM_PT = PION_MINIMUM_PT;
}

KstarFinder::~KstarFinder() {
}

void KstarFinder::clean(){
	kstar.clear();
	vtx.clear();
	kaon.clear();
	pion.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int KstarFinder::find(){
	clean();

	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();

	for (AA::PtlLst::const_iterator p1 = ptl_lst->begin(); p1 != ptl_lst->end(); ++p1)
		for (AA::PtlLst::const_iterator p2 = p1; p2 != ptl_lst->end(); ++p2){
			if ((*p1) == (*p2))
				continue;

			int qtot = ((*p1)->q()) * ((*p2)->q());

			if (qtot > 0)
				continue;

			if ((*p1)->pt() < MINIMUM_PT || (*p2)->pt() < MINIMUM_PT)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt vrt;
			if (!AA::findVrt((*p1), (*p2), &vrt))
				continue;

			//Construct new particle from these 2 tracks
			AA::Ptl* ptl = _boxp.newPtl();
			AA::Vrt* pve = _boxv.newVrt(vrt);
			if (!ptl->combine(pve, 0)) {
				_boxp.erase(ptl);
				_boxv.erase(pve);
				continue;
			}

			double masa_kp = ptl->mass(AA::K_PLUS, AA::PI_MINUS);
			double masa_pk = ptl->mass(AA::PI_PLUS, AA::K_MINUS);

//std::cout << MASS_KSTAR_MIN << '-' << MASS_KSTAR_MAX << " Mass (kp,pk):" << masa_kp << ' ' << masa_pk << std::endl;
//std::cout << "P1 pt " <<(*p1)->pt() << "  P2 pt " << (*p2)->pt() << std::endl;
			if (MASS_KSTAR_MIN < masa_kp && masa_kp < MASS_KSTAR_MAX &&
					(*p1)->pt() > KAON_MINIMUM_PT && (*p2)->pt() > PION_MINIMUM_PT){
//				std::cout << "Kaon Pion" << std::endl;
				kstar.push_back(ptl);
				vtx.push_back(pve);
				mass.push_back(masa_kp);
				kaon.push_back(*p1);
				pion.push_back(*p2);
			}

			if (MASS_KSTAR_MIN < masa_pk && masa_pk < MASS_KSTAR_MAX &&
					(*p1)->pt() > PION_MINIMUM_PT && (*p2)->pt() > KAON_MINIMUM_PT){
//  			       std::cout << "Pion Kaon" << std::endl;
				kstar.push_back(ptl);
				vtx.push_back(pve);
				mass.push_back(masa_pk);
				pion.push_back(*p1);
				kaon.push_back(*p2);
			}
		}
//std::cout << "Kstar: " << kstar.size() << std::endl;
	return kstar.size();
}

void KstarFinder::begin(){
	index = -1;
}

bool KstarFinder::next(){
	if (++index < kstar.size())
		return true;
	else
		index=-1;
	return false;
}

int KstarFinder::getIndex(){
	return index;
}

void KstarFinder::setIndex(int i){
	if ( i > kstar.size()  || i < 0){
		std::cout << "kstarFinder: Error seting index " << i << " max: " << kstar.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& KstarFinder::getKaon(){
	return *kaon[index];
}

AA::Ptl& KstarFinder::getPion(){
	return *pion[index];
}

AA::Ptl& KstarFinder::getKstar(){
	return *kstar[index];
}

AA::Vrt& KstarFinder::getVrt(){
	return *vtx[index];
}

double KstarFinder::getMass(){
	return mass[index];
}

void KstarFinder::fill(){
	kaon_saver.fill(getKaon());
	pion_saver.fill(getPion());
	kstar_saver.fill(getKstar());
	vrt_saver.fill(getVrt());
	masa = getMass();
}
