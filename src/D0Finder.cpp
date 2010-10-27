/*
 * D0Finder.cpp
 *
 *  Created on: Oct 22, 2010
 *      Author: magania
 */

#include "D0Finder.h"

D0Finder::D0Finder(TTree &tree) :
	kaon_saver("kaon",tree, true, false),
	pion_saver("pion",tree, true, false),
	d0_saver("d0",tree, false, false),
	vrt_saver("d0_vrt",tree)
{
	index = -1;
	tree.Branch("d0_index", &index, "d0_index/D");
	tree.Branch("d0_mass", &masa, "d0_mass/D");

	MINIMUM_PT = KAON_MINIMUM_PT;
	if (PION_MINIMUM_PT < MINIMUM_PT)
		MINIMUM_PT = PION_MINIMUM_PT;
}

D0Finder::~D0Finder() {
}

void D0Finder::clean(){
	d0.clear();
	vtx.clear();
	kaon.clear();
	pion.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int D0Finder::find(){
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
 
  		       if ((*p1)->nSMT() < MINIMUM_NSMT || (*p2)->nSMT() < MINIMUM_NSMT)
				continue;

  		       if ((*p1)->nCFT() < MINIMUM_NCFT || (*p2)->nCFT() < MINIMUM_NCFT)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt vrt;
			if (!AA::findVrt((*p1), (*p2), &vrt))
				continue;

                        if (vrt.chi2() < VRT_CHI2_MIN)
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
			if (MASS_D0_MIN < masa_kp && masa_kp < MASS_D0_MAX &&
					(*p1)->pt() > KAON_MINIMUM_PT && (*p2)->pt() > PION_MINIMUM_PT){
//				std::cout << "Kaon Pion" << std::endl;
				d0.push_back(ptl);
				vtx.push_back(pve);
				mass.push_back(masa_kp);
				kaon.push_back(*p1);
				pion.push_back(*p2);
			}

			if (MASS_D0_MIN < masa_pk && masa_pk < MASS_D0_MAX &&
					(*p1)->pt() > PION_MINIMUM_PT && (*p2)->pt() > KAON_MINIMUM_PT){
//  			       std::cout << "Pion Kaon" << std::endl;
				d0.push_back(ptl);
				vtx.push_back(pve);
				mass.push_back(masa_pk);
				pion.push_back(*p1);
				kaon.push_back(*p2);
			}
		}
//std::cout << "Kstar: " << kstar.size() << std::endl;
	return d0.size();
}

void D0Finder::begin(){
	index = -1;
}

bool D0Finder::next(){
	if (++index < d0.size())
		return true;
	else
		index=-1;
	return false;
}

int D0Finder::getIndex(){
	return index;
}

void D0Finder::setIndex(int i){
	if ( i > d0.size()  || i < 0){
		std::cout << "D0Finder: Error seting index " << i << " max: " << d0.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& D0Finder::getKaon(){
	return *kaon[index];
}

AA::Ptl& D0Finder::getPion(){
	return *pion[index];
}

AA::Ptl& D0Finder::getD0(){
	return *d0[index];
}

AA::Vrt& D0Finder::getVrt(){
	return *vtx[index];
}

double D0Finder::getMass(){
	return mass[index];
}

void D0Finder::fill(){
	kaon_saver.fill(getKaon());
	pion_saver.fill(getPion());
	d0_saver.fill(getD0());
	vrt_saver.fill(getVrt());
	masa = getMass();
}
