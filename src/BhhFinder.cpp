/*
 * BhhFinder.cpp
 *
 *  Created on: Aug 13, 2009
 *      Author: magania
 */

#include "BhhFinder.h"

BhhFinder::BhhFinder(TTree &tree) :
	h_plus_saver("h_plus",tree, true, true),
	h_minus_saver("h_minus",tree, true, true),
	b_saver("b",tree, false, false),
	vrt_saver("b_vrt",tree),
	pv_saver("b_pv",tree)
{
	index = -1;
	tree.Branch("b_index", &index, "b_index/D");
	tree.Branch("b_mass", &masa, "b_mass/D");

	tree.Branch("mc_hplus_ptr", &mc_hplus_ptr, "mc_hplus_ptr/I");
	tree.Branch("mc_hplus_id", &mc_hplus_id, "mc_hplus_id/I");
	tree.Branch("mc_hplus_p_ptr", &mc_hplus_p_ptr, "mc_hplus_p_ptr/I");
	tree.Branch("mc_hplus_p_id", &mc_hplus_p_id, "mc_hplus_p_id/I");
	tree.Branch("mc_hplus_pp_ptr", &mc_hplus_pp_ptr, "mc_hplus_pp_ptr/I");
	tree.Branch("mc_hplus_pp_id", &mc_hplus_pp_id, "mc_hplus_pp_id/I");

	tree.Branch("mc_hminus_ptr", &mc_hminus_ptr, "mc_hminus_ptr/I");
	tree.Branch("mc_hminus_id", &mc_hminus_id, "mc_hminus_id/I");
	tree.Branch("mc_hminus_p_ptr", &mc_hminus_p_ptr, "mc_hminus_p_ptr/I");
	tree.Branch("mc_hminus_p_id", &mc_hminus_p_id, "mc_hminus_p_id/I");
	tree.Branch("mc_hminus_pp_ptr", &mc_hminus_pp_ptr, "mc_hminus_pp_ptr/I");
	tree.Branch("mc_hminus_pp_id", &mc_hminus_pp_id, "mc_hminus_pp_id/I");
}

BhhFinder::~BhhFinder() {
}

void BhhFinder::clean(){
	b.clear();
	vtx.clear();
	vtx_pv.clear();
	h_plus.clear();
	h_minus.clear();
	mass.clear();
	_boxp.clear();
	_boxv.clear();
}

int BhhFinder::find(){
	clean();

	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();

	for (AA::PtlLst::const_iterator h1 = ptl_lst->begin(); h1 != ptl_lst->end(); ++h1)
		for (AA::PtlLst::const_iterator h2 = h1; h2 != ptl_lst->end(); ++h2){
			if ((*h1) == (*h2))
				continue;

			if (!(*h1)->primaryVertex() || (*h1)->primaryVertex() != (*h2)->primaryVertex())
					continue;

			int qtot = ((*h1)->q()) * ((*h2)->q());

			if (qtot > 0)
				continue;

			if ((*h1)->pt() < H_MINIMUM_PT || (*h2)->pt() < H_MINIMUM_PT)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt vrt;
			if (!AA::findVrt((*h1), (*h2), &vrt))
				continue;

			//Construct new particle from these 2 tracks
			AA::Ptl* ptl = _boxp.newPtl();
			AA::Vrt* pve = _boxv.newVrt(vrt);
			if (!ptl->combine(pve, 0)) {
				_boxp.erase(ptl);
				_boxv.erase(pve);
				continue;
			}

			double masa = ptl->mass(AA::PI_PLUS, AA::PI_MINUS);

			if (MASS_B_MIN > masa || masa > MASS_B_MAX)
				continue;

			b.push_back(ptl);
			vtx.push_back(pve);
			vtx_pv.push_back((*h1)->primaryVertex());
			mass.push_back(masa);
			if ((*h1)->q() > 0){
				h_plus.push_back(*h1);
				h_minus.push_back(*h2);
			} else {
				h_plus.push_back(*h2);
				h_minus.push_back(*h1);
			}
		}

	return b.size();
}

void BhhFinder::begin(){
	index = -1;
}

bool BhhFinder::next(){
	if (++index < b.size())
		return true;
	else
		index=-1;
	return false;
}

int BhhFinder::getIndex(){
	return index;
}

void BhhFinder::setIndex(int i){
	if ( i > b.size()  || i < 0){
		std::cout << "BhhFinder: Error seting index " << i << " max: " << b.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& BhhFinder::getHPlus(){
	return *h_plus[index];
}

AA::Ptl& BhhFinder::getHMinus(){
	return *h_minus[index];
}

AA::Ptl& BhhFinder::getB(){
	return *b[index];
}

AA::Vrt& BhhFinder::getVrt(){
	return *vtx[index];
}

const AA::Vrt& BhhFinder::getPVrt(){
	return *vtx_pv[index];
}

double BhhFinder::getMass(){
	return mass[index];
}

void BhhFinder::fill(){
	h_plus_saver.fill(getHPlus());
	h_minus_saver.fill(getHMinus());
	b_saver.fill(getB());
	vrt_saver.fill(getVrt());
	pv_saver.fill(getPVrt());
	masa = getMass();

	double chi;
	PtlMC* h_plus_mc, *h_minus_mc;
	h_plus_mc = DecayMC::imatch(&getHPlus(),chi);
	h_minus_mc = DecayMC::imatch(&getHMinus(),chi);
	//mc_hplus_ptr = (int)h_plus_mc; mc_hminus_ptr = (int)h_minus_mc;
	mc_hplus_id = DecayMC::getIdPdg(h_plus_mc); mc_hminus_id = DecayMC::getIdPdg(h_minus_mc);

	PtlMC* h_plus_p_mc, *h_minus_p_mc;
	h_plus_p_mc = DecayMC::getParent(h_plus_mc);
	h_minus_p_mc = DecayMC::getParent(h_minus_mc);
	//mc_hplus_p_ptr = (int) h_plus_p_mc; mc_hminus_p_ptr = (int)h_minus_p_mc;
	mc_hplus_p_id = DecayMC::getIdPdg(h_plus_p_mc); mc_hminus_p_id = DecayMC::getIdPdg(h_minus_p_mc);

	PtlMC* h_plus_pp_mc, *h_minus_pp_mc;
	h_plus_pp_mc = DecayMC::getParent(h_plus_p_mc);
	h_minus_pp_mc = DecayMC::getParent(h_minus_p_mc);
	//mc_hplus_pp_ptr = (int) h_plus_pp_mc; mc_hminus_pp_ptr = (int)h_minus_p_mc;
	mc_hplus_pp_id = DecayMC::getIdPdg(h_plus_pp_mc); mc_hminus_pp_id = DecayMC::getIdPdg(h_minus_p_mc);

}
