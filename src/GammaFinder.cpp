/*
 * GammaFinder.cpp
 *
 *  Created on: Aug 25, 2009
 *      Author: magania
 */

#include "GammaFinder.h"

GammaFinder::GammaFinder(TTree &tree, bool two) :
	e_plus_saver("e_plus",tree, true, false),
	e_minus_saver("e_minus",tree, true, false),
	gamma_saver("gamma", tree, false, false),
	vrt_saver("gamma_vrt", tree)
{
	_two = two;
	index = -1;
	tree.Branch("gamma_mass", &tree_gamma_mass, "gamma_mass/D");
	tree.Branch("gamma_index", &index, "gamma_index/D");

	if (_two){
		tree.Branch("gamma2_mass", &tree_gamma2_mass, "gamma2_mass/D");
		tree.Branch("gamma2_index", &index2, "gamma2_index/D");
		e2_plus_saver = new PtlSaver("e2_plus",tree, true, false);
		e2_minus_saver = new PtlSaver("e2_minus",tree, true, false);
		gamma2_saver = new PtlSaver("gamma2", tree, false, false);
		vrt2_saver = new VrtSaver("gamma2_vrt", tree);
	}
}

GammaFinder::~GammaFinder() {
	// TODO Auto-generated destructor stub
}

void GammaFinder::clean(){
	gamma_mass.clear();
	gamma.clear();
	e_plus.clear();
	e_minus.clear();
	vtx.clear();
	_boxp.clear();
	_boxv.clear();
}

int GammaFinder::find(){
	clean();

	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();

	for (AA::PtlLst::const_iterator e1 = ptl_lst->begin(); e1 != ptl_lst->end(); ++e1)
		for (AA::PtlLst::const_iterator e2 = e1; e2 != ptl_lst->end(); ++e2){
			if ((*e1) == (*e2))
				continue;

			int qtot = ((*e1)->q()) * ((*e2)->q());
			if (qtot > 0)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt vrt;
			if (!AA::findVrt((*e1), (*e2), &vrt))
				continue;

			//Construct new particle from these 2 tracks
			AA::Ptl* ptl = _boxp.newPtl();
			AA::Vrt* pve = _boxv.newVrt(vrt);
			if (!ptl->combine(pve, 0)) {
				_boxp.erase(ptl);
				_boxv.erase(pve);
				continue;
			}

			Double_t gamma_vtx_r = sqrt(vrt.x(1)*vrt.x(1) + vrt.x(2)*vrt.x(2));
			if (gamma_vtx_r < R_ORIGIN_MIN || gamma_vtx_r > R_ORIGIN_MAX)
				continue;

		    double masa = ptl->mass(AA::E_PLUS, AA::E_MINUS);

			if (masa > MASS_GAMMA_MAX)
				continue;

			gamma.push_back(ptl);
			vtx.push_back(pve);
			gamma_mass.push_back(masa);
			if ((*e1)->q() > 0){
				e_plus.push_back(*e1);
				e_minus.push_back(*e2);
			} else {
				e_plus.push_back(*e2);
				e_minus.push_back(*e1);
			}
	}
	return gamma.size();
}

void GammaFinder::begin(){
	index = -1;
}

bool GammaFinder::next(){
	if (++index < gamma.size())
		return true;
	else
		index=-1;
	return false;
}

int GammaFinder::getIndex(){
	return index;
}

void GammaFinder::setIndex(int i){
	if ( i > gamma.size()  || i < 0){
		std::cout << "GammaFinder: Error seting index " << i << " max: " << gamma.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& GammaFinder::getEPlus(){
	return *e_plus[index];
}

AA::Ptl& GammaFinder::getEMinus(){
	return *e_minus[index];
}

AA::Ptl& GammaFinder::getGamma(){
	return *gamma[index];
}

AA::Vrt& GammaFinder::getVrt(){
	return *vtx[index];
}

double GammaFinder::getMass(){
	return gamma_mass[index];
}

void GammaFinder::fill(){
	e_plus_saver.fill(getEPlus());
	e_minus_saver.fill(getEMinus());
	gamma_saver.fill(getGamma());
	tree_gamma_mass = getMass();
	vrt_saver.fill(getVrt());
}

void GammaFinder::fill(int i, int j){
	if (_two){
		setIndex(i);
		e_plus_saver.fill(getEPlus());
		e_minus_saver.fill(getEMinus());
		gamma_saver.fill(getGamma());
		tree_gamma_mass = getMass();
		vrt_saver.fill(getVrt());
		setIndex(j);
		e2_plus_saver->fill(getEPlus());
		e2_minus_saver->fill(getEMinus());
		gamma2_saver->fill(getGamma());
		tree_gamma2_mass = getMass();
		vrt2_saver->fill(getVrt());
		index2=j;
	} else {
		fill();
	}
}
