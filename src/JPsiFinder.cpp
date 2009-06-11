/*
 * JPsiFinder.cpp
 *
 *  Created on: Feb 21, 2009
 *      Author: magania
 */

#include "JPsiFinder.h"

JPsiFinder::JPsiFinder(TTree &tree) :
	mu_plus_saver("mu_plus",tree, true, true),
	mu_minus_saver("mu_minus",tree, true, true),
	jpsi_saver("jpsi", tree, false, false),
	vrt_saver("jpsi_vrt", tree)
{
	AA::findJPsiFlag = true;
	AA::selectJPsiFlag   = true;
	index = -1;
	tree.Branch("jpsi_mass", &tree_jpsi_mass, "jpsi_mass/D");
	tree.Branch("jpsi_index", &index, "jpsi_index/D");

	tree.Branch("mu_leading_pt", &mu_leading_pt, "mu_leading_pt/D");
	tree.Branch("mu_trailing_pt", &mu_trailing_pt, "mu_trailing_pt/D");
	tree.Branch("mu_leading_ptot", &mu_leading_ptot, "mu_leading_ptot/D");
	tree.Branch("mu_trailing_ptot", &mu_trailing_ptot, "mu_trailing_ptot/D");
}

JPsiFinder::~JPsiFinder() {
	// TODO Auto-generated destructor stub
}

void JPsiFinder::clean(){
	jpsi_mass.clear();
	jpsi.clear();
	mu_plus.clear();
	mu_minus.clear();
}

int JPsiFinder::find(){
	clean();

	// Get list of all J/psi candidates
	const AA::PtlLst* jpsifinder_list = AA::jPsiFinder.particles();

	for (AA::PtlLstCIt jpsi_l = jpsifinder_list->begin();
			jpsi_l!= jpsifinder_list->end(); ++jpsi_l) {

		/* Mass constraint on J/Psi */
		double mjpsi = (*jpsi_l)->mass(AA::MU_PLUS, AA::MU_MINUS);

		if (mjpsi < JPSI_WINDOW_MIN || mjpsi > JPSI_WINDOW_MAX)
			continue;

		//Pointers to 2 muons from J/psi -> mu+ mu-
		const std::list<AA::Chain*>* jpsi_children = (*jpsi_l)->children();
		AA::Ptl* mu1 = AA::particle(jpsi_children->front()->track());
		AA::Ptl* mu2 = AA::particle(jpsi_children->back()->track());

		if (mu1->q() > 0 && mu2->q() < 0){
			mu_plus.push_back(mu1);
			mu_minus.push_back(mu2);
			jpsi.push_back(*jpsi_l);
			jpsi_mass.push_back(mjpsi);
		}

		if (mu2->q() > 0 && mu1->q() < 0){
			mu_plus.push_back(mu2);
			mu_minus.push_back(mu1);
			jpsi.push_back(*jpsi_l);
			jpsi_mass.push_back(mjpsi);
		}
	}
	return jpsi.size();
}

void JPsiFinder::begin(){
	index = -1;
}

bool JPsiFinder::next(){
	if (++index < jpsi.size())
		return true;
	else
		index=-1;
	return false;
}

AA::Ptl& JPsiFinder::getMuPlus(){
	return *mu_plus[index];
}

AA::Ptl& JPsiFinder::getMuMinus(){
	return *mu_minus[index];
}

AA::Ptl& JPsiFinder::getJPsi(){
	return *jpsi[index];
}

double JPsiFinder::getMass(){
	return jpsi_mass[index];
}

void JPsiFinder::fill(){
	mu_plus_saver.fill(getMuPlus());
	mu_minus_saver.fill(getMuMinus());
	jpsi_saver.fill(getJPsi());
	tree_jpsi_mass = getMass();
	vrt_saver.fill(*getJPsi().decayVertex());

	if (getMuPlus().pt() > getMuMinus().pt()){
		mu_leading_pt = getMuPlus().pt();
		mu_trailing_pt = getMuMinus().pt();
	} else {
		mu_leading_pt = getMuMinus().pt();
		mu_trailing_pt = getMuPlus().pt();
	}
	if (getMuPlus().ptot() > getMuMinus().ptot()){
		mu_leading_ptot = getMuPlus().ptot();
		mu_trailing_ptot = getMuMinus().ptot();
	} else {
		mu_leading_ptot = getMuMinus().ptot();
		mu_trailing_ptot = getMuPlus().ptot();
	}
}
