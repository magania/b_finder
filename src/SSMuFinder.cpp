/*
 * SSMuFinder.cpp
 *
 *  Created on: Nov 12, 2011
 *      Author: magania
 */

#include "SSMuFinder.h"

SSMuFinder::SSMuFinder(TTree &tree) :
	mu1_saver("mu1",tree, true, true),
	mu2_saver("mu2",tree, true, true),
	vrt_saver("mu_vrt",tree)
{
	index = -1;
	tree.Branch("mu_index", &index, "mu_index/D");
}

SSMuFinder::~SSMuFinder() {
	// TODO Auto-generated destructor stub
}

void SSMuFinder::clean(){
	v_vtx.clear();
	v_mu1.clear();
	v_mu2.clear();
	_boxp.clear();
	_boxv.clear();
}

int SSMuFinder::find(){
	clean();

	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();

	for (AA::PtlLst::const_iterator mu1 = ptl_lst->begin(); mu1 != ptl_lst->end(); ++mu1){
  		AA::MuoGlb *gMu1 =(*mu1)->muon();
                if (!gMu1) continue;
		for (AA::PtlLst::const_iterator mu2 = mu1; mu2 != ptl_lst->end(); ++mu2){
  			AA::MuoGlb *gMu2 =(*mu2)->muon();
                	if (!gMu2) continue;

			if ((*mu1) == (*mu2))
				continue;

			int qtot = ((*mu1)->q()) * ((*mu2)->q());

			if (qtot < 0)
				continue;

			//Find vertex of these 2 tracks
			AA::Vrt *vrt = _boxv.newVrt();
			if (!AA::findVrt((*mu1), (*mu2), vrt))
				continue;

			v_vtx.push_back(vrt);
			v_mu1.push_back(*mu1);
			v_mu2.push_back(*mu2);
		}
	}

	return v_vtx.size();
}

void SSMuFinder::begin(){
	index = -1;
}

bool SSMuFinder::next(){
	if (++index < v_vtx.size())
		return true;
	else
		index=-1;
	return false;
}

int SSMuFinder::getIndex(){
	return index;
}

void SSMuFinder::setIndex(int i){
	if ( i > v_vtx.size()  || i < 0){
		std::cout << "SSMuFinder: Error seting index " << i << " max: " << v_vtx.size() << std::endl;
		exit(EXIT_FAILURE);
	} else {
		index = i;
	}
}

AA::Ptl& SSMuFinder::getMu1(){
	return *v_mu1[index];
}

AA::Ptl& SSMuFinder::getMu2(){
	return *v_mu2[index];
}

AA::Vrt& SSMuFinder::getVrt(){
	return *v_vtx[index];
}

void SSMuFinder::fill(){
	mu1_saver.fill(getMu1());
	mu2_saver.fill(getMu2());
	vrt_saver.fill(getVrt());
}
