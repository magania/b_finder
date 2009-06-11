/*
 * V0Finder.cpp
 *
 *  Created on: Feb 24, 2009
 *      Author: magania
 */

#include "V0Finder.h"

V0Finder::V0Finder(TTree &tree) :
	c_plus_saver("c_plus",tree, true, true),
	c_minus_saver("c_minus",tree, true, true),
	V_saver("v0",tree, false, false)
{
	AA::findJPsiFlag = true;
	AA::selectJPsiFlag   = true;
	index = -1;
	tree.Branch("V_index", &index, "V_index/D");
}

V0Finder::~V0Finder() {
	// TODO Auto-generated destructor stub
}

void V0Finder::clean(){
	V.clear();
	c_plus.clear();
	c_minus.clear();
}

int V0Finder::find(){
	clean();

	// Get list of all V0 candidates
	const AA::PtlLst* v0finder_list = AA::v0Finder.particles();

	for (AA::PtlLstCIt v_l = v0finder_list->begin();
			v_l!= v0finder_list->end(); ++v_l) {

		//Pointers to 2 particles
		const std::list<AA::Chain*>* v_children = (*v_l)->children();
		AA::Ptl* c1 = AA::particle(v_children->front()->track());
		AA::Ptl* c2 = AA::particle(v_children->back()->track());

		if (c1->q() > 0 && c2->q() < 0){
			c_plus.push_back(c1);
			c_minus.push_back(c2);
			V.push_back(*v_l);
		}

		if (c2->q() > 0 && c1->q() < 0){
			c_plus.push_back(c2);
			c_minus.push_back(c1);
			V.push_back(*v_l);
		}
	}
	return V.size();
}

void V0Finder::begin(){
	index = -1;
}

bool V0Finder::next(){
	if (++index < V.size())
		return true;
	else
		index=-1;
	return false;
}

AA::Ptl& V0Finder::getCPlus(){
	return *c_plus[index];
}

AA::Ptl& V0Finder::getCMinus(){
	return *c_minus[index];
}

AA::Ptl& V0Finder::getV0(){
	return *V[index];
}

void V0Finder::fill(){
	c_plus_saver.fill(getCPlus());
	c_minus_saver.fill(getCMinus());
	V_saver.fill(getV0());
}
