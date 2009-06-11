/*
 * PtlFinder.cpp
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#include "PtlFinder.h"

PtlFinder::PtlFinder(const char* name, TTree &tree, bool muon) :
	saver(name,	tree, true, muon)
{
		index = -1;

		myname = new char[20];
		myname_d = new char[20];
		strcpy(myname, name);
		const char *name_index = strcat(myname,"_index");
		strcpy(myname_d,myname);
		const char *name_index_d = strcat(myname_d,"/D");

		tree.Branch(name_index, &index, name_index_d);
}

PtlFinder::~PtlFinder() {
	delete myname;
	delete myname_d;
}

void PtlFinder::clean(){
	ptl_vector.clear();
}

int PtlFinder::find(){
	clean();

	const AA::PtlLst ptl_lst = *AA::ptlBox.particles();

	for (AA::PtlLstCIt p = ptl_lst.begin(); p != ptl_lst.end(); ++p)
		ptl_vector.push_back(*p);
	return ptl_vector.size();
}

void PtlFinder::begin(){
	index = -1;
}

bool PtlFinder::next(){
	if (++index < ptl_vector.size())
		return true;
	else
		index=-1;
	return false;
}

AA::Ptl& PtlFinder::getPtl(){
	return *ptl_vector[index];
}

void PtlFinder::fill(){
	saver.fill(getPtl());
}
