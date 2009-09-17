/*
 * PtlFinder.h
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#ifndef PTLFINDER_H_
#define PTLFINDER_H_

#include <PtlSaver.h>

#include <vector>
#include <list>

#include <Ptl.hpp>
#include <select.hpp>

#include <AA.hpp>

#include <TTree.h>

class PtlFinder {
public:
	PtlFinder(const char* name, TTree &tree, bool muon);
	virtual ~PtlFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getPtl();
	void fill();

private:
	void clean();

	char *myname;
	char *myname_d;
	PtlSaver saver;
	std::vector<AA::Ptl*> ptl_vector;
	int index;
};

#endif /* PTLFINDER_H_ */

