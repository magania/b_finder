	/*
 * V0Finder.h
 *
 *  Created on: Feb 24, 2009
 *      Author: magania
 */

#ifndef V0FINDER_H_
#define V0FINDER_H_

#include <PtlSaver.h>

#include <vector>
#include <list>

#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class V0Finder {
public:
	V0Finder(TTree &tree);
	virtual ~V0Finder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getCPlus();
	AA::Ptl& getCMinus();
	AA::Ptl& getV0();
	void fill();

private:
	void clean();

	PtlSaver c_plus_saver, c_minus_saver, V_saver;

	std::vector<AA::Ptl*> V, c_plus, c_minus;

	int index;
};

#endif /* V0FINDER_H_ */
