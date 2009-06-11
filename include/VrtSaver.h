/*
 * VrtSaver.h
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#ifndef VRTSAVER_H_
#define VRTSAVER_H_

#include <vector>

#include <Vrt.hpp>

#include <TTree.h>

class VrtSaver {
public:
	VrtSaver(const char *name, TTree &tree);
	virtual ~VrtSaver();

	void fill(const AA::Vrt &v);

private:
	const char* strccat(const char *a, const char *b);

	std::vector<char*> names;

	double x, y, z;
	double v11, v12, v13, v22, v23, v33;
	double chi2;

};

#endif /* VRTSAVER_H_ */
