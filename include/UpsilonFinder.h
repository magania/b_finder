/*
 * UpsilonFinder.h
 *
 *  Created on: Aug 26, 2009
 *      Author: magania
 */

#ifndef UPSILONFINDER_H_
#define UPSILONFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <AA.hpp>
#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>
#include <Muo.hpp>

#include <TTree.h>

class UpsilonFinder {
public:
	UpsilonFinder(TTree &tree);
	virtual ~UpsilonFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getMuPlus();
	AA::Ptl& getMuMinus();
	AA::Ptl& getUpsilon();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	double masa;

	PtlSaver mu_plus_saver, mu_minus_saver, upsilon_saver;
	VrtSaver vrt_saver;
	std::vector<AA::Ptl*> upsilon, mu_plus, mu_minus;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	const static double MASS_UPSILON_MIN = 9.4603 - 6 * 0.5;
	const static double MASS_UPSILON_MAX = 9.4603 + 6 * 0.5;

	const static double MU_MINIMUM_NSEG = 3;
	const static double MU_MINIMUM_PT = 1;
};

#endif /* UPSILONFINDER_H_ */
