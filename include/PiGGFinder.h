/*
 * PiGGFinder.h
 *
 *  Created on: Aug 28, 2009
 *      Author: magania
 */

#ifndef PIGGFINDER_H_
#define PIGGFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <AA.hpp>
#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

#include "DecayMC.h"
#include "GammaFinder.h"

class PiGGFinder {
public:
	PiGGFinder(TTree &tree, GammaFinder&);
	virtual ~PiGGFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getGamma1();
	AA::Ptl& getGamma2();
	AA::Ptl& getPi();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	GammaFinder *_gamma_finder;

	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	PtlSaver pi_saver;
	VrtSaver vrt_saver;
	std::vector<Ptl*> Pi;
	std::vector<int> gamma1, gamma2;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	double masa;

	const static double MASS_PI_MIN = 0;//0.1349766 - 3*;
	const static double MASS_PI_MAX = 0.8;//0.1349766 + 3*;
};

#endif /* PIGGFINDER_H_ */
