/*
 * XYGammaFinder.h
 *
 *  Created on: Aug 26, 2009
 *      Author: magania
 */

#ifndef XYGAMMAFINDER_H_
#define XYGAMMAFINDER_H_

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
#include "UpsilonFinder.h"
#include "GammaFinder.h"

class XYGammaFinder {
public:
	XYGammaFinder(TTree &tree, UpsilonFinder&, GammaFinder&);
	virtual ~XYGammaFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getUpsilon();
	AA::Ptl& getGamma();
	AA::Ptl& getX();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	GammaFinder *_gamma_finder;
	UpsilonFinder *_upsilon_finder;

	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	PtlSaver x_saver;
	VrtSaver vrt_saver, pv_saver;
	std::vector<Ptl*> X;
	std::vector<int> upsilon, gamma;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	double masa;

	const static double MASS_X_MIN = 9.4603 - 6 * 0.5;
	const static double MASS_X_MAX = 9.4603 + 6 * 0.5 + 3;
};

#endif /* XYGAMMAFINDER_H_ */
