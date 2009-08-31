/*
 * ChiJPsiGFinder.h
 *
 *  Created on: Aug 31, 2009
 *      Author: magania
 */

#ifndef CHIJPSIGFINDER_H_
#define CHIJPSIGFINDER_H_

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
#include "JPsiFinder.h"
#include "GammaFinder.h"

class ChiJPsiGFinder {
public:
	ChiJPsiGFinder(TTree &tree, JPsiFinder&, GammaFinder&);
	virtual ~ChiJPsiGFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getJPsi();
	AA::Ptl& getGamma();
	AA::Ptl& getChi();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	GammaFinder *_gamma_finder;
	JPsiFinder *_jpsi_finder;

	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	PtlSaver chi_saver;
	VrtSaver vrt_saver;
	std::vector<Ptl*> chi;
	std::vector<int> jpsi, gamma;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	double masa;

	const static double MASS_CHI_MIN = 3.096 - 5 * 0.0693458;
	const static double MASS_CHI_MAX = 3.096 + 5 * 0.0693458 + 3;
};

#endif /* CHIJPSIGFINDER_H_ */
