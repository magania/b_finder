/*
 * GammaFinder.h
 *
 *  Created on: Aug 25, 2009
 *      Author: magania
 */

#ifndef GAMMAFINDER_H_
#define GAMMAFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class GammaFinder {
public:
	GammaFinder(TTree &tree,bool two = false);
	virtual ~GammaFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getEPlus();
	AA::Ptl& getEMinus();
	AA::Ptl& getGamma();
	AA::Vrt& getVrt();
	double getMass();
	void fill();
	void fill(int, int);

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	bool _two;
	double tree_gamma_mass;
	double tree_gamma2_mass;
	int index2;
	PtlSaver e_plus_saver, e_minus_saver, gamma_saver;
    VrtSaver vrt_saver;
	PtlSaver *e2_plus_saver, *e2_minus_saver, *gamma2_saver;
    VrtSaver *vrt2_saver;

	std::vector<AA::Ptl*> gamma, e_plus, e_minus;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> gamma_mass;

	int index;

	static const double MASS_GAMMA_MAX = 0.100;

	static const double R_ORIGIN_MIN = 1.5;
	static const double R_ORIGIN_MAX = 15.0;
};

#endif /* GAMMAFINDER_H_ */
