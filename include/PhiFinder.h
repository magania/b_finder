/*
 * PhiFinder.h
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#ifndef PHIFINDER_H_
#define PHIFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <AA.hpp>
#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class PhiFinder {
public:
	PhiFinder(TTree &tree);
	virtual ~PhiFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getKPlus();
	AA::Ptl& getKMinus();
	AA::Ptl& getPhi();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	double masa;

	PtlSaver k_plus_saver, k_minus_saver, phi_saver;
	VrtSaver vrt_saver;
	std::vector<AA::Ptl*> phi, k_plus, k_minus;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	const static double MASS_PHI_MIN = 1.019455 - 6 * 0.005;
	const static double MASS_PHI_MAX = 1.019455 + 6 * 0.005;

	const static double K_MINIMUM_PT = 0.4;
};

#endif /* PHIFINDER_H_ */
