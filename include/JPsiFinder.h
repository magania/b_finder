/*
 * JPsiFinder.h
 *
 *  Created on: Feb 21, 2009
 *      Author: magania
 */

#ifndef JPSIFINDER_H_
#define JPSIFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class JPsiFinder {
public:
	JPsiFinder(TTree &tree);
	virtual ~JPsiFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getMuPlus();
	AA::Ptl& getMuMinus();
	AA::Ptl& getJPsi();
	double getMass();
	void fill();

private:
	void clean();

    double tree_jpsi_mass;
	PtlSaver mu_plus_saver, mu_minus_saver, jpsi_saver;
    VrtSaver vrt_saver;

    double mu_leading_pt, mu_trailing_pt;
    double mu_leading_ptot, mu_trailing_ptot;

	std::vector<AA::Ptl*> jpsi, mu_plus, mu_minus;
	std::vector<double> jpsi_mass;

	int index;

	static const double JPSI_WINDOW_MIN = 3.096 - 5 * 0.0693458; // Sigma from MC
	static const double JPSI_WINDOW_MAX = 3.096 + 5 * 0.0693458;
};

#endif /* JPSIFINDER_H_ */
