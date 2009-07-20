/*
 * KstarFinder.h
 *
 *  Created on: Jul 16, 2009
 *      Author: magania
 */

#ifndef KSTARFINDER_H_
#define KSTARFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <AA.hpp>
#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class KstarFinder {
public:
	KstarFinder(TTree &tree);
	virtual ~KstarFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getKaon();
	AA::Ptl& getPion();
	AA::Ptl& getKstar();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	double masa;

	PtlSaver kaon_saver, pion_saver, kstar_saver;
	VrtSaver vrt_saver;
	std::vector<AA::Ptl*> kstar, kaon, pion;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	const static double MASS_KSTAR_MIN = 0.89166 - 10 * 0.0356; // Obtained from MC
	const static double MASS_KSTAR_MAX = 0.89166 + 10 * 0.0356;

	const static double KAON_MINIMUM_PT = 0.4;
	const static double PION_MINIMUM_PT = 0.4;
	double MINIMUM_PT;
};

#endif /* KSTARFINDER_H_ */
