/*
 * BhhFinder.h
 *
 *  Created on: Aug 13, 2009
 *      Author: magania
 */

#ifndef BHHFINDER_H_
#define BHHFINDER_H_

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

class BhhFinder {
public:
	BhhFinder(TTree &tree);
	virtual ~BhhFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getHPlus();
	AA::Ptl& getHMinus();
	AA::Ptl& getB();
	AA::Vrt& getVrt();
	const AA::Vrt& getPVrt();
	double getMass();
	void fill();

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	PtlSaver h_plus_saver, h_minus_saver, b_saver;
	VrtSaver vrt_saver, pv_saver;
	std::vector<AA::Ptl*> b, h_plus, h_minus;
	std::vector<AA::Vrt*> vtx;
	std::vector<const AA::Vrt*> vtx_pv;
	std::vector<double> mass;

	int index;


	double masa;
	int mc_hplus_ptr, mc_hplus_id, mc_hminus_ptr, mc_hminus_id;
	int mc_hplus_p_ptr, mc_hplus_p_id, mc_hminus_p_ptr, mc_hminus_p_id;
	int mc_hplus_pp_ptr, mc_hplus_pp_id, mc_hminus_pp_ptr, mc_hminus_pp_id;

	const static double MASS_B_MIN = 4.5;
	const static double MASS_B_MAX = 6.2;

	const static double H_MINIMUM_PT = 0;
};

#endif /* BHHFINDER_H_ */
