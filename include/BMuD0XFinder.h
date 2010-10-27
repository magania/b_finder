/*
 * BMuD0XFinder.h
 *
 *  Created on: Oct 22, 2010
 *      Author: magania
 */

#ifndef BMUD0XFINDER_H_
#define BMUD0XFINDER_H_

#include <vector>
#include <list>

#include <EvtSaver.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <PtlFinder.h>
#include <D0Finder.h>
#include <TagSaver.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>

#include <AA.hpp>
#include <Trg.hpp>
#include <Muo.hpp>
#include <Ptl.hpp>
#include <DST.hpp>
#include <Field.hpp>
#include <Spot.hpp>
#include <select.hpp>
#include <Monitor.hpp>
#include <V0Finder.hpp>

#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include "TLorentzRotation.h"

class BMuD0XFinder {
public:
	BMuD0XFinder(TTree&, D0Finder&, PtlFinder&, PtlFinder&);
	virtual ~BMuD0XFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getB();
	AA::Ptl& getD0();
	AA::Ptl& getPion();
	AA::Ptl& getMu();
	double getDstarMass();
	void fill();

private:
	void clean();

        static const int MUON_NSEG_MIN = 3;
        static const double PT_MU_MIN = 2.0;
        static const double PTOT_MU_MIN = 3.0;
        static const double MUON_CHISQLOC_MIN = -0.5;
        static const int MUON_NSMT_MIN = 2;
        static const int MUON_NCFT_MIN = 2;

        static const double PI_K_CHI2VRT_MIN = 2.0;

        static const double COSXY_D0_MIN = 0.9;
        static const double COSXY_B_MIN = 0.95;
        static const double SIGMA_B_LXY = 9.0;
        static const double SIGMA_D0_LXY = 9.0;

	static const double MASS_B_MIN = 2.0;
	static const double MASS_B_MAX = 5.5;

        static const double CHI2_BVRT_MIN = 9.0;

	static const double MASS_DSmD0_MIN = 0.13;
	static const double MASS_DSmD0_MAX = 0.17;


        AA::PtlBox _boxp;
        AA::VrtBox _boxv;
	int index;

	D0Finder  *_d0_finder;
	PtlFinder *_pi_finder;
	PtlFinder *_mu_finder;

	PtlSaver b_saver;
	VrtSaver vb_saver;
	VrtSaver vp_saver;
	TagSaver tag_saver;

	std::vector<int> v_d0_index, v_mu_index;
	std::vector<AA::Ptl*> v_b;
	std::vector<double> v_mdstar;

	/* -- Info missing by the savers -- */
	double dstar_mass, b_dl, b_vdl, d0_dl, d0_vdl;
};

#endif /* BMUD0XFINDER_H_ */
