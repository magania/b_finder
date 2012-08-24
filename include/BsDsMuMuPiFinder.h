/*
 * BsDsMuMuPiFinder.h
 *
 *  Created on: Nov 12, 2011
 *      Author: magania
 */

#ifndef BSDSMUMUFINDER_H_
#define BSDSMUMUFINDER_H_

#include <vector>
#include <list>

#include <EvtSaver.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <PtlFinder.h>
#include <SSMuFinder.h>
#include <DsFinder.h>
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

class BsDsMuMuPiFinder {
public:
	BsDsMuMuPiFinder(TTree&, DsFinder&, SSMuFinder&, PtlFinder&);
	virtual ~BsDsMuMuPiFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getB();
	AA::Ptl& getDs();
	AA::Ptl& getMu1();
	AA::Ptl& getMu2();
	AA::Ptl& getPion();
	double getMass();
	void fill();

private:
	void clean();

	static const double PDG_BD_MASS = 5.3663;
	static const double PDG_KAON_MASS = 0.493677;
	static const double PDG_PION_MASS = 0.13957018;
	static const double PDG_MU_MASS = 0.1056583668;

	static const double MASS_B_MIN = 5.0;
	static const double MASS_B_LEFT = 5.16;
	static const double MASS_B_RIGHT = 5.56;
	static const double MASS_B_MAX = 5.8;

	static const double PT_PION_MIN = 0.5;

	static const bool blind = false;

    AA::PtlBox _boxp;
    AA::VrtBox _boxv;
	int index;

	DsFinder *_ds_finder;
	SSMuFinder *_ssmu_finder;
	PtlFinder *_pion_finder;

	PtlSaver b_saver;
	VrtSaver vb_saver;
	VrtSaver vp_saver;

	std::vector<int> v_ds_index, v_ssmu_index, v_pion_index;
	std::vector<AA::Ptl*> v_b;
	std::vector<double> v_mb, v_emb;

	/* -- Info missing by the savers -- */
	double b_mass, b_mass_error, b_dl;
};

#endif /* BSDSMUMUPIFINDER_H_ */
