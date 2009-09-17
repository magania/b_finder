/*
 * BChiKFinder.h
 *
 *  Created on: Sep 17, 2009
 *      Author: magania
 */

#ifndef BCHIKFINDER_H_
#define BCHIKFINDER_H_

#include <vector>
#include <list>

#include <EvtSaver.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <PtlFinder.h>
#include <ChiJPsiGFinder.h>
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

class BChiKFinder {
public:
	BChiKFinder(TTree&, ChiJPsiGFinder&, PtlFinder&);
	virtual ~BChiKFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getB();
	AA::Ptl& getChi();
	AA::Ptl& getKaon();
	double getMass();
	void fill();

private:
	void clean();

	static const double PDG_BD_MASS = 5.27950;
	static const double PDG_KAON_MASS = 0.493677;
	static const double PDG_PION_MASS = 0.13957018;
	static const double PDG_MU_MASS = 0.1056583668;

	static const double MASS_B_MIN = 5.0;
	static const double MASS_B_MAX = 5.8;

	static const double PT_KAON_MIN = 0.4;

    AA::PtlBox _boxp;
    AA::VrtBox _boxv;
	int index;

	ChiJPsiGFinder *_chi_finder;
	PtlFinder *_kaon_finder;

	PtlSaver b_saver;
	VrtSaver vb_saver;
	VrtSaver vp_saver;

	vector<int> v_chi_index, v_kaon_index;
	vector<Ptl*> v_b;
	vector<double> v_mb, v_emb;

	/* -- Info missing by the savers -- */
	double b_mass, b_mass_error, b_dl;
};

#endif /* BChiKFINDER_H_ */
