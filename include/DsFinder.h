/*
 * DsFinder.h
 *
 *  Created on: Nov 12, 2011
 *      Author: magania
 */

#ifndef DSFINDER_H_
#define DSFINDER_H_

#include <vector>
#include <list>

#include <EvtSaver.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <PtlFinder.h>
#include <PhiFinder.h>
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

class DsFinder {
public:
	DsFinder(TTree&, PhiFinder&, PtlFinder&);
	virtual ~DsFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getDs();
	AA::Ptl& getPhi();
	AA::Ptl& getKPlus();
	AA::Ptl& getKMinus();
	AA::Ptl& getPion();
	double getMass();
	void fill();

private:
	void clean();

	static const double PDG_DS_MASS = 1.96845;
	static const double PDG_KAON_MASS = 0.493677;
	static const double PDG_PION_MASS = 0.13957018;

	static const double MASS_DS_MIN = 1.96845 - 1;
	static const double MASS_DS_MAX = 1.96845 + 1;

	static const double PT_PION_MIN = 0.4;

    AA::PtlBox _boxp;
    AA::VrtBox _boxv;
	int index;

	PhiFinder *_phi_finder;
	PtlFinder *_pion_finder;

	PtlSaver ds_saver;

	std::vector<int> v_ds_index, v_phi_index, v_pion_index;
	std::vector<AA::Ptl*> v_ds;
	std::vector<double> v_mds, v_emds;

	/* -- Info missing by the savers -- */
	double ds_mass, ds_mass_error;
};

#endif /* DSFINDER_H_ */
