/*
 * BdJPsiKstarFinder.h
 *
 *  Created on: Jul 16, 2009
 *      Author: magania
 */

#ifndef BDJPSIKSTARFINDER_H_
#define BDJPSIKSTARFINDER_H_

#include <vector>
#include <list>

#include <BdJPsiKstarMCFinder.h>
#include <EvtSaver.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <KstarFinder.h>
#include <JPsiFinder.h>
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
#if defined(P17) || defined(P20) || defined(MC)
#include "IPcal.hpp"
#endif

#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include "TLorentzRotation.h"

class BdJPsiKstarFinder {
public:
	BdJPsiKstarFinder(JPsiFinder *jpsi, KstarFinder *kstar, TTree &tree, BdJPsiKstarMCFinder *mc = 0);
	virtual ~BdJPsiKstarFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getBs();
	AA::Ptl& getJPsi();
	AA::Ptl& getMuPlus();
	AA::Ptl& getMuMinus();
	AA::Ptl& getKstar();
	AA::Ptl& getKaon();
	AA::Ptl& getPion();
	double getMass();
	void fill();

private:
	void clean();
	/* -- Legendary threeAngles function:
	 * mu_p: positive muon
	 * mu_n: negative muon
	 * K_p : positive Kaon (Kaon for Bd)
	 * K_n : negative Kaon (Pion for Bd)
	 *
	 * w1 = phi, w2 = cos(theta), w3 = cos(psi)
	 *
	 */
	void threeAngles(TLorentzVector mu_p,
			           TLorentzVector mu_n,
			           TLorentzVector K_p,
			           TLorentzVector K_n,
			           Double_t &w1,Double_t &w2,Double_t &w3);

	static const double PDG_BD_MASS = 5.27950;
	static const double PDG_KAON_MASS = 0.493677;
	static const double PDG_PION_MASS = 0.13957018;
	static const double PDG_MU_MASS = 0.1056583668;

	static const double MASS_BD_MIN = 5.0;
	static const double MASS_BD_MAX = 5.8;

    AA::PtlBox _boxp;
    AA::VrtBox _boxv;
	int index;

	JPsiFinder *jpsi_finder;
	KstarFinder *kstar_finder;
	BdJPsiKstarMCFinder *mc_finder;
	TagSaver tag_saver; // Yorch tag saver.
	PtlSaver bd_saver;
	VrtSaver vb_saver;
	VrtSaver vp_saver;

    vector<HepVector> v_bd_mom;
	vector<int> v_jpsi_index, v_kstar_index;
	vector<Ptl*> v_bd, v_jpsi, v_muplus, v_muminus, v_kaon, v_pion;
	vector<Vrt*> v_bd_vrt, v_jpsi_kn_vrt, v_jpsi_pi_vrt, v_bd_pv;
	vector<double> v_mb, v_emb;
	vector<double> v_muplus_cpx, v_muplus_cpy, v_muplus_cpz;
	vector<double> v_muminus_cpx, v_muminus_cpy, v_muminus_cpz;
	vector<double> v_kaon_cpx, v_kaon_cpy, v_kaon_cpz;
	vector<double> v_pion_cpx, v_pion_cpy, v_pion_cpz;

	/* -- Info missing by the savers -- */
	double bd_mass, bd_mass_error, bd_pdl, bd_epdl, bd_ucpt, bd_ucptot;
    double mu_plus_dR, mu_minus_dR, kaon_dR, pion_dR;
	double mu_plus_cpx, mu_plus_cpy, mu_plus_cpz;
	double mu_minus_cpx, mu_minus_cpy, mu_minus_cpz;
	double kaon_cpx, kaon_cpy, kaon_cpz;
	double pion_cpx, pion_cpy, pion_cpz;
	double bd_iso, bd_iso_drmax, bd_iso_75;
	double bd_iso_pv, bd_iso_drmax_pv, bd_iso_75_pv;
	double bd_jpsikn_chi2, bd_jpsipi_chi2;
	double bd_angle_phi, bd_angle_ctheta, bd_angle_cpsi;
	double kstar_mass_corrected, kstar_mass_corrected_error;
	double mu_plus_cpt, mu_minus_cpt, kaon_cpt, pion_cpt, jpsi_cpt, kstar_cpt;

//    double tree_jpsi_mass;

//    double mu_leading_pt, mu_trailing_pt;
//    double mu_leading_ptot, mu_trailing_ptot;
};

#endif /* BDJPSIKstarFINDER_H_ */
