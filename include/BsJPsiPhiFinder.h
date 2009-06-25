/*
 * BsJPsiPhiFinder.h
 *
 *  Created on: Jun 22, 2009
 *      Author: magania
 */

#ifndef BSJPSIPHIFINDER_H_
#define BSJPSIPHIFINDER_H_

#include <vector>
#include <list>

#include <BsJPsiPhiMCFinder.h>
#include <EvtSaver.h>
#include <PtlSaver.h>
#include <VrtSaver.h>
#include <PhiFinder.h>
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

class BsJPsiPhiFinder {
public:
	BsJPsiPhiFinder(JPsiFinder *jpsi, PhiFinder *phi, TTree &tree, BsJPsiPhiMCFinder *mc = 0);
	virtual ~BsJPsiPhiFinder();

	int find();
	void begin();
	bool next();

	AA::Ptl& getBs();
	AA::Ptl& getJPsi();
	AA::Ptl& getMuPlus();
	AA::Ptl& getMuMinus();
	AA::Ptl& getPhi();
	AA::Ptl& getKPlus();
	AA::Ptl& getKMinus();
	double getMass();
	void fill();

private:
	void clean();
	/* -- Legendary threeAngles function:
	 * mu_p: positive muon
	 * mu_n: negative muon
	 * K_p : positive Kaon
	 * K_n : negative Kaon
	 *
	 * w1 = phi, w2 = cos(theta), w3 = cos(psi)
	 *
	 */
	void threeAngles(TLorentzVector mu_p,
			           TLorentzVector mu_n,
			           TLorentzVector K_p,
			           TLorentzVector K_n,
			           Double_t &w1,Double_t &w2,Double_t &w3);

	static const double PDG_BS_MASS = 5.36636;
	static const double PDG_KAON_MASS = 0.493677;
	static const double PDG_MU_MASS = 0.1056583668;

	static const double MASS_BS_MIN = 5.0;
	static const double MASS_BS_MAX = 5.8;

	int index;

	JPsiFinder *jpsi_finder;
	PhiFinder *phi_finder;
	BsJPsiPhiMCFinder *mc_finder;
	TagSaver tag_saver; // Yorch tag saver.
	PtlSaver bs_saver;
	VrtSaver vb_saver;
	VrtSaver vp_saver;

	vector<int> v_jpsi_index, v_phi_index;
	vector<Ptl*> v_bs, v_jpsi, v_muplus, v_muminus, v_kplus, v_kminus;
	vector<Vrt*> v_bs_vrt, v_jpsi_kp_vrt, v_jpsi_km_vrt, v_bs_pv;
	vector<double> v_lhtag, v_mb, v_emb;
	vector<double> v_muplus_cpx, v_muplus_cpy, v_muplus_cpz;
	vector<double> v_muminus_cpx, v_muminus_cpy, v_muminus_cpz;
	vector<double> v_kplus_cpx, v_kplus_cpy, v_kplus_cpz;
	vector<double> v_kminus_cpx, v_kminus_cpy, v_kminus_cpz;

	/* -- Info missing by the savers -- */
	double bs_mass, bs_mass_error, bs_lhtag, bs_pdl, bs_epdl;
	double mu_plus_cpx, mu_plus_cpy, mu_plus_cpz;
	double mu_minus_cpx, mu_minus_cpy, mu_minus_cpz;
	double k_plus_cpx, k_plus_cpy, k_plus_cpz;
	double k_minus_cpx, k_minus_cpy, k_minus_cpz;
	double bs_iso, bs_iso_drmax, bs_iso_75;
    double bs_iso_pv, bs_iso_drmax_pv, bs_iso_75_pv;
	double bs_jpsikp_chi2, bs_jpsikm_chi2;
	double bs_angle_phi, bs_angle_ctheta, bs_angle_cpsi;
	double phi_mass_corrected, phi_mass_corrected_error;

//    double tree_jpsi_mass;

    double mu_leading_pt, mu_trailing_pt;
    double mu_leading_ptot, mu_trailing_ptot;
};

#endif /* BSJPSIPHIFINDER_H_ */
