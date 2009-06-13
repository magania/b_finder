/*
 * bs_finder.cpp
 *
 *  Created on: Jun 10, 2009
 *      Author: magania
 */

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
#if defined(P17) || defined(P20)
#include "IPcal.hpp"
#endif

#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include "TLorentzRotation.h"

using namespace std;
using namespace AA;

double PDG_BS_MASS = 5.36636;
double PDG_KAON_MASS = 0.493677;
double PDG_MU_MASS = 0.1056583668;

double MASS_BS_MIN = 5.0;
double MASS_BS_MAX = 5.8;

/*
Ptl.hpp
bool pdl(PtlLst* plst, double& pdl, double& epdl);

Ptl.cpp
bool AA::Ptl::pdl(PtlLst* plst, double mass, double& pdl, double& e_pdl){
  //Exclude all particles in the list plst from the primary vertex
  //and compute the resulting decay length

  if(_pvp == 0) return false;
  if(_pve == 0) return false;

  //Create list of chains to be excluded
  list<const Chain*> lpch;
  for(PtlLstIt p = plst->begin(); p != plst->end(); ++p){
    lpch.push_back((*p)->chain());
  }
  HepVector x(3);
  HepSymMatrix vx(3);
  (const_cast<Vrt*>(_pvp))->exclude(lpch,x,vx);


  HepVector P(2);
  HepSymMatrix VP(2);
  P(1) = _mom(1); P(2) = _mom(2);
  VP(1,1) = _vmom(1,1);  VP(1,2) = _vmom(1,2);
  VP(2,2) = _vmom(2,2);

  HepVector PV(2);
  HepSymMatrix VPV(2);
  PV(1) = x(1); PV(2) = x(2);
  VPV(1,1) = vx(1,1);  VPV(1,2) = vx(1,2);
  VPV(2,2) = vx(2,2);

  HepVector SV(2);
  HepSymMatrix VSV(2);
  if(!_cnsApplied) {
     SV(1) = _pve->x(1);  SV(2) = _pve->x(2);
     VSV(1,1) = _pve->v(1,1);  VSV(1,2) = _pve->v(1,2);
     VSV(2,2) = _pve->v(2,2);
  } else {
     SV(1) = _cns(1);  SV(2) = _cns(2);
     VSV(1,1) = _vCns(1,1);  VSV(1,2) = _vCns(1,2);
     VSV(2,2) = _vCns(2,2);
  }

  HepVector L(2);
  L = SV - PV;
  double lxy = dot(L,P)/P.norm();
  pdl = lxy*mass/P.norm();

 / * Aij = PiPj/p2
   * Bij = LiLj/Lxy2 (Li = SVi - PVi)
   * EPij = Vij(P)/p2
   * ELij = Vij(L)/Lxy^2;
   * Cij = LiPj/(pLxy)
   * /

  HepSymMatrix VL(2);
  VL = VSV + VPV;

  HepSymMatrix A(2), B(2), EP(2), EL(2);
  HepMatrix C(2,2);
  double p2 = P.norm()*P.norm();
  A(1,1) = P(1)*P(1)/p2;  A(1,2) = P(1)*P(2)/p2;
  A(2,2) = P(2)*P(2)/p2;

  double lxy2 = lxy*lxy;
  B(1,1) = L(1)*L(1)/lxy2;  B(1,2) = L(1)*L(2)/lxy2;
  B(2,2) = L(2)*L(2)/lxy2;

  double plxy = P.norm()*lxy;
  C(1,1) = L(1)*P(1)/plxy;  C(1,2) = L(1)*P(2)/plxy;
  C(2,1) = L(2)*P(1)/plxy;  C(2,2) = L(2)*P(2)/plxy;

  EP = VP/p2;
  EL = VL/lxy2;

        //Computing Sigma Lxy
        // Sigma Lxy^2 = Tr{A*EL + (B + 4*A - 4*C)*E}
  HepMatrix SL(2,2);
  SL = A*EL + (B + 4*A - 4*C)*EP;
  e_pdl = fabs(pdl) * sqrt(SL(1,1)+SL(2,2));
  return true;
}

 */

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
		Double_t &w1,Double_t &w2,Double_t &w3){
  //--> Jpsi
  TLorentzVector Jpsi = mu_p; Jpsi+=mu_n;
  //--> Phi
  TLorentzVector Phi = K_p; Phi+=K_n;

  //--> Copies to be use later
  TLorentzVector Jpsi_r = Jpsi;
  TLorentzVector K_p_r = K_p;

  //--> Defining variables

  TVector3 pPhi3  = Phi.BoostVector(); // velocity of phi
  TVector3 pJpsi3  = Jpsi.BoostVector(); // velocity of phi

  Phi.Boost(-pJpsi3);
  K_p.Boost(-pJpsi3);
  mu_p.Boost(-pJpsi3);

  TVector3 xv = Phi.Vect().Unit();
  TVector3 pKpos = K_p.Vect().Unit();
  Double_t xv_dot_pKpos = xv.Dot(pKpos);
  TVector3 yv_tmp;  yv_tmp.SetXYZ(pKpos.X() - xv_dot_pKpos*xv.X(),pKpos.Y() - xv_dot_pKpos*xv.Y(),pKpos.Z() - xv_dot_pKpos*xv.Z());
  TVector3 yv = yv_tmp.Unit();
  TVector3 zv = xv.Cross(yv);
  TVector3 plpos = mu_p.Vect().Unit();

  Double_t sintheta_cospsi = plpos.Dot(xv);
  Double_t sintheta_sinpsi = plpos.Dot(yv);

  w1 = TMath::ATan2(sintheta_sinpsi,sintheta_cospsi);
  w2 = plpos.Dot(zv);

  K_p_r.Boost(-pPhi3);
  Jpsi_r.Boost(-pPhi3);
  TVector3 pKpos_r = K_p_r.Vect().Unit();
  TVector3 pJpsi_r = Jpsi_r.Vect().Unit();

  w3 = -pKpos_r.Dot(pJpsi_r);
  return;
}


void params() {
	//Parameters for Real Data and MC
	AA::debugLevel = 0;
	//  AA::firstEventIn = 6539424; // Optional first event to be processed
	AA::monitor.set("AA_ElmBox_Maps", 2, 2);
	AA::monitor.set("AA_State5_Chi2", 10, 2);
	//  AA::monitor.set("AA_DST_Field",100,1);

	AA::xBeamDefault = -0.030;
	AA::yBeamDefault = +0.030;

	AA::nTrkMinVrt = 3;
	AA::correctMomentum = false;

	//parameters for jet finding (see PYTHIA description)
	AA::mstu[46] = 5;
	AA::paru[44] = 15.;

	//definition of input blocks
	AA::blocksIn.clear();
	//definition of output blocks
	AA::blocksOut.clear();

	AA::tuneImpact = true;

	// :)
	findPVFlag = true;
	findV0Flag = true;
	findCascadeFlag = true;
	findJetsFlag = true;
	findJPsiFlag = true;
	findSVFlag = true;

	selectEleFlag = true;
}

void start_ana(int argc, char** argv) {
	if (argc != 3 || !(argv[1][0] == 'f' || argv[1][0] == 'e')) {
		std::cout << "USAGE: " << argv[0] << " f file" << std::endl;
		std::cout << "       " << argv[0] << " e elist" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string file = argv[2];

	AA::newJob();
	AA::newAnaJob();
	AA::newSelectJob();
	params();
	AA::displayParameters();

	bool ok = false;
	switch (argv[1][0]) {
	case 'e':
		std::cout << argv[0] << " II: Doing event lists " << std::endl;
		ok = AA::dst.eventLst(file);
		break;
	case 'f':
		std::cout << argv[0] << " II: Doing sequential mode" << std::endl;
		ok = AA::dst.fileLst(file);
		break;
	default:
		std::cout << argv[0] << " EE: Unknown switch " << argv[1][0]
				<< std::endl;
		exit(EXIT_FAILURE);
	}
	if (!ok) {
		std::cout << argv[0] << " EE: Input stream is not valid or is empty."
				<< std::endl;
		exit(EXIT_FAILURE);
	}
}

bool nextEvent() {
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();
	return AA::dst.getEvent();
}

int main(int argc, char** argv) {
	start_ana(argc, argv);
	/* -- We will save all in this tree --*/
	TFile root_file("bs.root", "RECREATE");
	TTree tree("tree", "all info.");

	/* -- Savers: Ricardo classes that save pretty much all :) -- */
	EvtSaver evt_saver(tree);
	TagSaver tag_saver(tree); // Yorch tag saver.
	JPsiFinder jpsi_finder(tree);
	PhiFinder phi_finder(tree);
	VrtSaver vb_saver("b_vrt", tree);
	VrtSaver vp_saver("b_pv_vrt", tree);

	/* -- Info missing by the savers -- */
	double bs_mass, bs_mass_error, bs_lhtag, bs_pdl, bs_epdl;
	double bs_iso, bs_iso_drmax, bs_iso_75;
	double bs_jpsikp_chi2, bs_jpsikm_chi2;
	double bs_angle_phi, bs_angle_ctheta, bs_angle_cpsi;

	tree.Branch("bs_mass", &bs_mass, "bs_mass/D");
	tree.Branch("bs_mass_error", &bs_mass_error, "bs_mass_error/D");
	tree.Branch("bs_lhtag", &bs_lhtag, "bs_lhtag/D");
	tree.Branch("bs_pdl", &bs_pdl, "bs_pdl/D");
	tree.Branch("bs_epdl", &bs_epdl, "bs_epdl/D");

	tree.Branch("bs_angle_phi", &bs_angle_phi, "bs_angle_phi/D");
	tree.Branch("bs_angle_ctheta", &bs_angle_ctheta, "bs_angle_ctheta/D");
	tree.Branch("bs_angle_cpsi", &bs_angle_cpsi, "bs_angle_cpsi/D");

	tree.Branch("bs_iso", &bs_iso, "bs_iso/D");
	tree.Branch("bs_iso_drmax", &bs_iso_drmax, "bs_iso_drmax/D");
	tree.Branch("bs_iso_75", &bs_iso_75, "bs_iso_75/D");

	tree.Branch("bs_jpsikp_chi2", &bs_jpsikp_chi2, "bs_jpsikp_chi2/D");
	tree.Branch("bs_jpsikm_chi2", &bs_jpsikm_chi2, "bs_jpsikm_chi2/D");

	/* -- Initilization of geometry, field and beam spot.
	 * Should be done after the fileLst or eventLst initilization -- */
	AA::det.input();
	AA::field.input();
	AA::spot.input();
#ifdef P17
	AA:ipcalib.input();
#endif
	AA::newEvent();
	AA::newAnaEvent();
	AA::newSelectEvent();

	while (nextEvent()) {
		AA::setField(); //define beam-spot and field value
		AA::spot.set();
		AA::analyse();
		AA::select(AA::TAG);

		if (!jpsi_finder.find())
			continue;
		if (!phi_finder.find())
			continue;

		jpsi_finder.begin();
		phi_finder.begin();
		while (jpsi_finder.next())
			while (phi_finder.next()) {
				Ptl* jpsi = &jpsi_finder.getJPsi();
				Ptl* muplus = &jpsi_finder.getMuPlus();
				Ptl* muminus = &jpsi_finder.getMuMinus();
				Ptl* kplus = &phi_finder.getKPlus();
				Ptl* kminus = &phi_finder.getKMinus();

				/* ---------------- Selection ------------- */
				if (muplus == kplus || muminus == kminus)
					continue;

			    /* -- b vertex -- */
			    TrkLst track_list;
			    track_list.push_back(muplus);
			    track_list.push_back(muminus);
			    track_list.push_back(kplus);
			    track_list.push_back(kminus);
			    Vrt b_vrt;
			    if (!b_vrt.fill(jpsi->decayVertex()->x(), &track_list))
			      continue;
			    if (!b_vrt.filter())
			      continue;
			    if (b_vrt.size() != 4)
			      continue;

			    /* -- jpsi + kplus vertex -- */
			    TrkLst jpsi_kp_track_list;
			    jpsi_kp_track_list.push_back(muplus);
			    jpsi_kp_track_list.push_back(muminus);
			    jpsi_kp_track_list.push_back(kplus);
			    Vrt jpsi_kp_vrt;
			    if (!jpsi_kp_vrt.fill(jpsi->decayVertex()->x(), &jpsi_kp_track_list) ||
			    		!jpsi_kp_vrt.filter() ||
			    		jpsi_kp_vrt.size() != 3)
			      continue;

			    /* -- jpsi + kminus vertex -- */
			    TrkLst jpsi_km_track_list;
			    jpsi_km_track_list.push_back(muplus);
			    jpsi_km_track_list.push_back(muminus);
			    jpsi_km_track_list.push_back(kminus);
			    Vrt jpsi_km_vrt;
			    if (!jpsi_km_vrt.fill(jpsi->decayVertex()->x(), &jpsi_km_track_list) ||
			    		!jpsi_km_vrt.filter() ||
			    		jpsi_km_vrt.size() != 3)
			        continue;

			    /* -- Construct new particle (B) from these 4 tracks -- */
			    Ptl b;
			    if (!b.combine(&b_vrt, 0))
			      continue;
			    if (!b.associateToVrt(&AA::vrtPBox))
			      continue;

			    /*-- Find primary vertex with minimal distance to B vertex --*/
			    Vrt* b_pv = 0;
			    if (!b_vrt.associateToVrt(&AA::vrtPBox, b_pv))
			      continue;

			    /* -- Guennadi Vertex Constraint gvc --*/
			    /* Define constraints */
			    CnsLst clst;
			    Constraint cns;

			    /* Constraint on J/psi mass */
			    cns.type = AA::CNS_MASS;
			    cns.lst.clear();
			    cns.lst.push_back(muplus);
			    cns.lst.push_back(muminus);
			    cns.value1 = AA::mass[AA::JPSI];
			    cns.masses.clear();
			    cns.masses.push_back(AA::mass[AA::MU_PLUS]);
			    cns.masses.push_back(AA::mass[AA::MU_MINUS]);
			    clst.push_back(cns);

			    /* Apply constraint */
			    if (!b.applyConstraints(&clst))
			      continue;

			    /* -- mass of Bs and uncertainty --*/
			    vector<PType> types(4);
			    types[0] = AA::MU_PLUS;
			    types[1] = AA::MU_MINUS;
			    types[2] = AA::K_PLUS;
			    types[3] = AA::K_MINUS;
			    double mb, vmb, emb;
			    b.mass(types, mb, vmb);
			    emb = sqrt(fabs(vmb));

			    /* -- cut on B mass -- */
			    if ( mb < MASS_BS_MIN || MASS_BS_MAX < mb)
			  	  continue;

			    /* Momentum of decay products and uncertainty
			     * (taking into account the applied constraints)
			     * The order of particles is:
			     * moms[0] - mu+
			     * moms[1] - mu-
			     * moms[2] - K+
			     * moms[3] - K-
			     */
			    vector<HepVector> moms;
			    vector<HepSymMatrix> vmoms;
			    if(!b.momChildren(moms,vmoms)) continue;
				/* ---------------- /Selection ------------- */
			    /* Selection generate this variables:
			     * Vrt jpsi_kp_vrt        JPsi + K+ vertex
			     * Vrt jpsi_km_vrt        JPsi + K- vertex
			     * Vrt b_vrt              Bs vertex
			     * Vrt* b_pv              Bs primary vertev
			     * Ptl b                  Bs particle
			     * double mb, vmb, emb    Bs mass, sigma(Bs mass)^2 , sigma(Bs_mass)
			     *
			     * clst,cns,types,track_list, track_kp_list, track_km_list
			     */

				/* --------------------- savers -------------------------*/
				evt_saver.fill();       // Run, Evt, Triggers.
				jpsi_finder.fill();     // Mu+, Mu-, Jpsi info.
				phi_finder.fill();      // K+, K-, Phi info.
				vb_saver.fill(b_vrt);  // B vertex info.
				vp_saver.fill(*b_pv);    // Primary vertex info.

				PtlLst particle_list;
				particle_list.push_back(muplus);
				particle_list.push_back(muminus);
				particle_list.push_back(kplus);
				particle_list.push_back(kminus);
				tag_saver.fill(b.mom(), b_pv, &particle_list);  // Tagging info.

				/* -- Info missing by the savers -- */
				bs_mass = mb;
				bs_mass_error = emb;

				/* -- lifetime & lifetime_error -- */
				b.pdl(&particle_list, PDG_BS_MASS, bs_pdl, bs_epdl);

				/* -- Calculate isolation of B -- */
				const AA::PtlLst* ptl_lst = AA::ptlBox.particles();
				double drmax = muminus->dR(b.mom());
				if ( muplus->dR(b.mom()) > drmax ) drmax = muplus->dR(b.mom());
				if ( kplus->dR(b.mom())  > drmax )  drmax = kplus->dR(b.mom());
				if ( kminus->dR(b.mom()) > drmax ) drmax = kminus->dR(b.mom());
				double sum = 0; double sum_drmax = 0; double sum_75 = 0;
				for (PtlLstCIt piso = ptl_lst->begin(); piso != ptl_lst->end(); ++piso) {
					Ptl* ptlIso = *piso;
					if (ptlIso == muplus || ptlIso == muminus || ptlIso == kplus || ptlIso == kminus)
						continue;
					double driso = ptlIso->dR(b.mom());
					if (driso <= 0.5)   sum += ptlIso->ptot();
					if (driso <= drmax) sum_drmax += ptlIso->ptot();
					if (driso <= 0.75)  sum_75 += ptlIso->ptot();
				}
				bs_iso       = b.ptot() / (b.ptot() + sum);
				bs_iso_drmax = b.ptot() / (b.ptot() + sum_drmax);
				bs_iso_75    = b.ptot() / (b.ptot() + sum_75);

				/* --  likelihood ratio --
				* Chi2 Difference between B vertex and J/psi vertex
				*/
				double chi2_b_minus_jpsi = b_vrt.chi2() - jpsi->decayVertex()->chi2();

				/* Define likelihood ratio discriminating variables */
				vector<double> lh_var;
				vector<double> lh_rat;
				lh_var.push_back(bs_iso);
				lh_var.push_back(chi2_b_minus_jpsi);
				lh_var.push_back(b.pt());
				lh_var.push_back(phi_finder.getMass());

				AA::combinedTagJPsiPhi(lh_var, lh_rat, bs_lhtag);

				/* -- Intermediate vertex chi2 -- */
				bs_jpsikp_chi2 = jpsi_kp_vrt.chi2();
				bs_jpsikm_chi2 = jpsi_km_vrt.chi2();

				/* -- Transversity Angles -- */
				TLorentzVector l_kplus, l_kminus, l_muplus, l_muminus;
				l_muplus.SetXYZM (moms[0](1), moms[0](2), moms[0](3), PDG_MU_MASS);
				l_muminus.SetXYZM(moms[1](1), moms[1](2), moms[1](3), PDG_MU_MASS);
				l_kplus.SetXYZM  (moms[2](1), moms[2](2), moms[2](3), PDG_KAON_MASS);
				l_kminus.SetXYZM (moms[3](1), moms[3](2), moms[3](3), PDG_KAON_MASS);

				threeAngles(l_muplus, l_muminus, l_kplus, l_kminus, bs_angle_phi, bs_angle_ctheta, bs_angle_cpsi);
				/* --------------------- /savers ------------------------*/
				tree.Fill();   // FILL
			}
	}//End while next event.

	tree.Write();
	root_file.Write();
	root_file.Close();
	std::cout << argv[0] << " II: bs_finder ended succesfully." << std::endl;
	std::exit(EXIT_SUCCESS);
}
