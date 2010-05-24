/*
 * BdJPsiKstarFinder.cpp
 *
 *  Created on: Jul 16, 2009
 *      Author: magania
 */

#include "BdJPsiKstarFinder.h"

BdJPsiKstarFinder::BdJPsiKstarFinder(JPsiFinder *jpsi, KstarFinder *kstar, TTree &tree, BdJPsiKstarMCFinder *mc) :
	tag_saver(tree),
	bd_saver("bd",tree, false, false),
	vb_saver("bd_vrt", tree),
	vp_saver("bd_pv_vrt", tree)
{
	index = -1;

	jpsi_finder = jpsi;
	kstar_finder = kstar;
	mc_finder = mc;

	tree.Branch("bd_mass", &bd_mass, "bd_mass/D");
	tree.Branch("bd_mass_error", &bd_mass_error, "bd_mass_error/D");
	tree.Branch("bd_mass_sw", &bd_mass_sw, "bd_mass_sw/D");
	tree.Branch("bd_mass_error_sw", &bd_mass_error_sw, "bd_mass_error_sw/D");
	tree.Branch("bd_pdl", &bd_pdl, "bd_pdl/D");
	tree.Branch("bd_epdl", &bd_epdl, "bd_epdl/D");

	tree.Branch("mu_plus_dR", &mu_plus_dR, "mu_plus_dR/D");
    tree.Branch("mu_minus_dR", &mu_minus_dR, "mu_minus_dR/D");
    tree.Branch("kaon_dR", &kaon_dR, "kaon_dR/D");
    tree.Branch("pion_dR", &pion_dR, "pion_dR/D");

	tree.Branch("mu_plus_cpx", &mu_plus_cpx, "mu_plus_cpx/D");
	tree.Branch("mu_plus_cpy", &mu_plus_cpy, "mu_plus_cpy/D");
	tree.Branch("mu_plus_cpz", &mu_plus_cpz, "mu_plus_cpz/D");
	tree.Branch("mu_minus_cpx", &mu_minus_cpx, "mu_minus_cpx/D");
	tree.Branch("mu_minus_cpy", &mu_minus_cpy, "mu_minus_cpy/D");
	tree.Branch("mu_minus_cpz", &mu_minus_cpz, "mu_minus_cpz/D");
	tree.Branch("kaon_cpx", &kaon_cpx, "kaon_cpx/D");
	tree.Branch("kaon_cpy", &kaon_cpy, "kaon_cpy/D");
	tree.Branch("kaon_cpz", &kaon_cpz, "kaon_cpz/D");
	tree.Branch("pion_cpx", &pion_cpx, "pion_cpx/D");
	tree.Branch("pion_cpy", &pion_cpy, "pion_cpy/D");
	tree.Branch("pion_cpz", &pion_cpz, "pion_cpz/D");

	tree.Branch("bd_angle_phi", &bd_angle_phi, "bd_angle_phi/D");
	tree.Branch("bd_angle_ctheta", &bd_angle_ctheta, "bd_angle_ctheta/D");
	tree.Branch("bd_angle_cpsi", &bd_angle_cpsi, "bd_angle_cpsi/D");
	
        tree.Branch("bd_angle_phi_sw", &bd_angle_phi_sw, "bd_angle_phi_sw/D");
	tree.Branch("bd_angle_ctheta_sw", &bd_angle_ctheta_sw, "bd_angle_ctheta_sw/D");
	tree.Branch("bd_angle_cpsi_sw", &bd_angle_cpsi_sw, "bd_angle_cpsi_sw/D");

	tree.Branch("bd_iso", &bd_iso, "bd_iso/D");
	tree.Branch("bd_iso_drmax", &bd_iso_drmax, "bd_iso_drmax/D");
	tree.Branch("bd_iso_75", &bd_iso_75, "bd_iso_75/D");
    tree.Branch("bd_iso_pv", &bd_iso_pv, "bd_iso_pv/D");
    tree.Branch("bd_iso_drmax_pv", &bd_iso_drmax_pv, "bd_iso_drmax_pv/D");
    tree.Branch("bd_iso_75_pv", &bd_iso_75_pv, "bd_iso_75_pv/D");

	tree.Branch("bd_jpsikn_chi2", &bd_jpsikn_chi2, "bd_jpsikn_chi2/D");
	tree.Branch("bd_jpsipi_chi2", &bd_jpsipi_chi2, "bd_jpsipi_chi2/D");

	tree.Branch("kstar_mass_corrected", &kstar_mass_corrected, "kstar_mass_corrected/D");
	tree.Branch("kstar_mass_corrected_error", &kstar_mass_corrected_error, "kstar_mass_corrected_error/D");
	tree.Branch("kstar_mass_corrected_sw", &kstar_mass_corrected_sw, "kstar_mass_corrected_sw/D");
	tree.Branch("kstar_mass_corrected_error_sw", &kstar_mass_corrected_error_sw, "kstar_mass_corrected_error_sw/D");

    tree.Branch("bd_ucpt",      &bd_ucpt,      "bd_ucpt/D");
    tree.Branch("bd_ucptot",      &bd_ucptot,      "bd_ucptot/D");
    tree.Branch("mu_plus_cpt",  &mu_plus_cpt,  "mu_plus_cpt/D");
    tree.Branch("mu_minus_cpt", &mu_minus_cpt, "mu_minus_cpt/D");
    tree.Branch("kaon_cpt",   &kaon_cpt,   "kaon_cpt/D");
    tree.Branch("pion_cpt",  &pion_cpt,  "pion_cpt/D");
    tree.Branch("jpsi_cpt",     &jpsi_cpt,     "jpsi_cpt/D");
    tree.Branch("kstar_cpt",      &kstar_cpt,      "kstar_cpt/D");
}

BdJPsiKstarFinder::~BdJPsiKstarFinder() {
}

void BdJPsiKstarFinder::clean(){
	_boxp.clear();
	_boxv.clear();

    v_jpsi_index.clear();
	v_kstar_index.clear();
    v_bd.clear();
    v_jpsi.clear();
	v_muplus.clear();
	v_muminus.clear();
	v_kaon.clear();
	v_pion.clear();
    v_bd_vrt.clear();
	v_jpsi_kn_vrt.clear();
	v_jpsi_pi_vrt.clear();
	v_bd_pv.clear();
	v_mb.clear();
	v_emb.clear();
    v_muplus_cpx.clear();
	v_muplus_cpy.clear();
	v_muplus_cpz.clear();
    v_muminus_cpx.clear();
	v_muminus_cpy.clear();
	v_muminus_cpz.clear();
    v_kaon_cpx.clear();
	v_kaon_cpy.clear();
	v_kaon_cpz.clear();
    v_pion_cpx.clear();
	v_pion_cpy.clear();
	v_pion_cpz.clear();
    v_bd_mom.clear();
}

int BdJPsiKstarFinder::find(){
    //std::cout << "Finding Bd .. " << std::endl;
	clean();
	jpsi_finder->begin();
	kstar_finder->begin();
	while (jpsi_finder->next())
		while (kstar_finder->next()) {
			Ptl* jpsi = &jpsi_finder->getJPsi();
			Ptl* muplus = &jpsi_finder->getMuPlus();
			Ptl* muminus = &jpsi_finder->getMuMinus();
			Ptl* kaon = &kstar_finder->getKaon();
			Ptl* pion = &kstar_finder->getPion();

			/* ---------------- Selection ------------- */
			if (muplus == kaon || muplus == pion || muminus == kaon || muminus == pion)
				continue;

		    /* -- b vertex -- */
		    TrkLst track_list;
		    track_list.push_back(muplus);
		    track_list.push_back(muminus);
		    track_list.push_back(kaon);
		    track_list.push_back(pion);
		    Vrt *bd_vrt = _boxv.newVrt();
		    if (!bd_vrt->fill(jpsi->decayVertex()->x(), &track_list))
		      continue;
		    if (!bd_vrt->filter())
		      continue;
		    if (bd_vrt->size() != 4)
		      continue;

		    /* -- jpsi + kaon vertex -- */
		    TrkLst jpsi_kn_track_list;
		    jpsi_kn_track_list.push_back(muplus);
		    jpsi_kn_track_list.push_back(muminus);
		    jpsi_kn_track_list.push_back(kaon);
		    Vrt* jpsi_kn_vrt = _boxv.newVrt();
		    if (!jpsi_kn_vrt->fill(jpsi->decayVertex()->x(), &jpsi_kn_track_list) ||
		    		!jpsi_kn_vrt->filter() ||
		    		jpsi_kn_vrt->size() != 3)
		      continue;
			    /* -- jpsi + pion vertex -- */
		    TrkLst jpsi_pi_track_list;
		    jpsi_pi_track_list.push_back(muplus);
		    jpsi_pi_track_list.push_back(muminus);
		    jpsi_pi_track_list.push_back(pion);
		    Vrt* jpsi_pi_vrt = _boxv.newVrt();
		    if (!jpsi_pi_vrt->fill(jpsi->decayVertex()->x(), &jpsi_pi_track_list) ||
		    		!jpsi_pi_vrt->filter() ||
		    		jpsi_pi_vrt->size() != 3)
		        continue;

		    /* -- Construct new particle (B) from these 4 tracks -- */
		    Ptl* bd = _boxp.newPtl();
		    if (!bd->combine(bd_vrt, 0))
		      continue;
		    if (!bd->associateToVrt(&AA::vrtPBox))
		      continue;

            /* -- Save the original bd pt -- */
		    HepVector tmp_bdmom = bd->mom();

		    /*-- Find primary vertex with minimal distance to B vertex --*/
		    Vrt* bd_pv = 0;
		    if (!bd_vrt->associateToVrt(&AA::vrtPBox, bd_pv))
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
		    if (!bd->applyConstraints(&clst))
		      continue;

		    /* -- mass of Bs and uncertainty --*/
		    vector<PType> types(4);
		    types[0] = AA::MU_PLUS;
		    types[1] = AA::MU_MINUS;
		    types[2] = AA::K_PLUS;
		    types[3] = AA::PI_PLUS;
		    double mb, vmb, emb;
		    bd->mass(types, mb, vmb);
		    emb = sqrt(fabs(vmb));

		    /* -- cut on B mass -- */
		    if ( mb < MASS_BD_MIN || MASS_BD_MAX < mb)
		  	  continue;

		    /* Momentum of decay products and uncertainty
		     * (taking into account the applied constraints)
		     * The order of particles is:
		     * moms[0] - mu+
		     * moms[1] - mu-
		     * moms[2] - kaon
		     * moms[3] - pi
		     */
		    vector<HepVector> moms;
		    vector<HepSymMatrix> vmoms;
		    if(!bd->momChildren(moms,vmoms)) continue;

            /* ---------------- /Selection ------------- */
		    /* Selection generate this variables:
		     * vector<Ptl*> v_bd, v_jpsi, v_muplus, v_muminus, v_kaon, v_pion;
			 * vector<Vrt*> v_b_vrt, v_jpsi_kp_vrt, v_jpsi_km_vrt, v_bd_pv;
             * vector<double> v_lhtag, v_mb, v_emb;
             * vector<double> v_muplus_cpx, v_muplus_cpy, v_muplus_cpz;
             * vector<double> v_muminus_cpx, v_muminus_cpy, v_muminus_cpz;
             * vector<double> v_kaon_cpx, v_kaon_cpy, v_kaon_cpz;
             * vector<double> v_pion_cpx, v_pion_cpy, v_pion_cpz;
		     */

            //std::cout << "Bs: " << v_bd.size()
            //          << " JPsi: " << jpsi_finder->getIndex()
            //          << " Kstar:" << kstar_finder->getIndex() << std::endl;
		    v_jpsi_index.push_back(jpsi_finder->getIndex());
		    v_kstar_index.push_back(kstar_finder->getIndex());

		    v_bd.push_back(bd);
		    v_jpsi.push_back(jpsi);
		    v_muplus.push_back(muplus);
		    v_muminus.push_back(muminus);
		    v_kaon.push_back(kaon);
		    v_pion.push_back(pion);
		    v_bd_vrt.push_back(bd_vrt);
		    v_jpsi_kn_vrt.push_back(jpsi_kn_vrt);
		    v_jpsi_pi_vrt.push_back(jpsi_pi_vrt);
		    v_bd_pv.push_back(bd_pv);
		    v_mb.push_back(mb);
		    v_emb.push_back(emb);
            v_bd_mom.push_back(tmp_bdmom);

		    v_muplus_cpx.push_back(moms[0](1)); v_muplus_cpy.push_back(moms[0](2)); v_muplus_cpz.push_back(moms[0](3));
		    v_muminus_cpx.push_back(moms[1](1));v_muminus_cpy.push_back(moms[1](2));v_muminus_cpz.push_back(moms[1](3));
		    v_kaon_cpx.push_back(moms[2](1));  v_kaon_cpy.push_back(moms[2](2));  v_kaon_cpz.push_back(moms[2](3));
		    v_pion_cpx.push_back(moms[3](1)); v_pion_cpy.push_back(moms[3](2)); v_pion_cpz.push_back(moms[3](3));
		}
  //std::cout << "Found: " << v_bd.size() << std::endl;
  return v_bd.size();
}

void BdJPsiKstarFinder::begin(){
	index = -1;
}

bool BdJPsiKstarFinder::next(){
	if (++index < v_bd.size())
		return true;
	else
		index=-1;
	return false;
}


void BdJPsiKstarFinder::fill(){
    //std::cout << "Filling: " << index << " JPsi:" << v_jpsi_index[index] << " Kstar:" << v_kstar_index[index] << std::endl;
	jpsi_finder->setIndex(v_jpsi_index[index]);
	kstar_finder->setIndex(v_kstar_index[index]);

	jpsi_finder->fill();     // Mu+, Mu-, Jpsi info.
	kstar_finder->fill();      // Kaon, Pion, Kstar info.
	bd_saver.fill(*v_bd[index]);  // Bs info.
	vb_saver.fill(*v_bd_vrt[index]);  // B vertex info.
	vp_saver.fill(*v_bd_pv[index]);    // Primary vertex info.
if (mc_finder)
	mc_finder->fill(v_muplus[index], v_muminus[index], v_kaon[index], v_pion[index]);

	PtlLst particle_list;
	particle_list.push_back(v_muplus[index]);
	particle_list.push_back(v_muminus[index]);
	particle_list.push_back(v_kaon[index]);
	particle_list.push_back(v_pion[index]);
	tag_saver.fill(v_bd[index]->mom(), v_bd_pv[index], &particle_list);  // Tagging info.

	/* -- Info missing by the savers -- */
	bd_mass = v_mb[index];
	bd_mass_error = v_emb[index];

        /* -- mass of Bs and uncertainty --*/
         vector<PType> types(4);
         types[0] = AA::MU_PLUS;
         types[1] = AA::MU_MINUS;
         types[2] = AA::PI_PLUS;
         types[3] = AA::K_PLUS;
         v_bd[index]->mass(types, bd_mass_sw, bd_mass_error_sw);
         bd_mass_error_sw = sqrt(fabs(bd_mass_error_sw));

    /* -- dR -- */
    mu_plus_dR = v_muplus[index]->dR(v_bd_mom[index]);
    mu_minus_dR =  v_muminus[index]->dR(v_bd_mom[index]);
    kaon_dR = v_kaon[index]->dR(v_bd_mom[index]);
    pion_dR =  v_pion[index]->dR(v_bd_mom[index]);

	/* -- Corrected momentums -- */
	mu_plus_cpx  = v_muplus_cpx[index];  mu_plus_cpy  = v_muplus_cpy[index];  mu_plus_cpz  = v_muplus_cpz[index];
	mu_minus_cpx = v_muminus_cpx[index]; mu_minus_cpy = v_muminus_cpy[index]; mu_minus_cpz = v_muminus_cpz[index];
	kaon_cpx   = v_kaon_cpx[index];   kaon_cpy   = v_kaon_cpy[index];   kaon_cpz   = v_kaon_cpz[index];
	pion_cpx  = v_pion_cpx[index];  pion_cpy  = v_pion_cpy[index];  pion_cpz  = v_pion_cpz[index];

	/* -- lifetime & lifetime_error -- */
    double ctau, vctau;
	v_bd[index]->decayLengthProper(PDG_BD_MASS, ctau, vctau, &particle_list);
	bd_pdl = ctau; bd_epdl=sqrt(fabs(vctau));

    /* -- uncorrected momentum -- */
    bd_ucpt = sqrt( v_bd_mom[index](1)*v_bd_mom[index](1) + v_bd_mom[index](2)*v_bd_mom[index](2) );
    bd_ucptot =  sqrt( v_bd_mom[index](1)*v_bd_mom[index](1) + v_bd_mom[index](2)*v_bd_mom[index](2)
                         + v_bd_mom[index](3)*v_bd_mom[index](3));

	/* -- Calculate isolation of B -- */
	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();
	double drmax = v_muminus[index]->dR(v_bd_mom[index]);
	if ( v_muplus[index]->dR(v_bd_mom[index]) > drmax ) drmax = v_muplus[index]->dR(v_bd_mom[index]);
	if ( v_kaon[index]->dR(v_bd_mom[index])  > drmax ) drmax = v_kaon[index]->dR(v_bd_mom[index]);
	if ( v_pion[index]->dR(v_bd_mom[index]) > drmax ) drmax = v_pion[index]->dR(v_bd_mom[index]);
	double sum = 0; double sum_drmax = 0; double sum_75 = 0;
    double sum_pv = 0; double sum_drmax_pv = 0; double sum_75_pv = 0;
	for (PtlLstCIt piso = ptl_lst->begin(); piso != ptl_lst->end(); ++piso) {
		Ptl* ptlIso = *piso;
		if (ptlIso == v_muplus[index] || ptlIso == v_muminus[index] || ptlIso == v_kaon[index] || ptlIso == v_pion[index])
			continue;
		double driso = ptlIso->dR(v_bd_mom[index]);
		if (driso <= 0.5)   sum += ptlIso->ptot();
		if (driso <= drmax) sum_drmax += ptlIso->ptot();
		if (driso <= 0.75)  sum_75 += ptlIso->ptot();
        if (v_jpsi[index]->primaryVertex() != ptlIso->primaryVertex()) continue;
        if (driso <= 0.5)   sum_pv += ptlIso->ptot();
        if (driso <= drmax) sum_drmax_pv += ptlIso->ptot();
        if (driso <= 0.75)  sum_75_pv += ptlIso->ptot();
	}
	bd_iso       = bd_ucptot / (bd_ucptot + sum);
	bd_iso_drmax = bd_ucptot / (bd_ucptot + sum_drmax);
	bd_iso_75    = bd_ucptot / (bd_ucptot + sum_75);
    bd_iso_pv       = bd_ucptot / (bd_ucptot + sum_pv);
    bd_iso_drmax_pv = bd_ucptot / (bd_ucptot + sum_drmax_pv);
    bd_iso_75_pv    = bd_ucptot / (bd_ucptot + sum_75_pv);

	/* -- Intermediate vertex chi2 -- */
	bd_jpsikn_chi2 = v_jpsi_kn_vrt[index]->chi2();
	bd_jpsipi_chi2 = v_jpsi_pi_vrt[index]->chi2();

	/* -- Transversity Angles -- */
	TLorentzVector l_kaon, l_pion, l_muplus, l_muminus, l_kaon_sw, l_pion_sw;
	l_muplus.SetXYZM (v_muplus_cpx[index], v_muplus_cpy[index], v_muplus_cpz[index], PDG_MU_MASS);
	l_muminus.SetXYZM(v_muminus_cpx[index], v_muminus_cpy[index], v_muminus_cpz[index], PDG_MU_MASS);
	l_kaon.SetXYZM  (v_kaon_cpx[index], v_kaon_cpy[index], v_kaon_cpz[index], PDG_KAON_MASS);
	l_pion.SetXYZM (v_pion_cpx[index], v_pion_cpy[index], v_pion_cpz[index], PDG_PION_MASS);
	l_kaon_sw.SetXYZM  (v_kaon_cpx[index], v_kaon_cpy[index], v_kaon_cpz[index], PDG_PION_MASS);
	l_pion_sw.SetXYZM (v_pion_cpx[index], v_pion_cpy[index], v_pion_cpz[index], PDG_KAON_MASS);

	threeAngles(l_muplus, l_muminus, l_kaon, l_pion, bd_angle_phi, bd_angle_ctheta, bd_angle_cpsi);
	threeAngles(l_muplus, l_muminus, l_pion_sw, l_kaon_sw, bd_angle_phi_sw, bd_angle_ctheta_sw, bd_angle_cpsi_sw);
        
	/* -- Kstar corrected mass -- */
	vector<PType> kstar_types(2);
	vector<PType> kstar_types_sw(2);
	PtlLst kstar_list;
	kstar_list.push_back(v_kaon[index]);
	kstar_list.push_back(v_pion[index]);
	kstar_types[0] = AA::K_PLUS;
	kstar_types[1] = AA::PI_MINUS;
	kstar_types_sw[0] = AA::PI_PLUS;
	kstar_types_sw[1] = AA::K_MINUS;
	v_bd[index]->mass(kstar_types,kstar_mass_corrected,kstar_mass_corrected_error,&kstar_list);
	v_bd[index]->mass(kstar_types_sw,kstar_mass_corrected_sw,kstar_mass_corrected_error_sw,&kstar_list);
	kstar_mass_corrected_error = sqrt(fabs(kstar_mass_corrected_error));
	kstar_mass_corrected_error_sw = sqrt(fabs(kstar_mass_corrected_error_sw));

	/* -- corrected pt and bd uncorrected pt -- */
    mu_plus_cpt = sqrt(v_muplus_cpx[index]*v_muplus_cpx[index] + v_muplus_cpy[index]*v_muplus_cpy[index]);
    mu_minus_cpt = sqrt(v_muminus_cpx[index]*v_muminus_cpx[index] + v_muminus_cpy[index]*v_muminus_cpy[index]);
    kaon_cpt = sqrt(v_kaon_cpx[index]*v_kaon_cpx[index] + v_kaon_cpy[index]*v_kaon_cpy[index]);
    pion_cpt = sqrt(v_pion_cpx[index]*v_pion_cpx[index] + v_pion_cpy[index]*v_pion_cpy[index]);
    jpsi_cpt = sqrt( (v_muplus_cpx[index]+v_muminus_cpx[index])*(v_muplus_cpx[index]+v_muminus_cpx[index])
                     + (v_muplus_cpy[index]+v_muminus_cpy[index])*(v_muplus_cpy[index]+v_muminus_cpy[index]) );
    kstar_cpt  = sqrt( (v_kaon_cpx[index]+v_pion_cpx[index])*(v_kaon_cpx[index]+v_pion_cpx[index])
                     + (v_kaon_cpy[index]+v_pion_cpy[index])*(v_kaon_cpy[index]+v_pion_cpy[index]) );
}

/* -- Legendary threeAngles function:
 * mu_p: positive muon
 * mu_n: negative muon
 * K_p : positive Kaon
 * K_n : negative Kaon
 *
 * w1 = phi, w2 = cos(theta), w3 = cos(psi)
 *
 */
void BdJPsiKstarFinder::threeAngles(TLorentzVector mu_p,
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
