/*
 * BsJPsiFinder.cpp
 *
 *  Created on: Jun 22, 2009
 *      Author: magania
 */

#include "BsJPsiPhiFinder.h"

BsJPsiPhiFinder::BsJPsiPhiFinder(JPsiFinder *jpsi, PhiFinder *phi, TTree &tree, BsJPsiPhiMCFinder *mc) :
	tag_saver(tree),
	bs_saver("bs",tree, false, false),
	vb_saver("bs_vrt", tree),
	vp_saver("bs_pv_vrt", tree)
{
	index = -1;

	jpsi_finder = jpsi;
	phi_finder = phi;
	mc_finder = mc;

	tree.Branch("bs_mass", &bs_mass, "bs_mass/D");
	tree.Branch("bs_mass_error", &bs_mass_error, "bs_mass_error/D");
	tree.Branch("bs_lhtag", &bs_lhtag, "bs_lhtag/D");
	tree.Branch("bs_pdl", &bs_pdl, "bs_pdl/D");
	tree.Branch("bs_epdl", &bs_epdl, "bs_epdl/D");

	tree.Branch("mu_plus_cpx", &mu_plus_cpx, "mu_plus_cpx/D");
	tree.Branch("mu_plus_cpy", &mu_plus_cpy, "mu_plus_cpy/D");
	tree.Branch("mu_plus_cpz", &mu_plus_cpz, "mu_plus_cpz/D");
	tree.Branch("mu_minus_cpx", &mu_minus_cpx, "mu_minus_cpx/D");
	tree.Branch("mu_minus_cpy", &mu_minus_cpy, "mu_minus_cpy/D");
	tree.Branch("mu_minus_cpz", &mu_minus_cpz, "mu_minus_cpz/D");
	tree.Branch("k_plus_cpx", &k_plus_cpx, "k_plus_cpx/D");
	tree.Branch("k_plus_cpy", &k_plus_cpy, "k_plus_cpy/D");
	tree.Branch("k_plus_cpz", &k_plus_cpz, "k_plus_cpz/D");
	tree.Branch("k_minus_cpx", &k_minus_cpx, "k_minus_cpx/D");
	tree.Branch("k_minus_cpy", &k_minus_cpy, "k_minus_cpy/D");
	tree.Branch("k_minus_cpz", &k_minus_cpz, "k_minus_cpz/D");

	tree.Branch("bs_angle_phi", &bs_angle_phi, "bs_angle_phi/D");
	tree.Branch("bs_angle_ctheta", &bs_angle_ctheta, "bs_angle_ctheta/D");
	tree.Branch("bs_angle_cpsi", &bs_angle_cpsi, "bs_angle_cpsi/D");

	tree.Branch("bs_iso", &bs_iso, "bs_iso/D");
	tree.Branch("bs_iso_drmax", &bs_iso_drmax, "bs_iso_drmax/D");
	tree.Branch("bs_iso_75", &bs_iso_75, "bs_iso_75/D");
    tree.Branch("bs_iso_pv", &bs_iso_pv, "bs_iso_pv/D");
    tree.Branch("bs_iso_drmax_pv", &bs_iso_drmax_pv, "bs_iso_drmax_pv/D");
    tree.Branch("bs_iso_75_pv", &bs_iso_75_pv, "bs_iso_75_pv/D");

	tree.Branch("bs_jpsikp_chi2", &bs_jpsikp_chi2, "bs_jpsikp_chi2/D");
	tree.Branch("bs_jpsikm_chi2", &bs_jpsikm_chi2, "bs_jpsikm_chi2/D");

	tree.Branch("phi_mass_corrected", &phi_mass_corrected, "phi_mass_corrected/D");
	tree.Branch("phi_mass_corrected_error", &phi_mass_corrected_error, "phi_mass_corrected_error/D");

        tree.Branch("bs_ucpt",      &bs_ucpt,      "bs_ucpt/D");
        tree.Branch("mu_plus_cpt",  &mu_plus_cpt,  "mu_plus_cpt/D");
        tree.Branch("mu_minus_cpt", &mu_minus_cpt, "mu_minus_cpt/D");
        tree.Branch("k_plus_cpt",   &k_plus_cpt,   "k_plus_cpt/D");
        tree.Branch("k_minus_cpt",  &k_minus_cpt,  "k_minus_cpt/D");
        tree.Branch("jpsi_cpt",     &jpsi_cpt,     "jpsi_cpt/D");
        tree.Branch("phi_cpt",      &phi_cpt,      "phi_cpt/D");
}

BsJPsiPhiFinder::~BsJPsiPhiFinder() {
	// TODO Auto-generated destructor stub
}

void BsJPsiPhiFinder::clean(){
        v_jpsi_index.clear();
	v_phi_index.clear();
        v_bs.clear();
        v_jpsi.clear();
	v_muplus.clear();
	v_muminus.clear();
	v_kplus.clear();
	v_kminus.clear();
        v_bs_vrt.clear();
	v_jpsi_kp_vrt.clear();
	v_jpsi_km_vrt.clear();
	v_bs_pv.clear();
        v_lhtag.clear();
	v_mb.clear();
	v_emb.clear();
        v_muplus_cpx.clear();
	v_muplus_cpy.clear();
	v_muplus_cpz.clear();
        v_muminus_cpx.clear();
	v_muminus_cpy.clear();
	v_muminus_cpz.clear();
        v_kplus_cpx.clear();
	v_kplus_cpy.clear();
	v_kplus_cpz.clear();
        v_kminus_cpx.clear();
	v_kminus_cpy.clear();
	v_kminus_cpz.clear();
}

int BsJPsiPhiFinder::find(){
        //std::cout << "Finding Bs .. " << std::endl;
	clean();
	jpsi_finder->begin();
	phi_finder->begin();
	while (jpsi_finder->next())
		while (phi_finder->next()) {
			Ptl* jpsi = &jpsi_finder->getJPsi();
			Ptl* muplus = &jpsi_finder->getMuPlus();
			Ptl* muminus = &jpsi_finder->getMuMinus();
			Ptl* kplus = &phi_finder->getKPlus();
			Ptl* kminus = &phi_finder->getKMinus();

			/* ---------------- Selection ------------- */
			if (muplus == kplus || muminus == kminus)
				continue;

		    /* -- b vertex -- */
		    TrkLst track_list;
		    track_list.push_back(muplus);
		    track_list.push_back(muminus);
		    track_list.push_back(kplus);
		    track_list.push_back(kminus);
		    Vrt *bs_vrt = new Vrt();  //***
		    if (!bs_vrt->fill(jpsi->decayVertex()->x(), &track_list))
		      continue;
		    if (!bs_vrt->filter())
		      continue;
		    if (bs_vrt->size() != 4)
		      continue;

		    /* -- jpsi + kplus vertex -- */
		    TrkLst jpsi_kp_track_list;
		    jpsi_kp_track_list.push_back(muplus);
		    jpsi_kp_track_list.push_back(muminus);
		    jpsi_kp_track_list.push_back(kplus);
		    Vrt* jpsi_kp_vrt = new Vrt(); //***
		    if (!jpsi_kp_vrt->fill(jpsi->decayVertex()->x(), &jpsi_kp_track_list) ||
		    		!jpsi_kp_vrt->filter() ||
		    		jpsi_kp_vrt->size() != 3)
		      continue;
			    /* -- jpsi + kminus vertex -- */
		    TrkLst jpsi_km_track_list;
		    jpsi_km_track_list.push_back(muplus);
		    jpsi_km_track_list.push_back(muminus);
		    jpsi_km_track_list.push_back(kminus);
		    Vrt* jpsi_km_vrt = new Vrt(); //***
		    if (!jpsi_km_vrt->fill(jpsi->decayVertex()->x(), &jpsi_km_track_list) ||
		    		!jpsi_km_vrt->filter() ||
		    		jpsi_km_vrt->size() != 3)
		        continue;

		    /* -- Construct new particle (B) from these 4 tracks -- */
		    Ptl* bs = new Ptl(); //***
		    if (!bs->combine(bs_vrt, 0))
		      continue;
		    if (!bs->associateToVrt(&AA::vrtPBox))
		      continue;
                      
                    /* -- Save the original bs pt -- */
                    double tmp_ucpt = bs->pt();

		    /*-- Find primary vertex with minimal distance to B vertex --*/
		    Vrt* bs_pv = 0;
		    if (!bs_vrt->associateToVrt(&AA::vrtPBox, bs_pv))
			      continue;

		    /* --  likelihood ratio --
			 * Chi2 Difference between B vertex and J/psi vertex
			 */
		    double sum_iso_lh = 0;
			const AA::PtlLst* ptl_lst = AA::ptlBox.particles();
		    for( PtlLstCIt piso = ptl_lst->begin(); piso != ptl_lst->end(); ++piso){
		    	Ptl *p = *piso;
		        if(jpsi->primaryVertex() != p->primaryVertex()) continue;
		        if(p == muplus || p == muminus || p==kplus || p==kminus) continue;
		        if(p->dR(bs->mom()) > 0.5) continue;
		        sum_iso_lh += p->ptot();
		    }
	        double bs_iso_lh = bs->ptot()/(bs->ptot() + sum_iso_lh);

	        double chi2_b_minus_jpsi = bs_vrt->chi2() - jpsi->decayVertex()->chi2();
			/* Define likelihood ratio discriminating variables */
			vector<double> lh_var;
			vector<double> lh_rat;
			lh_var.push_back(bs_iso_lh);
			lh_var.push_back(chi2_b_minus_jpsi);
			lh_var.push_back(bs->pt());
			lh_var.push_back(phi_finder->getMass());

			double lhtag;
			AA::combinedTagJPsiPhi(lh_var, lh_rat, lhtag);

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
		    if (!bs->applyConstraints(&clst))
		      continue;

		    /* -- mass of Bs and uncertainty --*/
		    vector<PType> types(4);
		    types[0] = AA::MU_PLUS;
		    types[1] = AA::MU_MINUS;
		    types[2] = AA::K_PLUS;
		    types[3] = AA::K_MINUS;
		    double mb, vmb, emb;
		    bs->mass(types, mb, vmb);
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
		    if(!bs->momChildren(moms,vmoms)) continue;

            /* ---------------- /Selection ------------- */
		    /* Selection generate this variables:
		     * vector<Ptl*> v_bs, v_jpsi, v_muplus, v_muminus, v_kplus, v_kminus;
			 * vector<Vrt*> v_b_vrt, v_jpsi_kp_vrt, v_jpsi_km_vrt, v_bs_pv;
             * vector<double> v_lhtag, v_mb, v_emb;
             * vector<double> v_muplus_cpx, v_muplus_cpy, v_muplus_cpz;
             * vector<double> v_muminus_cpx, v_muminus_cpy, v_muminus_cpz;
             * vector<double> v_kplus_cpx, v_kplus_cpy, v_kplus_cpz;
             * vector<double> v_kminus_cpx, v_kminus_cpy, v_kminus_cpz;
		     */

                    //std::cout << "Bs: " << v_bs.size() 
                    //          << " JPsi: " << jpsi_finder->getIndex() 
                    //          << " Phi:" << phi_finder->getIndex() << std::endl;
		    v_jpsi_index.push_back(jpsi_finder->getIndex());
		    v_phi_index.push_back(phi_finder->getIndex());

		    v_bs.push_back(bs);
		    v_jpsi.push_back(jpsi);
		    v_muplus.push_back(muplus);
		    v_muminus.push_back(muminus);
		    v_kplus.push_back(kplus);
		    v_kminus.push_back(kminus);
		    v_bs_vrt.push_back(bs_vrt);
		    v_jpsi_kp_vrt.push_back(jpsi_kp_vrt);
		    v_jpsi_km_vrt.push_back(jpsi_km_vrt);
		    v_bs_pv.push_back(bs_pv);
		    v_lhtag.push_back(lhtag);
		    v_mb.push_back(mb);
		    v_emb.push_back(emb);
                    v_bs_ucpt.push_back(tmp_ucpt);

		    v_muplus_cpx.push_back(moms[0](1)); v_muplus_cpy.push_back(moms[0](2)); v_muplus_cpz.push_back(moms[0](3));
		    v_muminus_cpx.push_back(moms[1](1));v_muminus_cpy.push_back(moms[1](2));v_muminus_cpz.push_back(moms[1](3));
		    v_kplus_cpx.push_back(moms[2](1));  v_kplus_cpy.push_back(moms[2](2));  v_kplus_cpz.push_back(moms[2](3));
		    v_kminus_cpx.push_back(moms[3](1)); v_kminus_cpy.push_back(moms[3](2)); v_kminus_cpz.push_back(moms[3](3));
		}
  //std::cout << "Found: " << v_bs.size() << std::endl;
  return v_bs.size();
}

void BsJPsiPhiFinder::begin(){
	index = -1;
}

bool BsJPsiPhiFinder::next(){
	if (++index < v_bs.size())
		return true;
	else
		index=-1;
	return false;
}


void BsJPsiPhiFinder::fill(){
        //std::cout << "Filling: " << index << " JPsi:" << v_jpsi_index[index] << " Phi:" << v_phi_index[index] << std::endl;
	jpsi_finder->setIndex(v_jpsi_index[index]);
	phi_finder->setIndex(v_phi_index[index]);

	jpsi_finder->fill();     // Mu+, Mu-, Jpsi info.
	phi_finder->fill();      // K+, K-, Phi info.
	bs_saver.fill(*v_bs[index]);  // Bs info.
	vb_saver.fill(*v_bs_vrt[index]);  // B vertex info.
	vp_saver.fill(*v_bs_pv[index]);    // Primary vertex info.
if (mc_finder)
	mc_finder->fill(v_muplus[index], v_muminus[index], v_kplus[index], v_kminus[index]);

	PtlLst particle_list;
	particle_list.push_back(v_muplus[index]);
	particle_list.push_back(v_muminus[index]);
	particle_list.push_back(v_kplus[index]);
	particle_list.push_back(v_kminus[index]);
	tag_saver.fill(v_bs[index]->mom(), v_bs_pv[index], &particle_list);  // Tagging info.

	/* -- Info missing by the savers -- */
	bs_mass = v_mb[index];
	bs_mass_error = v_emb[index];

	/* -- Corrected momentums -- */
	mu_plus_cpx  = v_muplus_cpx[index];  mu_plus_cpy  = v_muplus_cpy[index];  mu_plus_cpz  = v_muplus_cpz[index];
	mu_minus_cpx = v_muminus_cpx[index]; mu_minus_cpy = v_muminus_cpy[index]; mu_minus_cpz = v_muminus_cpz[index];
	k_plus_cpx   = v_kplus_cpx[index];   k_plus_cpy   = v_kplus_cpy[index];   k_plus_cpz   = v_kplus_cpz[index];
	k_minus_cpx  = v_kminus_cpx[index];  k_minus_cpy  = v_kminus_cpy[index];  k_minus_cpz  = v_kminus_cpz[index];

	/* -- lifetime & lifetime_error -- */
    double ctau, vctau;
	v_bs[index]->decayLengthProper(PDG_BS_MASS, ctau, vctau, &particle_list);
	bs_pdl = ctau; bs_epdl=sqrt(fabs(vctau));

        /* -- guennadi likelihood ratio -- */
        bs_lhtag = v_lhtag[index];

	/* -- Calculate isolation of B -- */
	const AA::PtlLst* ptl_lst = AA::ptlBox.particles();
	double drmax = v_muminus[index]->dR(v_bs[index]->mom());
	if ( v_muplus[index]->dR(v_bs[index]->mom()) > drmax ) drmax = v_muplus[index]->dR(v_bs[index]->mom());
	if ( v_kplus[index]->dR(v_bs[index]->mom())  > drmax ) drmax = v_kplus[index]->dR(v_bs[index]->mom());
	if ( v_kminus[index]->dR(v_bs[index]->mom()) > drmax ) drmax = v_kminus[index]->dR(v_bs[index]->mom());
	double sum = 0; double sum_drmax = 0; double sum_75 = 0;
    double sum_pv = 0; double sum_drmax_pv = 0; double sum_75_pv = 0;
	for (PtlLstCIt piso = ptl_lst->begin(); piso != ptl_lst->end(); ++piso) {
		Ptl* ptlIso = *piso;
		if (ptlIso == v_muplus[index] || ptlIso == v_muminus[index] || ptlIso == v_kplus[index] || ptlIso == v_kminus[index])
			continue;
		double driso = ptlIso->dR(v_bs[index]->mom());
		if (driso <= 0.5)   sum += ptlIso->ptot();
		if (driso <= drmax) sum_drmax += ptlIso->ptot();
		if (driso <= 0.75)  sum_75 += ptlIso->ptot();
        if(v_jpsi[index]->primaryVertex() != ptlIso->primaryVertex()) continue;
        if (driso <= 0.5)   sum_pv += ptlIso->ptot();
        if (driso <= drmax) sum_drmax_pv += ptlIso->ptot();
        if (driso <= 0.75)  sum_75_pv += ptlIso->ptot();
	}
	bs_iso       = v_bs[index]->ptot() / (v_bs[index]->ptot() + sum);
	bs_iso_drmax = v_bs[index]->ptot() / (v_bs[index]->ptot() + sum_drmax);
	bs_iso_75    = v_bs[index]->ptot() / (v_bs[index]->ptot() + sum_75);
    bs_iso_pv       = v_bs[index]->ptot() / (v_bs[index]->ptot() + sum_pv);
    bs_iso_drmax_pv = v_bs[index]->ptot() / (v_bs[index]->ptot() + sum_drmax_pv);
    bs_iso_75_pv    = v_bs[index]->ptot() / (v_bs[index]->ptot() + sum_75_pv);

	/* -- Intermediate vertex chi2 -- */
	bs_jpsikp_chi2 = v_jpsi_kp_vrt[index]->chi2();
	bs_jpsikm_chi2 = v_jpsi_km_vrt[index]->chi2();

	/* -- Transversity Angles -- */
	TLorentzVector l_kplus, l_kminus, l_muplus, l_muminus;
	l_muplus.SetXYZM (v_muplus_cpx[index], v_muplus_cpy[index], v_muplus_cpz[index], PDG_MU_MASS);
	l_muminus.SetXYZM(v_muminus_cpx[index], v_muminus_cpy[index], v_muminus_cpz[index], PDG_MU_MASS);
	l_kplus.SetXYZM  (v_kplus_cpx[index], v_kplus_cpy[index], v_kplus_cpz[index], PDG_KAON_MASS);
	l_kminus.SetXYZM (v_kminus_cpx[index], v_kminus_cpy[index], v_kminus_cpz[index], PDG_KAON_MASS);

	threeAngles(l_muplus, l_muminus, l_kplus, l_kminus, bs_angle_phi, bs_angle_ctheta, bs_angle_cpsi);
	/* -- Phi corrected mass -- */
	vector<PType> phi_types(2);
	double mkkC,vmkkC;
	PtlLst phi_list;
	phi_list.push_back(v_kplus[index]);
	phi_list.push_back(v_kminus[index]);
	phi_types[0] = AA::K_PLUS;
	phi_types[1] = AA::K_MINUS;
	v_bs[index]->mass(phi_types,mkkC,vmkkC,&phi_list);
	phi_mass_corrected = mkkC;
	phi_mass_corrected_error = sqrt(fabs(vmkkC));

	/* -- corrected pt and bs uncorrected pt -- */
        bs_ucpt = v_bs_ucpt[index];
        mu_plus_cpt = sqrt(v_muplus_cpx[index]*v_muplus_cpx[index] + v_muplus_cpy[index]*v_muplus_cpy[index]);
        mu_minus_cpt = sqrt(v_muminus_cpx[index]*v_muminus_cpx[index] + v_muminus_cpy[index]*v_muminus_cpy[index]);
        k_plus_cpt = sqrt(v_kplus_cpx[index]*v_kplus_cpx[index] + v_kplus_cpy[index]*v_kplus_cpy[index]);
        k_minus_cpt = sqrt(v_kminus_cpx[index]*v_kminus_cpx[index] + v_kminus_cpy[index]*v_kminus_cpy[index]);
        jpsi_cpt = sqrt( (v_muplus_cpx[index]+v_muminus_cpx[index])*(v_muplus_cpx[index]+v_muminus_cpx[index])
                         + (v_muplus_cpy[index]+v_muminus_cpy[index])*(v_muplus_cpy[index]+v_muminus_cpy[index]) );
        phi_cpt  = sqrt( (v_kplus_cpx[index]+v_kminus_cpx[index])*(v_kplus_cpx[index]+v_kminus_cpx[index])  
                         + (v_kplus_cpy[index]+v_kminus_cpy[index])*(v_kplus_cpy[index]+v_kminus_cpy[index]) );
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
void BsJPsiPhiFinder::threeAngles(TLorentzVector mu_p,
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
