/*
 * PtlSaver.cpp
 *
 *  Created on: Feb 21, 2009
 *      Author: magania
 */

#include "PtlSaver.h"
#include <stdio.h>
#include <string.h>
#include <Muo.hpp>

const char* PtlSaver::strccat(const char *a, const char *b){
	char *temp = new char[20];
	names.push_back(temp);
	strcpy(temp, a);
	const char *result = strcat(temp,b);
	return result;
}

PtlSaver::PtlSaver(const char *name, TTree &tree, bool _real, bool _muon) {
	real = _real;
	muon = _muon;

	tree.Branch(strccat(name, "_pt"), &pt, strccat(name, "_pt/D"));
	tree.Branch(strccat(name, "_ptot"), &ptot, strccat(name, "_ptot/D"));
	tree.Branch(strccat(name, "_phi"), &phi, strccat(name, "_phi/D"));
	tree.Branch(strccat(name, "_eta"), &eta, strccat(name, "_eta/D"));
	tree.Branch(strccat(name, "_ipA"), &ipA, strccat(name, "_ipA/D"));
	tree.Branch(strccat(name, "_vipA"), &vipA, strccat(name,"_vipA/D"));
	tree.Branch(strccat(name, "_ipS"), &ipS, strccat(name, "_ipS/D"));

	tree.Branch(strccat(name, "_vipS"), &vipS, strccat(name,"_vipS/D"));

	tree.Branch(strccat(name, "_px"), &px, strccat(name, "_px/D"));
	tree.Branch(strccat(name, "_py"), &py, strccat(name, "_py/D"));
	tree.Branch(strccat(name, "_pz"), &pz, strccat(name, "_pz/D"));

	tree.Branch(strccat(name, "_vp11"), &vp11, strccat(name, "_vp11/D"));
	tree.Branch(strccat(name, "_vp12"), &vp12, strccat(name, "_vp12/D"));
	tree.Branch(strccat(name, "_vp13"), &vp13, strccat(name, "_vp13/D"));
	tree.Branch(strccat(name, "_vp22"), &vp22, strccat(name, "_vp22/D"));
	tree.Branch(strccat(name, "_vp23"), &vp23, strccat(name, "_vp23/D"));
	tree.Branch(strccat(name, "_vp33"), &vp33, strccat(name, "_vp33/D"));

	if (muon){
		tree.Branch(strccat(name, "_nseg"), &nseg, strccat(name, "_nseg/I"));
		tree.Branch(strccat(name,"_calesig"), &calesig, strccat(name,"_calesig/D"));
	}

	if (real) {
		tree.Branch(strccat(name, "_chi2"), &chi2, strccat(name, "_chi2/D"));
		tree.Branch(strccat(name, "_q"), &q, strccat(name, "_q/D"));
		tree.Branch(strccat(name, "_dedx"), &dedx, strccat(name, "_dedx/D"));
		tree.Branch(strccat(name, "_nSMT"), &nSMT, strccat(name, "_nSMT/I"));
		tree.Branch(strccat(name, "_nCFT"), &nCFT, strccat(name, "_nCFT/I"));
	}
}

PtlSaver::~PtlSaver() {
	// TODO Auto-generated destructor.
}

void PtlSaver::fill(AA::Ptl &p) {
	pt = p.pt();
	ptot = p.ptot();
	phi = p.phi();
	eta = p.eta();
	ipA = p.imp(1);
	vipA = p.vimp(1, 1);
	ipS = p.imp(2);
	vipS = p.vimp(2, 2);

	px = p.mom(1);
	py = p.mom(2);
	pz = p.mom(3);

	vp11 = p.vmom(1,1);
	vp12 = p.vmom(1,2);
	vp13 = p.vmom(1,3);
	vp22 = p.vmom(2,2);
	vp23 = p.vmom(2,3);
	vp33 = p.vmom(3,3);

	if (muon){
		AA::MuoGlb *m = p.muon();
		if (m) {
			nseg = m->nseg();
			calesig = m->calEsig();
		} else {
			nseg = -1;
			calesig = -1;
		}
	}

	if (real) {
		chi2 = p.track()->chi2();
		q = p.q();
		dedx = p.track()->dedx();
		nSMT = p.nSMT();
		nCFT = p.nCFT();
	}
}
