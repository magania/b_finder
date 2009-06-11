/*
 * PtlSaver.h
 *
 *  Created on: Feb 21, 2009
 *      Author: magania
 */

#ifndef PTLSAVER_H_
#define PTLSAVER_H_

#include <vector>

#include <Ptl.hpp>

#include <TTree.h>

class PtlSaver {
public:
	PtlSaver(const char *name, TTree &tree, bool _real = true, bool _muon = false);
	virtual ~PtlSaver();

	void fill(AA::Ptl &p);

private:
	const char* strccat(const char *a, const char *b);

	std::vector<char*> names;

	bool real, muon;

	double pt, ptot, phi, eta;
	int nseg, nSMT, nCFT;
	double chi2, q, dedx;
	double ipA, vipA, ipS, vipS;
	double calesig;
	double px,py,pz;
	double vp11, vp12, vp13, vp22, vp23, vp33;
};

#endif /* PTLSAVER_H_ */
