/*
 * SSMuFinder.h
 *
 *  Created on: Nov 12, 2011
 *      Author: magania
 */

#ifndef SSMUFINDER_H_
#define SSMUFINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <AA.hpp>
#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class SSMuFinder {
public:
	SSMuFinder(TTree &tree);
	virtual ~SSMuFinder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getMu1();
	AA::Ptl& getMu2();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	PtlSaver mu1_saver, mu2_saver;
	VrtSaver vrt_saver;
	std::vector<AA::Ptl*> v_mu1, v_mu2;
	std::vector<AA::Vrt*> v_vtx;

	int index;
};

#endif /* SSMUFINDER_H_ */
