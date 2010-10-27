/*
 * D0Finder.h
 *
 *  Created on: Oct 22, 2010
 *      Author: magania
 */

#ifndef D0FINDER_H_
#define D0FINDER_H_

#include <PtlSaver.h>
#include <VrtSaver.h>

#include <vector>
#include <list>

#include <AA.hpp>
#include <Ptl.hpp>
#include <V0Finder.hpp>
#include <select.hpp>

#include <TTree.h>

class D0Finder {
public:
	D0Finder(TTree &tree);
	virtual ~D0Finder();

	int find();
	void begin();
	bool next();
	int getIndex();
	void setIndex(int i);

	AA::Ptl& getKaon();
	AA::Ptl& getPion();
	AA::Ptl& getD0();
	AA::Vrt& getVrt();
	double getMass();
	void fill();

private:
	AA::PtlBox _boxp;
	AA::VrtBox _boxv;

	void clean();

	double masa;

	PtlSaver kaon_saver, pion_saver, d0_saver;
	VrtSaver vrt_saver;
	std::vector<AA::Ptl*> d0, kaon, pion;
	std::vector<AA::Vrt*> vtx;
	std::vector<double> mass;

	int index;

	const static double MASS_D0_MIN = 1.4;
	const static double MASS_D0_MAX = 2.2;

	const static double KAON_MINIMUM_PT = 0.7;
	const static double PION_MINIMUM_PT = 0.7;
	double MINIMUM_PT;

	const static double MINIMUM_D0_PT = 5.0;

	const static double VRT_CHI2_MIN = 9.0;

	const static int MINIMUM_NSMT = 2;
	const static int MINIMUM_NCFT = 2;

};

#endif /* D0FINDER_H_ */
