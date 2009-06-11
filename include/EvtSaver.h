/*
 * EvtSaver.h
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#ifndef EVTSAVER_H_
#define EVTSAVER_H_

#include <vector>
#include <string>
#include <list>

#include <AA.hpp>
#include <Trg.hpp>

#include <TTree.h>

class EvtSaver {
public:
	EvtSaver(TTree &tree);
	~EvtSaver();

	void fill();

private:
//	const char* strccat(const char *a, const char *b);

//	std::vector<char*> names;

	int run, event;
	int n_triggers;
	int len_triggers;
	char triggers[10000];
};

#endif /* PTLSAVER_H_ */
