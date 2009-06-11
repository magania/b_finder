/*
 * EvtSaver.cpp
 *
 *  Created on: Mar 16, 2009
 *      Author: magania
 */

#include <EvtSaver.h>

EvtSaver::EvtSaver(TTree &tree){
	tree.Branch("run", &run, "run/I");
	tree.Branch("event", &event, "event/I");
	tree.Branch("n_triggers", &n_triggers, "n_triggers/I");
	tree.Branch("len_triggers", &len_triggers, "len_triggers/I");
	tree.Branch("triggers", &triggers, "triggers[len_triggers]/C");
}

EvtSaver::~EvtSaver(){
	return;
}

void EvtSaver::fill(){
		run = AA::runNumber;
		event = AA::evtNumber;

		std::string AllTrg = "";
		std::string tmptrig;
		bool debug = false;
		bool showall = false;

		n_triggers = 0;
		if (debug)
			std::cout << AA::runNumber << " " << AA::evtNumber << " ";
		if (AA::trigger.triggers()->size() == 0) {
			if (debug && showall)
				std::cout << "***** No trigger found, recheking..." << std::endl;
			len_triggers = sprintf(triggers, "%s", AllTrg.c_str());;
			return;
		}

		// L3 Trigger names for this event
		for (std::list<std::string>::const_iterator p = AA::trigger.triggers()->begin(); p
				!= AA::trigger.triggers()->end(); ++p) {
			tmptrig = *p;
			tmptrig += " ";
			//if((AllTrg.size() + tmptrig.size())<10000)
			AllTrg += tmptrig;
			if (debug)
				std::cout << " " << *p;
			n_triggers++;
		}

		if (debug)
			std::cout << std::endl;

		len_triggers = sprintf(triggers, "%s", AllTrg.c_str());

		if (debug && showall)
			std::cout << " Triggers :  " << AllTrg << std::endl;
		if (debug && showall)
			std::cout << " BANA INFO: Triggers :  " << triggers << std::endl;
}
