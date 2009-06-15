#include <DecayMC.h>

/* Regresa verdadero si p decae en id_c1 y id_c2
 *
 * Busca si los hijos de p tienen PDG ID id_c1 y id_c2, de ser asi
 * guarda los hijos en c1 y c2 respectivamente.
 */
bool DecayMC::asignChildren(PtlMC* p, int id_c1, PtlMC** c1, int id_c2,
		PtlMC** c2) {
	if (!p->vrtEnd())
		return false;
	PtlMCLst child_list = p->vrtEnd()->children();
	//    cout << "count " << child_list.size()  << endl;
	if (child_list.size() != 2)
		return false;
	*c1 = 0;
	*c2 = 0;
	for (PtlMCLstCIt c = child_list.begin(); c != child_list.end(); ++c) {
		//cout << "id " << (*c)->idPDG()  << endl;
		if (*c1 == 0 && (*c)->idPDG() == id_c1) {
			*c1 = *c;
			continue;
		}
		if (*c2 == 0 && (*c)->idPDG() == id_c2) {
			*c2 = *c;
			continue;
		}
		return false;
	}
	//    if (( *c1 && *c2 ) && *c1==0)
	//      cout << (*c1) << ' ' << (*c2)  << ' ' << (*c1 && *c2) << endl;
	return *c1 && *c2;
}

/* Busca decaimientos primary -> c1(c1c1 c1c2) c2(c2c1 c2c1)
 *
 * Recibe el PDG ID de la particulas buscadas y asigna los PtlMC. Regresa el
 * numero de decaimientos del tipo encontrados aunque solo asigna el ultimo
 */
int DecayMC::findAndAsignDecay(int id_primary, PtlMC** primary,
		int id_c1, PtlMC** c1, int id_c2, PtlMC** c2, int id_c1c1,
		PtlMC** c1c1, int id_c1c2, PtlMC** c1c2, int id_c2c1, PtlMC** c2c1,
		int id_c2c2, PtlMC** c2c2) {
	PtlMC* tc1, *tc2, *tc1c1, *tc1c2, *tc2c1, *tc2c2;
	int found = 0;
	for (PtlMCLstCIt p = ptlMCBox.particles()->begin(); p
			!= ptlMCBox.particles()->end(); ++p) {
		if ((*p)->idPDG() == id_primary) {
			*primary = *p;

			if (!asignChildren(*primary, id_c1, &tc1, id_c2, &tc2))
				continue;
			if (!asignChildren(tc1, id_c1c1, &tc1c1, id_c1c2, &tc1c2))
				continue;
			if (!asignChildren(tc2, id_c2c1, &tc2c1, id_c2c2, &tc2c2))
				continue;
			*c1 = tc1;
			*c2 = tc2;
			*c1c1 = tc1c1;
			*c1c2 = tc1c2;
			*c2c1 = tc2c1;
			*c2c2 = tc2c2;
			found++;
		}
	}

	//    if ( *c2c1==0 && found >0)
	//      cout << "Que Pasa: "<< *c2 << ' ' << id_c2c1 << endl;
	return found;
}

/* Ahora con decaimientos primary -> c1(c1c1 c1c2) c2
 */
int DecayMC::findAndAsignDecay(int id_primary, PtlMC** primary,
		int id_c1, PtlMC** c1, int id_c2, PtlMC** c2, int id_c1c1,
		PtlMC** c1c1, int id_c1c2, PtlMC** c1c2) {
	PtlMC* tc1, *tc2, *tc1c1, *tc1c2;
	int found = 0;
	for (PtlMCLstCIt p = ptlMCBox.particles()->begin(); p
			!= ptlMCBox.particles()->end(); ++p) {
		if ((*p)->idPDG() == id_primary) {
			*primary = *p;

			if (!asignChildren(*primary, id_c1, &tc1, id_c2, &tc2))
				continue;
			if (!asignChildren(tc1, id_c1c1, &tc1c1, id_c1c2, &tc1c2))
				continue;
			*c1 = tc1;
			*c2 = tc2;
			*c1c1 = tc1c1;
			*c1c2 = tc1c2;
			found++;
		}
	}

	return found;
}

/* Busca decaimientos primary -> left right
 *
 * Recibe el PDG ID de la particulas buscadas y asigna los PtlMC. Regresa el
 * numero de decaimientos del tipo encontrados aunque solo asigna el ultimo
 */
int DecayMC::findAndAsignDecay(int id_primary, vector<PtlMC*> *primary,
		int id_left, vector<PtlMC*> *left, int id_right, vector<PtlMC*> *right) {
	primary->clear();
	left->clear();
	right->clear();
	int found = 0;
	for (PtlMCLstCIt p = ptlMCBox.particles()->begin(); p
			!= ptlMCBox.particles()->end(); ++p) {
		//      cout << (*p)->idPDG() << endl;
		if ((*p)->idPDG() == id_primary) {
			PtlMC *left_mc, *right_mc;
			if (!asignChildren(*p, id_left, &left_mc, id_right, &right_mc))
				continue;
			primary->push_back(*p);
			left->push_back(left_mc);
			right->push_back(right_mc);
			found++;
		}
	}
	return found;
}

/* Hace el matching entre pmc y lo guarda en p
 *
 * Calcula chi2 entre pmc y todas las particulas del evento, la que regrese
 * el menor chi2 (siempre que sea mayor a chi) la guarda en p.
 */
bool DecayMC::match(PtlMC* pmc, Ptl** p, double chi) {
	PtlLst lst_P = *AA::ptlBox.particles();
	TrkLst lst;
	for (PtlLstCIt p1 = lst_P.begin(); p1 != lst_P.end(); ++p1)
		lst.push_back(*p1);

	Trk* p_t;

	double chi_t = 1000;
	pmc->associate(lst, p_t, chi_t);

	if (chi_t > chi)
		return false;
	*p = (Ptl*) p_t;
	return true;
}

bool DecayMC::match(PtlMC* pmc_1, Ptl** p_1, PtlMC* pmc_2, Ptl** p_2,
		double chi) {
	PtlLst lst_P = *AA::ptlBox.particles();
	TrkLst lst;
	for (PtlLstCIt p1 = lst_P.begin(); p1 != lst_P.end(); ++p1)
		lst.push_back(*p1);

	Trk* p1_t;
	Trk* p2_t;

	double chi_t1 = 1000;
	double chi_t2 = 1000;
	pmc_1->associate(lst, p1_t, chi_t1);
	pmc_2->associate(lst, p2_t, chi_t2);

	if (chi_t1 > chi || chi_t2 > chi)
		return false;
	*p_1 = (Ptl*) p1_t;
	*p_2 = (Ptl*) p2_t;

	return true;
}

/* Regresa si la particula es un Meson B */
bool DecayMC::isB(PtlMC* p) {
	if (!p)
		return false;
	int id = abs(p->idPDG());
	if (id > 500 && id < 600)
		return true;
	if (id > 10500 && id < 10600)
		return true;
	if (id > 20500 && id < 20600)
		return true;
	return false;
}

/* Rergesa el primer padre del vertice de creacion de la particula */
PtlMC* DecayMC::getParent(PtlMC* p) {
	if (!p)
		return 0;
	if (!p->vrtStart())
		return 0;
	if (p->vrtStart()->parents().size() > 0)
		return p->vrtStart()->parents().front();
	return 0;
}

/* Busca la particula en el MC a nivel generacion mas parecida a p*/
PtlMC* DecayMC::imatch(Ptl* p) {
	TrkLst lst;
	lst.push_back(p->track());

	Trk* p_t;
	double chi = 100;
	double chi_min = 100;
	PtlMC* m = 0;
	for (PtlMCLstCIt pmc = ptlMCBox.particles()->begin(); pmc
			!= ptlMCBox.particles()->end(); ++pmc) {
		(*pmc)->associate(lst, p_t, chi);
		if (chi < chi_min) {
			chi = chi_min;
			m = (*pmc);
		}
	}
	return m;
}

int DecayMC::getIdPdg(PtlMC* p) {
	if (!p)
		return 0;
	return p->idPDG();
}

int DecayMC::getIndex(PtlMC* p) {
	if (!p)
		return 0;
	return p->index();
}

/* Hace lo mismo que match pero para 4 particulas
 */
bool DecayMC::match(PtlMC* p1_mc, Ptl**pt1, PtlMC* p2_mc, Ptl** pt2,
		PtlMC* p3_mc, Ptl**pt3, PtlMC* p4_mc, Ptl** pt4, int chi) {
	PtlLst lst_P = *AA::ptlBox.particles();
	TrkLst lst;
	for (PtlLstCIt p1 = lst_P.begin(); p1 != lst_P.end(); ++p1)
		lst.push_back(*p1);

	Trk* p1_t, *p2_t, *p3_t, *p4_t;

	double chi_p1 = 1000;
	double chi_p2 = 1000;
	double chi_p3 = 1000;
	double chi_p4 = 1000;
	p1_mc->associate(lst, p1_t, chi_p1);
	p2_mc->associate(lst, p2_t, chi_p2);
	p3_mc->associate(lst, p3_t, chi_p3);
	p4_mc->associate(lst, p4_t, chi_p4);

	if (chi_p1 > chi || chi_p2 > chi || chi_p3 > chi || chi_p4 > chi)
		return false;

	*pt1 = (Ptl*) p1_t;
	*pt2 = (Ptl*) p2_t;
	*pt3 = (Ptl*) p3_t;
	*pt4 = (Ptl*) p4_t;

	return true;
}

/* Hace lo mismo que match pero para 3 particulas
 */
bool DecayMC::match(PtlMC* p1_mc, Ptl**pt1, PtlMC* p2_mc, Ptl** pt2,
		PtlMC* p3_mc, Ptl**pt3, double chi) {
	PtlLst lst_P = *AA::ptlBox.particles();
	TrkLst lst;
	for (PtlLstCIt p1 = lst_P.begin(); p1 != lst_P.end(); ++p1)
		lst.push_back(*p1);

	Trk* p1_t, *p2_t, *p3_t;

	double chi_p1 = 1000;
	double chi_p2 = 1000;
	double chi_p3 = 1000;
	p1_mc->associate(lst, p1_t, chi_p1);
	p2_mc->associate(lst, p2_t, chi_p2);
	p3_mc->associate(lst, p3_t, chi_p3);

	if (chi_p1 > chi || chi_p2 > chi || chi_p3 > chi)
		return false;

	*pt1 = (Ptl*) p1_t;
	*pt2 = (Ptl*) p2_t;
	*pt3 = (Ptl*) p3_t;

	return true;
}
