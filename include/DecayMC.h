#ifndef _DECAYMC_H
#define _DECAYMC_H

#include <stdlib.h>
#include <AA.hpp>
#include <PtlMC.hpp>
#include <Ptl.hpp>

using namespace std;
using namespace AA;

class DecayMC {
public:
  /* Regresa verdadero si p decae en id_c1 y id_c2
   *
   * Busca si los hijos de p tienen PDG ID id_c1 y id_c2, de ser asi
   * guarda los hijos en c1 y c2 respectivamente.
   */
  static bool asignChildren(PtlMC* p, int id_c1, PtlMC** c1, int id_c2, PtlMC** c2);

  /* Busca decaimientos primary -> c1(c1c1 c1c2) c2(c2c1 c2c1)
   *
   * Recibe el PDG ID de la particulas buscadas y asigna los PtlMC. Regresa el
   * numero de decaimientos del tipo encontrados aunque solo asigna el ultimo
   */
  static int findAndAsignDecay(int id_primary, PtlMC** primary,
			       int id_c1,   PtlMC** c1,   int id_c2,   PtlMC** c2,
			       int id_c1c1, PtlMC** c1c1, int id_c1c2, PtlMC** c1c2,
			       int id_c2c1, PtlMC** c2c1, int id_c2c2, PtlMC** c2c2);

  /* Ahora con decaimientos primary -> c1(c1c1 c1c2) c2
   */
  static int findAndAsignDecay(int id_primary, PtlMC** primary,
			       int id_c1,   PtlMC** c1,   int id_c2,   PtlMC** c2,
			       int id_c1c1, PtlMC** c1c1, int id_c1c2, PtlMC** c1c2);

  /* Busca decaimientos primary -> left right
   *
   * Recibe el PDG ID de la particulas buscadas y asigna los PtlMC. Regresa el
   * numero de decaimientos del tipo encontrados aunque solo asigna el ultimo
   */
  static int findAndAsignDecay(int id_primary, vector<PtlMC*> *primary,
			       int id_left, vector<PtlMC*> *left,
			       int id_right, vector<PtlMC*> *right);

  /* Hace el matching entre pmc y lo guarda en p
   *
   * Calcula chi2 entre pmc y todas las particulas del evento, la que regrese
   * el menor chi2 (siempre que sea mayor a chi) la guarda en p.
   */
  static  bool match(PtlMC* pmc, Ptl** p, double chi = 99.0);

  static  bool match(PtlMC* pmc_1, Ptl** p_1, PtlMC* pmc_2, Ptl** p_2,  double chi = 99.0);

  /* Regresa si la particula es un Meson B */
  static bool isB(PtlMC* p);

  /* Rergesa el primer padre del vertice de creacion de la particula */
  static PtlMC* getParent(PtlMC* p);

  /* Busca la particula en el MC a nivel generacion mas parecida a p*/
  static PtlMC* imatch(Ptl* p);

  static int getIdPdg(PtlMC* p);
  static int getIndex(PtlMC* p);

  /* Hace lo mismo que match pero para 4 particulas
   */
  static  bool match(PtlMC* p1_mc, Ptl**pt1, PtlMC* p2_mc, Ptl** pt2,
	     PtlMC* p3_mc, Ptl**pt3, PtlMC* p4_mc, Ptl** pt4, int chi = 99.0);

  /* Hace lo mismo que match pero para 3 particulas
   */
  static  bool match(PtlMC* p1_mc, Ptl**pt1, PtlMC* p2_mc, Ptl** pt2,
	     PtlMC* p3_mc, Ptl**pt3, double chi = 99.0);

};

#endif
