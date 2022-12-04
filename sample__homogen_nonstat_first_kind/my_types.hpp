/*

*/
#ifndef __MY_TYPES__
#define __MY_TYPES__
typedef unsigned uns;

struct constPhysProperties {
  double ro, c, k;
  //
};

//Дискретный аналог фазового пространства
//В нём также хранятся параметры разбиения:
// (шаг между гранями, интервал между узловыми точками, шаг по времени)
struct phaseSpace {
  double len; //
  uns Nv, Np, Nt; //
  double *Faces; //
  double *Nodes; //
  double *FacesIntervals; //
  double *NodesIntervals; //
  double T_end; //
  double *cur_t; //
  double *dt; //
};

struct errsArrays {
  uns Nt, Np;
  double *abs; //Абсолютные значения погрешности для текущего временного слоя
  double maxAbs; //
  double *rel; //
  double maxRel; //
};

#endif
