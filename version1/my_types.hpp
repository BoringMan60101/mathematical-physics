/*

*/
#ifndef __MY_TYPES__
#define __MY_TYPES__
typedef unsigned uns;

struct physConsts {
  double ro, c, k;
};

//Дискретный аналог фазового пространства
//В нём также хранятся параметры разбиения
struct discrMesh {
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

struct errors {
  uns Nt, Np;
  double *abs; //Абсолютные значения погрешности для текущего временного слоя
  double maxAbs; //
  double *rel; //
  double maxRel; //
};

#endif
