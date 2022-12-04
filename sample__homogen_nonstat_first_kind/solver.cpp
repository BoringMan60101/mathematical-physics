#include <cstdlib>
#include <cstdio>
#include "solver.hpp" //для проверки соответсвия прототипа функции (в заголовочнике) и её описания (здесь)
#include "my_types.hpp"
#include <cmath>

double ** getExactSolution(double *x, uns Np, double *cur_t, uns Nt) {
  double ** T_exact = new double *[Nt+1]; //Массив хранит значения аналитического решения во всех рассчётных точках
  CheckMemAlloc(T_exact, "Mem Error!\n In function \'T_exact\'\n");
  for(uns t = 0; t <= Nt; t++) {
    T_exact[t] = new double [Np];
    CheckMemAlloc(T_exact[t], "Mem Error!\n In function \'T_exact\'\n");
    T_exact[t][0] = 0.0f;
    T_exact[t][Np-1] = 0.0f;
    for(uns i = 1; i <= Np-2; i++)
      T_exact[t][i] = exp(-cur_t[t])*sin(x[i]); //!!!!!!!!!!!!!!!!!!!!!!!
  }
  return T_exact;
}
double * TDMA(double *T0, phaseSpace discrSpace, uns t, constPhysProperties Consts) {
  double *aW, *aP, *aE, *bP, *P, *Q, *T;
  uns Np =  discrSpace.Np;
  aW = new double [Np];
  aP = new double [Np];
  aE = new double [Np];
  bP = new double [Np];
  T = new double [Np];
  P = new double [Np];
  Q = new double [Np];
  if(!(aW && aP && aE && bP && P && Q && T)) {
    fprintf(stderr, "%s\n", "Mem Error!\n In function \'TDMA\'\n");
    exit(EXIT_FAILURE);
  }
  aW[0] = NAN; aP[0] = 1.0f; aE[0] = 0.0f;
  P[0] = 0.0f;
  T[0] = bP[0] = Q[0] = 0.0f; //Задано однородное граничное условие первого рода
  //Прямой ход
  for(uns i = 1; i <= Np-2; i++) {
    aW[i] = (discrSpace.dt[t]*0.5f*Consts.k) / (discrSpace.Nodes[i] - discrSpace.Nodes[i-1]);
    aE[i] = (discrSpace.dt[t]*0.5f*Consts.k) / (discrSpace.Nodes[i+1] - discrSpace.Nodes[i]);
    aP[i] = Consts.ro*Consts.c*(discrSpace.Faces[i+1] - discrSpace.Faces[i]) + aW[i] + aE[i];
    bP[i] = Consts.ro*Consts.c*(discrSpace.Faces[i+1] - discrSpace.Faces[i])*T0[i] + aE[i]*(T0[i+1] - T0[i]) - aW[i]*(T0[i] - T0[i-1]);
    P[i] = aE[i]/(aP[i] - aW[i]*P[i-1]);
    Q[i] = (aW[i]*Q[i-1] + bP[i])/(aP[i] - aW[i]*P[i-1]);
  }
  aW[Np-1] = 0.0f; aP[Np-1] = 1.0f; aE[Np-1] = NAN;
  P[Np-1] = 0.0f;
  T[Np-1] = bP[Np-1] = Q[Np-1] = 0.0f; //Задано однородное граничное условие первого рода
  //Обратный ход
  for(uns i = Np-2; i >= 1; i--)
    T[i] = P[i]*T[i+1] + Q[i];

  delete [] aW; delete [] aP; delete [] aE; delete [] P; delete [] Q;
  return T;
}
