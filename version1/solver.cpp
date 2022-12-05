#include <cstdlib>
#include <cstdio>
#include "solver.hpp" //для проверки соответсвия прототипа функции (в заголовочнике) и её описания (здесь)
#include "my_types.hpp"
#include <cmath>

double ** getExactSolution(const discrMesh &mesh) {
  double ** T_exact = new double *[mesh.Nt+1]; //Массив хранит значения аналитического решения во всех рассчётных точках
  checkMemAlloc(T_exact, "Mem Error!\n In function \'T_exact\'\n");
  for(uns t = 0; t <= mesh.Nt; t++) {
    T_exact[t] = new double [mesh.Np];
    checkMemAlloc(T_exact[t], "Mem Error!\n In function \'T_exact\'\n");
    T_exact[t][0] = 0.0f;
    T_exact[t][mesh.Np-1] = 0.0f;
    for(uns i = 1; i <= mesh.Np-2; i++)
      T_exact[t][i] = (1.0f - exp(-mesh.cur_t[t]))*sin(mesh.Nodes[i]);
  }
  return T_exact;
}

double * TDMA(double *T0, const discrMesh &mesh, const physConsts &consts, uns t) {
  double aW, aP, aE, bP, *P, *Q, *T;
  double S; //Вклад стационарного источника
  uns Np =  mesh.Np;
  //aW = new double [Np];
  //aP = new double [Np];
  //aE = new double [Np];
  //bP = new double [Np];
  T = new double [Np];
  P = new double [Np];
  Q = new double [Np];
  if(!(P && Q && T)) {
    fprintf(stderr, "%s\n", "Mem Error!\n In function \'TDMA\'\n");
    exit(EXIT_FAILURE);
  }

  //aW[0] = NAN; aP[0] = 1.0f; aE[0] = 0.0f;
  P[0] = 0.0f;
  T[0] = Q[0] = 0.0f; //Задано однородное граничное условие первого рода
  //T[0] = bP[0] = Q[0] = 0.0f; //Задано однородное граничное условие первого рода
  //Прямой ход
  for(uns i = 1; i <= Np-2; i++) {
    aW = (mesh.dt[t]*0.5f*consts.k) / (mesh.Nodes[i] - mesh.Nodes[i-1]);
    aE = (mesh.dt[t]*0.5f*consts.k) / (mesh.Nodes[i+1] - mesh.Nodes[i]);
    aP = consts.ro*consts.c*(mesh.Faces[i+1] - mesh.Faces[i]) + aW + aE;
    //Выражение для bP[i] сделать более разборчивым
    S = sin(mesh.Nodes[i])*consts.ro*consts.c*(mesh.Faces[i+1] - mesh.Faces[i])*mesh.dt[i];
    bP = consts.ro*consts.c*(mesh.Faces[i+1] - mesh.Faces[i])*T0[i];
    bP += aE*(T0[i+1] - T0[i]) - aW*(T0[i] - T0[i-1]) + S;
    P[i] = aE/(aP - aW*P[i-1]);
    Q[i] = (aW*Q[i-1] + bP) / (aP - aW*P[i-1]);
  }
  //aW[Np-1] = 0.0f; aP[Np-1] = 1.0f; //aE[Np-1] = NAN;
  P[Np-1] = 0.0f;
  T[Np-1] = Q[Np-1] = 0.0f; //Задано однородное граничное условие первого рода
  //T[Np-1] = bP[Np-1] = Q[Np-1] = 0.0f; //Задано однородное граничное условие первого рода
  //Обратный ход
  for(uns i = Np-2; i >= 1; i--)
    T[i] = P[i]*T[i+1] + Q[i];

  delete [] P; delete [] Q;
  return T;
}
