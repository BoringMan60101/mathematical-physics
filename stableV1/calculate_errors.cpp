/* В этом модуле реализованы функции для определения погрешностей */
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "my_types.hpp"
void checkMemAlloc(void *pointer, const char *ErrMsg);

//Определение абсолютной погрешности (применяется для каждого временного слоя)
double * getAbsErr(double *curT_exact, double *curT, uns Np) {
  double *curAbsErr = new double [Np];
  checkMemAlloc(curAbsErr, "Mem Error!\n In function \'getAbsErr\'\n");
  for(uns i = 0; i < Np; i++)
    curAbsErr[i] = fabs(curT_exact[i] - curT[i]);
  return curAbsErr;
}

//Определение максимальной абсолютной погрешности (одно значение с каждого слоя)
double getMaxAbsErr(double *curAbsErr, uns Np) {
  double maxAbsErr = curAbsErr[0];
  for(uns i = 1; i < Np; i++)
    if(curAbsErr[i] > maxAbsErr)
      maxAbsErr = curAbsErr[i];
  return maxAbsErr;
}

//Определение относительной погрешности (применяется для каждого временного слоя)
double * getRelErr(double *curT_exact, double *curT, const discrMesh &mesh, uns time) {
  //Интегрирование точного решения на заданном временном слое
  double integral = -cos(mesh.Nodes[mesh.Np-1]) + cos(mesh.Nodes[0]);

  //Получение усреднённого значения точного решения
  double Texact_Avg = fabs(1.0f - exp(-mesh.cur_t[time])) * (integral/mesh.len);
  
  double *curRelErr = new double [mesh.Np];
  checkMemAlloc(curRelErr, "Mem Error!\n In function \'getRelErr\'\n");
  for(uns i = 0; i < mesh.Np; i++)
    curRelErr[i] = 100.0f*fabs(curT_exact[i] - curT[i]) / Texact_Avg;
  return curRelErr;
}

//Определение максимальной относительной погрешности (одно значение с каждого слоя)
double getMaxRelErr(double *curRelErr, uns Np) {
  double maxRel = curRelErr[0];
  for(uns i = 1; i < Np; i++)
    if(curRelErr[i] > maxRel)
      maxRel = curRelErr[i];
  return 100.f*maxRel;
}
