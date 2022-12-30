/* В этом модуле реализованы функции для определения погрешностей */
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "my_types.hpp"
void checkMemAlloc(void *pointer, const char *ErrMsg);

//Определение абсолютной погрешности
//(применяется для каждого временного слоя)
double *getAbsErr(double *curT_exact, double *curT, uns Np) {
  double *curAbsErr = new double [Np];
  checkMemAlloc(curAbsErr, "Mem Error!\n In func \'getAbsErr\'\n");
  for(uns i = 0; i < Np; i++)
    curAbsErr[i] = fabs(curT_exact[i] - curT[i]);
  return curAbsErr;
}

//Определение максимальной абсолютной погрешности.
//Получается по одному значению с каждого временного слоя.
double getMaxAbsErr(double *curAbsErr, uns Np) {
  double maxAbsErr = curAbsErr[0];
  for(uns i = 1; i < Np; i++)
    if(curAbsErr[i] > maxAbsErr)
      maxAbsErr = curAbsErr[i];
  return maxAbsErr;
}

//Определение относительной погрешности.
//Применяется для каждого временного слоя.
double *getRelErr(double *curT_exact,
                  double *curT, const discrMesh &mesh, uns time) {
  //Получение усреднённого значения точного решения
  double Texact_Avg=-cos(mesh.Nodes[mesh.Np-1])+cos(mesh.Nodes[0]);
  Texact_Avg = (1.0f - exp(-mesh.cur_t[time])) / mesh.len;

  double *curRelErr = new double [mesh.Np];
  checkMemAlloc(curRelErr, "Mem Error!\n In func \'getRelErr\'\n");
  for(uns i = 0; i < mesh.Np; i++)
    curRelErr[i] = 100.0f*(fabs(curT_exact[i]-curT[i])/Texact_Avg);
  return curRelErr;
}

//Определение максимальной относительной погрешности.
//Получается по одному значению с каждого временного слоя.
double getMaxRelErr(double *curRelErr, uns Np) {
  double maxRel = curRelErr[0];
  for(uns i = 1; i < Np; i++)
    if(curRelErr[i] > maxRel)
      maxRel = curRelErr[i];
  return maxRel;
}
