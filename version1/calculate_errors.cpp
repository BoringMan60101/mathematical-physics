#include <cmath>
#include "my_types.hpp"
void checkMemAlloc(void *pointer, const char *ErrMsg);
double * getAbsErr(double *curT_exact, double *curT, uns Np) {
  double *curAbsErr = new double [Np];
  checkMemAlloc(curAbsErr, "Mem Error!\n In function \'getAbsErr\'\n");
  for(uns i = 0; i < Np; i++)
    curAbsErr[i] = fabs(curT_exact[i] - curT[i]);
  return curAbsErr;
}
double getMaxAbsErr(double *curAbsErr, uns Np) {
  double maxAbsErr = curAbsErr[0];
  for(uns i = 1; i < Np; i++)
    if(curAbsErr[i] > maxAbsErr)
      maxAbsErr = curAbsErr[i];
  return maxAbsErr;
}
double T_exactFormula(const double time_layer, const double x) {
  return (1.0f - exp(-time_layer))*sin(x);
}
double getAverageOfCurT_exact(double time_layer, double len) {
  //Считаем среднее значение T по выбранному временное слою t
  uns n = 2;
  double I1, I2, dx;
  double integralError = 1e-12;
  do {
    I1 = I2 = 0.0f;
    dx = len / n;
    for(uns i = 0; i <= n; i++)
      I1 += T_exactFormula(time_layer, i*dx);
    I1 = I1 * dx;
    n *= 2;
    dx = len / n;
    for(uns i = 0; i <= n; i++)
      I2 += T_exactFormula(time_layer, i*dx);
    I2 = I2 * dx;
  } while(fabs(I2 - I1) >= integralError);
  return I2/len;
}
double * getRelErr(double *curT_exact, double *curT, double time_layer, uns Np, double len) {
  //Получаем усреднённое значение точного решения для данного временного слоям
  double T_exactAvg = getAverageOfCurT_exact(time_layer, len);

  double *curRelErr = new double [Np];
  checkMemAlloc(curRelErr, "Mem Error!\n In function \'getRelErr\'\n");
  for(uns i = 0; i < Np; i++)
    curRelErr[i] = fabs((curT_exact[i] - curT[i])/T_exactAvg);
  return  curRelErr;
}
double getMaxRelErr(double *curRelErr, uns Np) {
  double maxRel = curRelErr[0];
  for(uns i = 1; i < Np; i++)
    if(curRelErr[i] > maxRel)
      maxRel = curRelErr[i];
  return maxRel;
}
