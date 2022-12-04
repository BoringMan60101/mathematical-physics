#include <cmath>
#include "my_types.hpp"
void CheckMemAlloc(void *pointer, const char *ErrMsg);
double * getAbsErr(double **T_exact, double *T, uns Np, uns t) {
  double *TabsErr = new double [Np];
  CheckMemAlloc(TabsErr, "Mem Error!\n In function \'getAbsErr\'\n");
  for(uns i = 0; i < Np; i++)
    TabsErr[i] = fabs(T_exact[t][i] - T[i]);
  return TabsErr;
}
double getMaxAbsErr(double *absErr, uns Np) {
  double maxAbsErr = absErr[0];
  for(uns i = 1; i < Np; i++)
    if(absErr[i] > maxAbsErr)
      maxAbsErr = absErr[i];
  return maxAbsErr;
}
double getT_exactAvg(double time_layer, double len) {
  //Считаем среднее значение T по выбранному временное слою t
  uns n = 2;
  double I1 = 0.0f, I2 = 0.0f, dx = len / n, eps = 1e-12;
  do {
    I1 = I2 = 0.0f;
    dx = len / n;
    for(uns i = 0; i <= n; i++)
      I1 += sin(i*dx);
    I1 = I1 * dx;
    n *= 2;
    dx = len / n;
    for(uns i = 0; i <= n; i++)
      I2 += sin(i*dx);
    I2 = I2 * dx;
  } while(fabs(I2 - I1) > eps);
  return (exp(-time_layer)*I2)/len;
}
double * getRelErr(double **T_exact, double *T, uns Np, double time_layer, uns t, double len) {
  //Получаем усреднённое значение точного решения для данного временного слоям
  double T_exactAvg = getT_exactAvg(time_layer, len);
  double *TrelErr = new double [Np];
  CheckMemAlloc(TrelErr, "Mem Error!\n In function \'getRelErr\'\n");
  for(uns i = 0; i < Np; i++)
    TrelErr[i] = fabs((T_exact[t][i] - T[i])/T_exactAvg);
  return  TrelErr;
}
double getMaxRelErr(double *relErr, uns Np) {
  double maxRel = relErr[0];
  for(uns i = 1; i < Np; i++)
    if(relErr[i] > maxRel)
      maxRel = relErr[i];
  return maxRel;
}
