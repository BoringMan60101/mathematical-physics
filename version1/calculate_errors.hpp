#ifndef __CALCULATE_ERRORS__
#define __CALCULATE_ERRORS__
#include "my_types.hpp"
#include <cmath>
void checkMemAlloc(void *pointer, const char *ErrMsg);
double * getAbsErr(double *curT_exact, double *curT, uns Np);
double getMaxAbsErr(double *curAbsErr, uns Np);
double T_exactFormula(const double time_layer, const double x);
double getAverageOfCurT_exact(double time_layer, double len);
double * getRelErr(double *curT_exact, double *curT, double time_layer, uns Np, double len);
double getMaxRelErr(double *curRelErr, uns Np);
#endif
