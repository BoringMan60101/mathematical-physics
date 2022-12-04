#ifndef __CALCULATE_ERRORS__
#define __CALCULATE_ERRORS__
#include "my_types.hpp"
double * getAbsErr(double **T_exact, double *T, uns Np, uns t);
double getMaxAbsErr(double *absErr, uns Np);
double getT_exactAvg(double time_layer, double len);
double * getRelErr(double **T_exact, double *T, uns Np, double time_layer, uns t, double len);
double getMaxRelErr(double *relErr, uns Np);
#endif
