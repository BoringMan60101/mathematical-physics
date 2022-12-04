//В дальнейшем возможности могут быть расширены. Например коэффициенты могут перестать юыть константами.
//double * TDMA(double *T0, struct Q mesh, struct constPhysProperties consts, unsigned t); //Возможно на будущее
#ifndef __CV_SOLVER__
#define __CV_SOLVER__
#include "my_types.hpp"
void CheckMemAlloc(void *pointer, const char *ErrMsg);
double ** getExactSolution(double *x, uns Np, double *cur_t, uns Nt);
double * TDMA(double *T0, phaseSpace discrSpace, uns t, constPhysProperties Consts); //Этот пока испо
//double * TDMA(double *T0, double *xf, double *x, uns Np, double *dt, uns t, double ro, double c, double k); //Этот пока используется
#endif
