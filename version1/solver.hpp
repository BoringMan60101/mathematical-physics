//В дальнейшем возможности могут быть расширены. Например коэффициенты могут перестать юыть константами.
//double * TDMA(double *T0, struct Q mesh, struct constPhysProperties consts, unsigned t); //Возможно на будущее
#ifndef __CV_SOLVER__
#define __CV_SOLVER__
#include "my_types.hpp"
void checkMemAlloc(void *pointer, const char *ErrMsg);
double ** getExactSolution(const discrMesh &mesh);
double * TDMA(double *T0, const discrMesh &mesh, const physConsts &constsm, uns t);
#endif
