/* В этом модуле объявляны функции,
выполняющие вычислительную часть программы. */
#ifndef __CV_SOLVER__
#define __CV_SOLVER__
#include "my_types.hpp"
void checkMemAlloc(void *pointer, const char *ErrMsg);

//Вычисляет сеточные значения точного решения.
double **getExactSolution(const discrMesh &mesh);

//Реализация алгоритма ТДМА
double *TDMA(double *T0,const discrMesh &mesh,const physConsts &consts);
#endif
