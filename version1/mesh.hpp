/*
Заголовочный файл (содержит только объявления объектов).
Данный файл содержит объявления (прототипы) функций,
которые используются для построения дискретного аналога фазового пространства.
(одномерного в случае со стрежнем)
*/
#ifndef __MESH_HPP__
#define __MESH_HPP__
#include "my_types.hpp"
void checkMemAlloc(void *pointer, const char *ErrMsg);
double * getUniformFaces(double len, uns Nv);
double * getFacesIntervals(double *xf, uns Np);
double * getUniformNodes(double *xf, double len, uns Np);
double * getNodesIntervals(double *x, uns Np);
double * getUniformArray_dt(double T_end, uns Nt);
double * getArray_cur_t(double *dt, uns Nt);
void getUnformMesh(discrMesh &mesh);
void freeMesh(discrMesh &mesh);
#endif
