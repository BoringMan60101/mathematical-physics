//Реализация функций, использующихся для разбиения [0; L] по координатам и по времени (от 0 до T_end)
/*
В данном файле записаны реализации функций, необходимых для формирования
одномерного дискретного аналога фазового пространства. (Разбиение области на конт. объёмы и временные слои)
*/
#include "mesh.hpp" //Подключение объявлений функций, чьи реализации идут ниже в файле
#include "my_types.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>

//Смещение по индексу на 1
double * getUniformFaces(double len, uns Nv) {
  double *xf = new double [Nv+2];
  checkMemAlloc(xf, "Error with memory allocation (in function - \'getUniform_xf\')");
  double dx = len/Nv;
  xf[0] = NAN; xf[Nv+1] = len;
  for(uns i = 1; i <= Nv; i++)
    xf[i] = (i-1)*dx;
  xf[Nv+1] = len;
  return xf;
}

//Смещение по индексу на 2
double * getFacesIntervals(double *xf, uns Np) {
  double *dxf = new double [Np];
  checkMemAlloc(xf, "Error with memory allocation (in function - \'get_dxf\')");
  dxf[0] = dxf[1] = NAN;
  for(uns i = 2; i < Np; i++)
    dxf[i] = xf[i] - xf[i-1];
  return dxf;
}

//Массив узловых точек
double * getUniformNodes(double *xf, double len, unsigned Np) {
  double *x = new double [Np];
  checkMemAlloc(xf, "Error with memory allocation (in function - \'getUniform_x\')");
  x[0] = 0.0f; x[Np-1] = len;
  for(uns i = 1; i <= Np-2; i++)
    x[i] = (xf[i+1] + xf[i])/2.0f;
  return x;
}

//Смещение на 1 по индексу
double * getNodesIntervals(double *x, uns Np) {
  double *dx = new double [Np];
  checkMemAlloc(dx, "Error with memory allocation (in function - \'get_dx\')");
  dx[0] = NAN;
  for(uns i = 1; i < Np; i++)
    dx[i] = x[i] - x[i-1];
  return dx;

}

//Получаем равномерный (пока) массивы временных интервалов
double * getUniformArray_dt(double T_end, uns Nt) {
  double *dt = new double [Nt+1];
  checkMemAlloc(dt, "Error with memory allocation (in function - \'getUniform_dt\')\n");
  double UnsformTimeStep = T_end/Nt;
  for(uns i = 0; i <= Nt; i++)
    dt[i] = UnsformTimeStep;
  return dt;
}

//Формирует массив, хранящий значения времени
double * getArray_cur_t(double *dt, uns Nt) {
  double *cur_t = new double [Nt+1];
  checkMemAlloc(cur_t, "Error with memory allocation (in function - \'get_cur_t\')\n");
  for(uns i = 0; i <= Nt; i++)
    cur_t[i] = i*dt[i];
  return cur_t;
}

void getUnformMesh(discrMesh &mesh) {
  mesh.Faces = getUniformFaces(mesh.len, mesh.Nv);
  mesh.FacesIntervals = getFacesIntervals(mesh.Faces, mesh.Np);
  mesh.Nodes = getUniformNodes(mesh.Faces, mesh.len, mesh.Np);
  mesh.NodesIntervals = getNodesIntervals(mesh.Nodes, mesh.Np);
  mesh.dt = getUniformArray_dt(mesh.T_end, mesh.Nt);
  mesh.cur_t = getArray_cur_t(mesh.dt, mesh.Nt);
}

void freeMesh(discrMesh &mesh) {
  delete [] mesh.Faces;
  delete [] mesh.FacesIntervals;
  delete [] mesh.Nodes;
  delete [] mesh.NodesIntervals;
  delete [] mesh.cur_t;
  delete [] mesh.dt;
}
