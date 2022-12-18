/* В данном файле выполнены реализации функций, необходимых для формирования
одномерного дискретного аналога фазового пространства. (Расчётной сетки)
В такой реализации сетка равномерная, как по координатам, так и по времени.

Все массивы относящиеся к пространству (координаты и расстояния между точками)
имеют одинаковую длину равную Np - общему числу узловых точек.
В массивах Faces, FacesIntervals, dFaces, dNodes - фактических значений меньше
Поэтому неиспользуемым ячейкам присваивается значение NAN. */

#include "mesh.hpp"
#include "my_types.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>

//Заполняет массив с координатами граней контрольных объёмов.
//Первое значение не используется
double * getUniformFaces(double len, uns Nv) {
  double *xf = new double [Nv+2];
  checkMemAlloc(xf, "Error with memory allocation (in function - \'getUniform_xf\')");
  double dx = len/Nv;
  xf[0] = NAN;
  xf[Nv+1] = len;
  for(uns i = 1; i <= Nv; i++)
    xf[i] = (i-1)*dx;
  xf[Nv+1] = len;
  return xf;
}


double * getFacesIntervals(double *xf, uns Np) {
  //В массиве размеров контрольных объёмов
  //сделано смещение по индексу на 2.
  double *dxf = new double [Np];
  checkMemAlloc(xf, "Error with memory allocation (in function - \'get_dxf\')");
  dxf[0] = dxf[1] = NAN;  //неиспользуемые ячейки
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
    x[i] = (xf[i+1] + xf[i])/2.0f; //Узел помещается в посередине между гранями
  return x;
}

//Смещение на 1 по индексу
double * getNodesIntervals(double *x, uns Np) {
  double *dx = new double [Np];
  checkMemAlloc(dx, "Error with memory allocation (in function - \'get_dx\')");
  dx[0] = NAN; //Спец. значение (не используется)
  for(uns i = 1; i < Np; i++)
    dx[i] = x[i] - x[i-1];
  return dx;

}

//Создаёт равномерный массив временных интервалов
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

//В этой функции заполняются все массивы, представляющие рассчётную сетку
void getUniformMesh(discrMesh &mesh) {
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
