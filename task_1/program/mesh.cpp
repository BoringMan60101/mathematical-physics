/* В данном файле выполнены реализации функций,
необходимых для формирования одномерного дискретного
аналога фазового пространства - расчётной сетки.
В такой реализации сетка равномерная.
(Как по координатам, так и по времени.)

Все массивы относящиеся к пространству имеют одинаковую длину,
которая равна Np - общему числу узловых точек.
В массивах:
Faces, FacesIntervals, dFaces, dNodes - фактических значений меньше
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
  checkMemAlloc(xf, "Error in function - getUniform_xf\n");
  double dx = len / Nv;
  xf[0] = NAN;
  xf[Nv+1] = len;
  for(uns i = 1; i <= Nv; i++)
    xf[i] = (i-1)*dx;
  xf[Nv+1] = len;
  return xf;
}


//Массив узловых точек
//Узел помещается в посередине между гранями
double * getUniformNodes(double *xf, double len, unsigned Np) {
  double *x = new double [Np];
  checkMemAlloc(xf, "Error in function - getUniform_x\n");
  x[0] = 0.0f; x[Np-1] = len;
  for(uns i = 1; i <= Np-2; i++)
    x[i] = (xf[i+1] + xf[i]) / 2.0f;
  return x;
}

//Создаёт равномерный массив временных интервалов
double * getUniformArray_dt(double T_end, uns Nt) {
  double *dt = new double [Nt+1];
  checkMemAlloc(dt, "Error in function getUniform_dt\n");
  double UnsformTimeStep = T_end / Nt;
  for(uns i = 0; i <= Nt; i++)
    dt[i] = UnsformTimeStep;
  return dt;
}

//Формирует массив, хранящий значения времени
double * getArray_cur_t(double T_end, uns Nt) {
  double *cur_t = new double [Nt+1];
  checkMemAlloc(cur_t, "Error in function get_cur_t\n");
  double h = T_end / Nt;
  for(uns i = 0; i <= Nt; i++)
    cur_t[i] = i*h;
  return cur_t;
}

//В этой функции заполняются все массивы, которые
//представляют рассчётную сетку
void getUniformMesh(discrMesh &mesh) {
  mesh.Faces = getUniformFaces(mesh.len, mesh.Nv);
  mesh.Nodes = getUniformNodes(mesh.Faces, mesh.len, mesh.Np);;
  mesh.dt = mesh.T_end / mesh.Nt;
  mesh.cur_t = getArray_cur_t(mesh.T_end, mesh.Nt);
}

void freeMesh(discrMesh &mesh) {
  delete [] mesh.Faces;
  delete [] mesh.Nodes;
  delete [] mesh.cur_t;
}
