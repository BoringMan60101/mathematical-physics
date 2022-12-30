/* Данный модуль отвечает за вычислительную часть программы.
В нём реализован алгоритм трёдиагональной матрицы для случая
однородной задачи Дирихле, когда известно аналитическое выражение
стационарного источника. */

#include <cstdlib>
#include <cstdio>
#include "solver.hpp"
#include "my_types.hpp"
#include <cmath>

void checkMemAlloc(void *pointer, const char *ErrMsg);

double **getExactSolution(const discrMesh &mesh) {
//Массив хранит знач. аналитич. решения во всех рассчётных точках
  double ** T_exact = new double *[mesh.Nt+1];
  checkMemAlloc(T_exact, "Erorr in function T_exact\n");
  for(uns t = 0; t <= mesh.Nt; t++) {
    T_exact[t] = new double [mesh.Np];
    checkMemAlloc(T_exact[t], "Erorr in function T_exact\n");
    T_exact[t][0] = 0.0f; //Однородные граничные условия
    T_exact[t][mesh.Np-1] = 0.0f; //Однородные граничные условия
    for(uns i = 1; i <= mesh.Np-2; i++)
      T_exact[t][i] = (1.0f-exp(-mesh.cur_t[t]))*sin(mesh.Nodes[i]);
  }
  return T_exact;
}

double *TDMA(double *T0,const discrMesh &mesh,const physConsts &consts)
{
  //mesh - Структура, где хранятся все данные расчётной сетки
  //consts - Структура, где хранятся физические константы ro, c, k
  //t - Значение очередного времееного слоя
  double dxW; //Расстояние между узловыми точками (P-W)
  double dxE; //Расстояние между узловыми точками (E-P)
  double dFaces; //Размер очередного контр. объема (e-w)
  double fW, fE; //Координаты западной и восточной граней
  double aW, aP, aE, bP; //Коэф. уравнений на кажд. врем. слое
  //Прямой ход метода ТДМА заполняет массивы P, Q.
  //На этапе обратного хода с их помощью вычисляется температура
  //в каждой узловой точке
  double *P, *Q;
  double *T;
  T = new double [mesh.Np];
  P = new double [mesh.Np];
  Q = new double [mesh.Np];
  if(!(P && Q && T)) {
    fprintf(stderr, "%s\n", "Mem Error!\n In function \'TDMA\'\n");
    exit(EXIT_FAILURE);
  }

  //Расстояние между узловыми точками на равномерной сетке одинаково
  //В данном случае можно вычислить коэффициенты один раз.
  dxW = dxE = mesh.Nodes[1] - mesh.Nodes[0];
  aW = aE = (0.5f*consts.k) / dxW;
  //Прямой ход
  T[0] = P[0] = Q[0] = 0.0f; //Граничное условие первого рода
  for(uns i = 1; i <= mesh.Np-2; i++) {
    //Вычисляем P[i], Q[i] для внутренних объёмов
    //dFaces - Расст. между гранями очередного контрольного объема
    dFaces = mesh.Faces[i+1] - mesh.Faces[i];
    aP = (consts.ro*consts.c*dFaces)/mesh.dt + aW + aE;
    fW = mesh.Faces[i]; //Западная грань
    fE = mesh.Faces[i+1]; //Восточная грань
    //bP - Вклад из левой части
    bP = (consts.ro*consts.c*dFaces*T0[i])/mesh.dt;

//Вклад после интегрирования конвективного члена
//При интегрировании по времени исп. - схема Кранка-Николсона
    bP += aE*(T0[i+1] - T0[i]) - aW*(T0[i] - T0[i-1]);
    bP += cos(fW) - cos(fE); //Вклад стационарного источника
    P[i] = aE / (aP - aW*P[i-1]); //Заполнение массивов P, Q
    Q[i] = (aW*Q[i-1] + bP) / (aP - aW*P[i-1]);
  }
  //Задано однородное граничное условие первого рода
  T[mesh.Np-1] = P[mesh.Np-1] = Q[mesh.Np-1] = 0.0f;

  //Обратный ход (определение температуры во внутренних точках)
  for(uns i = mesh.Np-2; i >= 1; i--)
    T[i] = P[i]*T[i+1] + Q[i];

  delete [] P; delete [] Q;
  return T;
}
