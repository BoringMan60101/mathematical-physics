#define _USE_MATH_DEFINES
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <cstring>
#include <fcntl.h>
#include <sstream>
#include <sys/stat.h>
#include <iostream>

#include "mesh.hpp"
#include "solver.hpp"
#include "calculate_errors.hpp"
#include "my_types.hpp"

void checkMemAlloc(void *pointer, const char *ErrMsg);
void mkDirSafely(const char *dirName) {
  if(mkdir(dirName, 0777) == -1) {
    printf("Error with making directory \'%s\'\n", dirName);
    delete dirName;
    exit(EXIT_FAILURE);
  }
}
void getMeshParametrs(discrMesh &mesh);
void getPhysicsParamets(physConsts &Consts);
void printDataToFile(const discrMesh &mesh, double *curT_exact, double *curT);
void printErrorsToFile(const discrMesh &mesh, const errors &errs, uns t);
void printMaxErrorsToFile(const discrMesh &mesh, const errors &errs, uns t);

int main(const int argc, const char **argv) {
  //Подготовка директорий для вывода результатов
  if(argc != 2) {
    printf("Error! Give on one argument - a name for directtory to save results\n");
    return 1;
  }
  //Создание директории для вывода результатов (температура и погрешности)
  mkDirSafely(argv[1]);
  if(chdir(argv[1]) == -1) { //Переход в созданную директорию
    printf("Error with changing directory to \"%s\"\n", argv[1]);
    return 1;
  }

  //Блок объявления переменных.
  //Описание структур в заголовочном файле my_types.hpp
  discrMesh mesh; //Структура для хранения всех параметров сетки
  physConsts consts; //Структура для хранения физических констант (ro, c, k)
  errors errs; //Структура для хранения массивов погрешностей

  //Формирование сетки; получение физических констант
  getMeshParametrs(mesh);
  getUnformMesh(mesh);
  getPhysicsParamets(consts);

  //T_exact -
  double **T_exact, *T0 = new double [mesh.Np], *T = new double [mesh.Np];
  checkMemAlloc(T0, "Error with memory allocation \"(T0) main\"\n");
  checkMemAlloc(T, "Error with memory allocation \"(T) main\"\n");
  T_exact = getExactSolution(mesh); //Заполнили сеточные значения точного решения
  for(uns i = 0; i < mesh.Np; i++)
    T0[i] = T_exact[0][i];

  //Основной цикл:
  //1) Рассчёт температуры на каждом слое;
  //2) Вывод данных и погрешностей в файлы
  for(uns t = 1; t <= mesh.Nt; t++) {
    //Вычисление температуры на текущем временном слое
    T = TDMA(T0, mesh, consts, t);

    //Получение погрешностей значения температуры для данного временного слоя
    errs.abs = getAbsErr(T_exact[t], T, mesh.Np);
    errs.maxAbs = getMaxAbsErr(errs.abs, mesh.Np);
    //Время для вычисления интеграла передаётся
    errs.rel = getRelErr(T_exact[t], T, mesh.cur_t[t], mesh.Np, mesh.len);
    errs.maxRel = getMaxRelErr(errs.rel, mesh.Np);

    //Вывод значений и их погрешностей в соответсвующие файлы
    printDataToFile(mesh, T_exact[t], T);
    printErrorsToFile(mesh, errs, t);
    printMaxErrorsToFile(mesh, errs, t);

    //Обновление значений - переход к новому слою
    for(uns i = 0; i < mesh.Np; i++)
      T0[i] = T[i];

    delete [] T; delete [] errs.abs; delete [] errs.rel;
  }

  //Очистка всей  остальной выделенной памяти
  for(uns t = 0; t < mesh.Nt; t++)
    delete [] T_exact[t];
  delete [] T_exact;
  delete [] T0;
  freeMesh(mesh);
  return 0;
}


void printDataToFile(const discrMesh &mesh, double *curT_exact, double *curT) {
  FILE * fp = fopen("data", "a+");
  if(fp == nullptr) {
    printf("Erorr with opening file to write data\n");
    exit(EXIT_FAILURE);
  }
  for(uns i = 0; i < mesh.Np; i++)
    fprintf(fp, "%lf %lf %lf\n", mesh.Nodes[i], curT_exact[i], curT[i]);
  fprintf(fp, "\n\n"); //Разделение между слоями в две пустых строки
  fclose(fp);
}
void printErrorsToFile(const discrMesh &mesh, const errors &errs, uns t) {
  FILE * fp = fopen("errors", "a+");
  if(fp == nullptr) {
    printf("Erorr with opening file to write errors\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "Time %lf \n", mesh.cur_t[t]);
  for(uns i = 0; i < mesh.Np; i++)
    fprintf(fp, "%lf %lf %lf\n", mesh.Nodes[i], errs.abs[i], errs.rel[i]);
  fprintf(fp, "\n"); //Разделение между слоями в две пустых строки
  fclose(fp);
}
void printMaxErrorsToFile(const discrMesh &mesh, const errors &errs, uns t) {
    FILE * fp = fopen("MaxErrors", "a+");
    if(fp == nullptr) {
      printf("Erorr with opening file to write MaxErrors\n");
      exit(EXIT_FAILURE);
    }
    fprintf(fp, "Time %lf MaxAbs %lf MaxRel %lf\n", mesh.cur_t[t], errs.maxAbs, errs.maxRel);
    fclose(fp);
}
void getMeshParametrs(discrMesh &mesh) {
  double len = M_PI; //Длина исследуемого стержня
  //printf("Enter length of stick --> "); scanf("%lf", &(mesh.len));
  mesh.len = len;
  printf("Enter number of control volumes --> ");
  scanf("%u", &(mesh.Nv));
  printf("Enter time boundary (T_end) --> ");
  scanf("%lf", &(mesh.T_end));
  printf("Enter number of time layers --> ");
  scanf("%u", &(mesh.Nt));
  mesh.Np = mesh.Nv + 2;
}
void getPhysicsParamets(physConsts &Consts) {
  printf("Enter constant value of \'ro\' --> ");
  scanf("%lf", &(Consts.ro));
  printf("Enter constant value of \'c\' --> ");
  scanf("%lf", &(Consts.c));
  printf("Enter constant value of \'k\' --> ");
  scanf("%lf", &(Consts.k));
}
void checkMemAlloc(void *pointer, const char *ErrMsg) {
  if(nullptr == pointer) {
    fprintf(stderr, "%s\n", ErrMsg);
    exit(EXIT_FAILURE);
  }
}
