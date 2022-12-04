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


void CheckMemAlloc(void *pointer, const char *ErrMsg);
void mkDirSafely(const char *dirName) {
  if(mkdir(dirName, 0777) == -1) {
    printf("Error with making directory \'%s\'\n", dirName);
    delete dirName;
    exit(EXIT_FAILURE);
  }
}
void freeMemory(phaseSpace &phaseSpace, double **T_exact);

int main() {
  //Подготовка директорий для вывода результатов
  const unsigned strMaxSize = 128;
  char *caseDirName = new char [strMaxSize];
  const char *dataDirName = "data";
  const char *errorsDirName = "errors";
  CheckMemAlloc(caseDirName, "Error with memory allocation in function \'main\'");
  printf("Enter a name of a directory, where results will be saved\n");
  scanf("%s", caseDirName);
  mkDirSafely(caseDirName); //Создание основной директории
  if(chdir(caseDirName) == -1) {
    printf("Unexpected error with changing directory in function \'main\'\n");
    return 1;
  }
  mkDirSafely(dataDirName);
  mkDirSafely(errorsDirName);

  //Основной блок объявления и инициализации переменных
  //В этой структуре хранятся физические параметры задачи (ro, c, k) - пока как константы
  unsigned Nv, Np, Nt = 5;
  double len = M_PI; //Длина исследуемого стержня
  double T_end; //T_end - конечное время исследования.
  constPhysProperties consts;
  printf("Enter number of control volumes --> "); scanf("%u", &Nv);
  printf("Enter time boundary (T_end) --> "); scanf("%lf", &T_end);
  printf("Enter number of time layers --> "); scanf("%u", &Nt);
  printf("Enter constant value of \'ro\' --> "); scanf("%lf", &(consts.ro));
  printf("Enter constant value of \'c\' --> "); scanf("%lf", &(consts.c));
  printf("Enter constant value of \'k\' --> "); scanf("%lf", &(consts.k));
  Np = Nv + 2;
  //printf("Enter length of stick --> "); scanf("%lf", &len);

  //инициализация дискретного аналога фазового пространства
  //(массивы: граней, узловых точек, временных словёв)
  phaseSpace discrSpace;
  getUnformMesh(discrSpace, len, Nv, T_end, Nt);


  errsArrays errorsArrays; //В этой структуре хранятся все массивы нужны для хранения ошибок вычислений
  //T0 - массив температуры на текущем временном слое, T-на следующем,
  //T_exact - массив значений точного решения во всех точках дискретной фазазовой области.
  double *T0 = new double [Np]; double *T = new double [Np], **T_exact = new double *[Nt+1];
  CheckMemAlloc(T0, "Mem Error! In function \'Main\'(array-T0)\n");
  CheckMemAlloc(T, "Mem Error! In function \'Main\'(array-T)\n");
  CheckMemAlloc(T_exact, "Mem Error! In function \'Main\'(array-T)\n");
  T_exact = getExactSolution(discrSpace.Nodes, discrSpace.Np, discrSpace.cur_t, discrSpace.Nt);
  for(uns i = 0; i < Np; i++)
     T0[i] = T_exact[0][i];

  //Вычисление значения температуры в контрольных точках на каждом временном слое от 0 до T_end с шагом dt
  //Основной цикл рассчёта температур в контрольных точках по всем временным слоям
  for(uns t = 1; t <= Nt; t++) {
    //Вычисление значений темепературы на текущем слое.
    //СЛАУ полученная после дискретизации дифференциального уравнения решается методом ТДМА
    T = TDMA(T0, discrSpace, t, consts);

    for(uns i = 1; i <= Np - 2; i++)  //Смена значений ("переход" к новому слою)
      T0[i] = T[i];

    //Подсчёт погрешностей на данном слое
    errorsArrays.abs = getAbsErr(T_exact, T, discrSpace.Np, t);
    errorsArrays.maxAbs = getMaxAbsErr(errorsArrays.abs, discrSpace.Np);
    errorsArrays.rel = getRelErr(T_exact, T, discrSpace.Np, discrSpace.cur_t[t], t, discrSpace.len);
    errorsArrays.maxRel = getMaxRelErr(errorsArrays.rel, discrSpace.Np);

    //Вывод данных о температуре в файлы (каждому временному слою одноименный файл).
    char *fname = new char [strMaxSize]; //имя файла куда будут записаны данные (соответсвтует временному слою)
    char *data = new char [strMaxSize]; //преобразованное в текст вещественное число равное значению темепературы в узловой точке
    CheckMemAlloc(fname, "Mem Error! In function \'Main\'\n");
    CheckMemAlloc(data, "Mem Error! In function \'Main\'\n");

    //переход в поддиреторию data
    FILE *fp = nullptr;
    if(chdir(dataDirName) != -1) {
      //Вывод данных в файлы о температуре на каждом временном слое
      //На один слой - один файл
      sprintf(fname, "%lf", discrSpace.cur_t[t]); //Получено имя для очередного файла
      fp = fopen(fname, "w");
      if(fp == nullptr) {
        printf("Error with opening \"%s\"\n", fname);
        return 1;
      }
      fprintf(fp, " X    T[x]    T_exact[x]\n");
      for(uns j = 0; j < discrSpace.Np; j++) {
        sprintf(data, "%lf %lf %lf", discrSpace.Nodes[j], T[j], T_exact[t][j]); //Перевели несколько числовых значений в строку
        fprintf(fp, "%s\n", data);
      }
      fclose(fp);
      //возврат обратно в главную директорю, потом переход в поддиректорю errors для записи погрешностей
      if(chdir("..") == -1) {
        printf("Error with changing directory (from \"%s\" to main \"%s\")\n", dataDirName, caseDirName);
        freeMemory(discrSpace, T_exact);
        delete [] caseDirName;
        return 1;
      }
    }
    else {
      printf("Error with changing directory from main \"%s\" to \"%s\"\n", caseDirName, dataDirName);
      delete [] T; delete [] errorsArrays.abs; delete [] errorsArrays.rel; delete [] fname; delete [] data;
        freeMemory(discrSpace, T_exact);
        delete [] caseDirName;
      return 1;
    }


    if(chdir(errorsDirName) != -1) { //вывод погрешностей по каждому временному слою
      //На один слой-один файл
      sprintf(fname, "%lf", discrSpace.cur_t[t]); //Получено имя для очередного файла
      fp = fopen(fname, "w");
      if(fp == nullptr) {
        printf("Error with opening file \"%s\"\n", fname);
        return 1;
      }
      fprintf(fp, "X    AbsErr[x]    RelErr[x]\n");
      for(uns j = 0; j < discrSpace.Np; j++) {
        sprintf(data, "%lf %lf %lf", discrSpace.Nodes[j], errorsArrays.abs[j], errorsArrays.rel[j]); //Перевели несколько числовых значений в строку
        fprintf(fp, "%s\n", data);
      }
      fclose(fp);

      //Отдельный файл, где для каждого слоя по одному числу - наибольшей погрешности (в обоих видах)
      fp = fopen("MaxErrors", "a+");
      if(fp == nullptr) {
        printf("Error with opening \"MaxErrors\"\n");
        return 1;
      }
      if(t == 1) fprintf(fp, "time  MaxAbsErr  MaxRelErr\n");
      sprintf(data, "%lf %lf %lf", discrSpace.cur_t[t], errorsArrays.maxAbs, errorsArrays.maxRel); //Перевели несколько числовых значений в строку
      fprintf(fp, "%s\n", data);
      fclose(fp);

      //возврат обратно в главную директорю, для следующего цикла записей
      if(chdir("..") == -1) {
        printf("Error with changing directory from main \"%s\" to main \"%s\"\n", errorsDirName, caseDirName);
        freeMemory(discrSpace, T_exact);
        delete [] caseDirName;
        return 1;
      }
    }
    else {
      printf("Error with changing directory from main \"%s\" to \"%s\"\n", caseDirName, errorsDirName);
      delete [] T; delete [] errorsArrays.abs; delete [] errorsArrays.rel; delete [] fname; delete [] data;
      freeMemory(discrSpace, T_exact);
      delete [] caseDirName;
      return 1;
    }

    //Освобождение памяти, поскольку данные в этих массивах больше не нужны
    delete [] T; delete [] errorsArrays.abs; delete [] errorsArrays.rel;
  }
  freeMemory(discrSpace, T_exact);
  delete [] caseDirName;
  return 0;
}

void freeMemory(phaseSpace &phaseSpace, double **T_exact) {
  delete [] phaseSpace.Faces;
  delete [] phaseSpace.FacesIntervals;
  delete [] phaseSpace.Nodes;
  delete [] phaseSpace.NodesIntervals;
  delete [] phaseSpace.cur_t;
  delete [] phaseSpace.dt;
  for(uns t = 0; t <= phaseSpace.Nt; t++)
    delete [] T_exact[t];
  delete T_exact;
}
void mkDirsSafely(const char *dirName) {
  if(mkdir(dirName, 0777) == -1) {
    printf("Error with making directory for saving data files\n");
    delete dirName;
    exit(EXIT_FAILURE);
  }
}
void CheckMemAlloc(void *pointer, const char *ErrMsg) {
  if(nullptr == pointer) {
    fprintf(stderr, "%s\n", ErrMsg);
    exit(EXIT_FAILURE);
  }
}
