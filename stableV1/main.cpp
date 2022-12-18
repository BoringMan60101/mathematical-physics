#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <sstream>
#include <sys/stat.h>

#include "my_types.hpp" //Описание используемых структур данных

//Функции для обработки входных параметров задачи
//и для вывода результатов в файлы
#include "IO_functions.hpp"

#include "mesh.hpp" //Процедуры для создания расчётной сетки

//Реализация ТДМА (для однродной задачи Дирихле со стационарным источником),
//Функция заполняющая сеточные значения эталонного решения
#include "solver.hpp"

#include "calculate_errors.hpp" //Функции для вычисления погрешностей

//Проверяет корректность динамического выделения памяти
//В случае ошибки выводит "отладочное сообщение" ErrMsg
//И завершает работу программы
void checkMemAlloc(void *pointer, const char *ErrMsg) {
  if(nullptr == pointer) {
    fprintf(stderr, "%s\n", ErrMsg);
    exit(EXIT_FAILURE);
  }
}

//Создаёт в текущей директории каталог с именем, указанным пользователем
//в качестве первого аргумента при запуске программы
//В созданном каталоге будут текстовые файлы с результатами работы программы
//Они содержат данные для визуализации, значения погрешностей.
void mkDirSafely(const char *dirName) {
  if(mkdir(dirName, 0777) == -1) {
    printf("Error with making directory \'%s\'\n", dirName);
    delete dirName;
    exit(EXIT_FAILURE);
  }
}

int main(const int argc, const char **argv) {
  //Проверка правильного запуска программы
  //Обязательные аргументы при запуске: имя директории, имя файла со входными данными
  if(argc != 3) {
    printf("Error! Give two arguments: outputDirName, 2-initValuesFileName\n");
    return 1;
  }

  //Описание структур в заголовочном файле my_types.hpp
  discrMesh mesh; //Хранит все параметры сетки
  physConsts consts; //Хранит физические константы (ro, c, k)
  errors errs; //Хранит массивы с погрешностями

  //Получение исходных данных для формирования сетки, и физических констант (ro, c, k)
  getCaseParameters(argv[2], consts, mesh);

  //Создание сетки с заданными параметрами (формируются массивы координат и временных слоёв)
  getUniformMesh(mesh);

  //Создание директории для вывода результатов расчёта
  mkDirSafely(argv[1]);

  //Переход в созданную директорию
  if(chdir(argv[1]) == -1) {
    printf("Error with changing directory to \"%s\"\n", argv[1]);
    return 1;
  }

  //T_exact - набор значений аналитеского решения в узловых точках по всем временным слоям
  //T0 - массив температур (известных) в узловых точках на текущем временном слое
  //T - массив температур (неизвестных) в узловых точках на следующем временном слое
  double **T_exact, *T0 = new double [mesh.Np], *T = new double [mesh.Np];
  checkMemAlloc(T0, "Error with memory allocation \"(T0) main\"\n");
  checkMemAlloc(T, "Error with memory allocation \"(T) main\"\n");
  T_exact = getExactSolution(mesh);
  for(uns i = 0; i < mesh.Np; i++)
    T0[i] = T_exact[0][i];

  //Основной цикл состоит из двух этапов:
  //1) Рассчёт температуры на каждом слое;
  //2) Вывод данных и погрешностей в файлы
  for(uns t = 1; t <= mesh.Nt; t++) {
    //Вычисление температуры на текущем временном слое с помощью алгоритма TDMA
    T = TDMA(T0, mesh, consts, t);

    //Вычисление погрешностей для данного временного слоя
    errs.abs = getAbsErr(T_exact[t], T, mesh.Np);
    errs.maxAbs = getMaxAbsErr(errs.abs, mesh.Np);
    errs.rel = getRelErr(T_exact[t], T, mesh, t);
    errs.maxRel = getMaxRelErr(errs.rel, mesh.Np);

    //Вывод значений температуры и погрешностей в соответсвующие файлы
    printDataToFile(mesh, T_exact[t], T);
    printErrorsOnTimeLayerToFile(mesh, errs, t);
    printMaxErrorsToFile(mesh, errs, t);

    //Обновление значений - переход к новому слою
    for(uns i = 0; i < mesh.Np; i++)
      T0[i] = T[i];
    delete [] T; delete [] errs.abs; delete [] errs.rel;
  }

  //Очистка выделенной памяти
  for(uns t = 0; t < mesh.Nt; t++)
    delete [] T_exact[t];
  delete [] T_exact;
  delete [] T0;
  freeMesh(mesh);
  return 0;
}
