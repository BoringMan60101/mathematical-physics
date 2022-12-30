#define _USE_MATH_DEFINES //Для доступа к константе - числу Пи
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <sys/stat.h>

//Подключение пользовательских структур данных
#include "my_types.hpp" //mesh, consts, errors, сокр. для unsigned

//Функции для обработки входных параметров задачи
//и для вывода результатов в файлы
#include "IO_functions.hpp"

#include "mesh.hpp" //Процедуры для создания расчётной сетки

//Реализация алгоритма трёх. диаг. матрицы - ТДМА
//В данном случае для однород. задачи Дирихле со стац. источником
//И вспомогат. функ. заполн. сеточ. знач. эталонного решения
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

//Создаёт в тек. директории каталог с именем, заданным польз.
//в качестве первого аргумента при запуске программы
//В созд. каталоге будут текст. файлы с результатами.
//Они содержат данные для визуализации, значения погрешностей.
void mkdirSafely(const char *newDir) {
  if(mkdir(newDir, 0777) == -1) {
    printf("Error with making directory \'%s\'\n", newDir);
    exit(EXIT_FAILURE);
  }
}

//Безопсный переход в директорию для записи результатов
void chdirSafely(const char *destDir) {
  if(chdir(destDir) == -1) {
    printf("Error with changing directory to \"%s\"\n", destDir);
    exit(EXIT_FAILURE);
  }
}

int main(const int argc, const char **argv) {
  //Проверка правильного запуска программы
  //Обязательные аргументы при запуске:
  //имя директории, имя файла со входными данными задачи:
  //Это длина, число объёмов, конеч. время, ro, c, k
  //Именно в таком порядке они должны быть в конфигурац. файле.
  if(argc != 3) {
    printf("Error! Need two arguments: \n");
    printf("1 --- Name of file with initial values\n");
    printf("2 --- Name of directory to save results\n");
    return 1;
  }

  discrMesh mesh; //Хранит все параметры сетки
  physConsts consts; //Хранит физические константы (ro, c, k)
  errors errs; //Хранит массивы с погрешностями

  //Получение исходных данных для формирования сетки:
  //длину, число контр об., кон. время, число врем. слоём
  //В этом же файле задаются изические константы (ro, c, k)
  getCaseParameters(argv[1], consts, mesh);

  //Создание сетки с заданными параметрами:
  //формируются массивы координат и временных слоёв.
  getUniformMesh(mesh);

  //Создание директории для вывода результатов расчёта
  mkdirSafely(argv[2]);

  //Переход в созданную директорию
  chdirSafely(argv[2]);

//T_exact - набор знач. аналит. реш. в узловых т. по всем врем. слоям
//T0 - массив темп. (известных) в узловых т. на тек. врем. слое
//T - массив темп. (неизвестных) в узловых т. на след. врем. слое
//Заполняем таблицу сеточных значений точного решения
  double **T_exact = getExactSolution(mesh);
  double *T0 = new double [mesh.Np];
  double *T = new double [mesh.Np];
  checkMemAlloc(T0, "Error with memory allocation \"(T0) main\"\n");
  checkMemAlloc(T, "Error with memory allocation \"(T) main\"\n");
  for(uns i = 0; i < mesh.Np; i++)
    T0[i] = T_exact[0][i];

  //Основной цикл состоит из двух этапов:
  //1) Рассчёт температуры на каждом слое;
  //2) Вывод данных и погрешностей в файлы
  for(uns t = 1; t <= mesh.Nt; t++) {
  //Вычисл. темп. на тек. врем. слое с помощью алгоритма TDMA
    T = TDMA(T0, mesh, consts);
    
    //Вычисление погрешностей для данного временного слоя
    errs.abs = getAbsErr(T_exact[t], T, mesh.Np);
    errs.maxAbs = getMaxAbsErr(errs.abs, mesh.Np);
    errs.rel = getRelErr(T_exact[t], T, mesh, t);
    errs.maxRel = getMaxRelErr(errs.rel, mesh.Np);

    //Вывод значений темп. и погр. в соответсвующие файлы
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
