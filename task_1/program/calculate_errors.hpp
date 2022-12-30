/* В этом модуле записаны прототипы фукнций, используемых для
вычисления погрешностей.
Получается 4 вида информации о погрешности вычисления.
1)Абсолютная погрешность для каждого временного слоя.
2)Относительная погрешность для каждого временного слоя.
3)Значение максимальной абсолютной погрешности по временному слою.
4)Значение максимальной относительной погрешности по временному слою.
Эти величины записываются в специальные текстовые файлы.
1,2 - в файл errors, 3,4 - в файл MaxErrors. */

#ifndef __CALCULATE_ERRORS__
#define __CALCULATE_ERRORS__
#include "my_types.hpp"
#include <cmath>
void checkMemAlloc(void *pointer, const char *ErrMsg);

//Абсолютная погрешность для каждого временного слоя.
double * getAbsErr(double *curT_exact, double *curT, uns Np);

//Значение максимальной абсолютной погрешности по временному слою.
double getMaxAbsErr(double *curAbsErr, uns Np);

//Относительная погрешность для каждого временного слоя.
//Погр. выч. с исп. усреднённого значения точного решения.
double *getRelErr(double *curT_exact, double *curT,
                  const discrMesh &mesh, uns time);

//Знач. максимальной относительной погрешности по временному слою.
double getMaxRelErr(double *curRelErr, uns Np);
#endif
