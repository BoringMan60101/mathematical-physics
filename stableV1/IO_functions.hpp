/* Данный Заголовочный файл содержит объявления (прототипы) функций,
которые используются для получения исходных данных задачи
из файла, заданного пользователем. А также функции вывода
результатов расчёта и погрешностей в отдельные файлы */

#ifndef __IO_FUNCTIONS__
#define __IO_FUNCTIONS__
#include "my_types.hpp"

//Проверяет корректность выделения блока памяти
void checkMemAlloc(void *pointer, const char *ErrMsg);

//Вспомогательная функция, исп. при обработке файла с исх. данными задачи (кейса)
bool isNotDigitOrSign(const char ch);

//Получения значений физ. констант и параметров сетки (Nv, Tend, Nt, len)
//Файл с исходными данными задачи передаётся пользователем как аргумент
void getCaseParameters(const char *fileName, physConsts & consts, discrMesh & mesh);

//Выводит рассчитанные и точные значения температуры в файл "data" (по всем слоям в один)
//Между каждым временным слоем разделение в виде двух пустых строк
void printDataToFile(const discrMesh &mesh, double *curT_exact, double *curT);

//Выводит значения для каждого погрешностей в файл "errors" (по всем слоям в один)
void printErrorsOnTimeLayerToFile(const discrMesh &mesh, const errors &errs, uns t);

//Выводит максимальные значения погрешностей в файл "MaxErrors"
//(по одному значению с каждого временного слоя)
void printMaxErrorsToFile(const discrMesh &mesh, const errors &errs, uns t);
#endif
