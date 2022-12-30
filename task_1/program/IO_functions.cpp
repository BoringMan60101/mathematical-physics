/* В данном файле реализованы функции для обработки
файла с иходными данными задачи.
А также функции для вывода результатов в файлы */

#include "my_types.hpp"
#include <cstdlib>
#include <cstdio>

void checkMemAlloc(void *pointer, const char *ErrMsg);

//Вспомогательная функция, используемая при обработке вход.о файла.
//Она определяет чем является симвлом, цифрой или нет.
bool isDigitOrSign(const char ch) {
  char digits[13] = { '-', '+', '.', '0', '1',
                      '2', '3', '4', '5', '6',
                      '7', '8', '9'};
  bool result = false;
  for(int i = 0; i < 13; i++)
    if(ch == digits[i]) {
      result = true;
      break;
    }
  return result;
}

//В этой процедуре реализовано заполнение входных парам. задачи
//из заданного пользователем файла.
void getCaseParameters(const char *fileName,
                       physConsts &consts, discrMesh &mesh) {
  uns maxSize = 256; //Максимальная длина для считывания строки
  uns linesCounter = 0; //Счётчик числа считанных строк
  uns paramsNum = 7; //Сколько параметров должно быть файле

  FILE *fp = fopen(fileName, "r");
  if(fp == nullptr) {
    printf("Error with opening file \"%s\"\n", fileName);
    exit(EXIT_FAILURE);
  }

  //Обработка файла происходит построчно
  char *line = new char [maxSize];
  checkMemAlloc(line, "Error in function getCaseParameters\n");
  char *numStr; //В неё заносится число в виде набора символов


  while(1) {
    //fgets() - возвращает nullptr в двух случаях!
    //1)если содерж. файла закончилось и считывать больше нечего
    //2)если возникла ошибка считывания
    if(fgets(line, maxSize-1, fp) == nullptr) {
      //Считывание прервано
      //Либо файл закончился, либо неправильный файл
      if(linesCounter != paramsNum) {
      //Если в файле не верное число строк,
      //то в этом случае программа завершается.
        printf("Error with reading file \"%s\"!\n", fileName);
        exit(EXIT_FAILURE);
      }
      break;
    }

    linesCounter += 1; //Увеличение счётчика считанных строк
    int j = 0; //Индекс символа в текущей считанной строке

//В массив numStr будут записаны численые знач. параметров посимвольно
    numStr = new char [maxSize];
    checkMemAlloc(numStr, "Error in func \"getCaseParameters\"\n");

    //Отсечение лишних символов, до первой цифры или знака (+/-/.)
    while(isDigitOrSign(line[j]) == false)
      j++;


    int i = 0; //Индекс символа в массиве numStr
    //j - индекс символа в считанной строке line
    while(line[j]) {
      if(line[j] == '\n') {
        numStr[i] = '\0';
        break;
      }
      numStr[i] = line[j];
      i++; j++;
    }

    switch (linesCounter) {
//Заполнение соответсвующих полей структур.
//В такой реализации входной файл должен быть оформлен в спец. виде.
//Чтобы записанные в нём значения вводились в подходящие поля.
      case 1: { mesh.len = atof(numStr); break; }
      case 2: {
        mesh.Nv = atoi(numStr);
        mesh.Np = mesh.Nv + 2;
        break;
      }
      case 3: { mesh.T_end = atof(numStr); break; }
      case 4: { mesh.Nt = atoi(numStr); break; }
      case 5: { consts.ro = atof(numStr); break; }
      case 6: { consts.c = atof(numStr); break; }
      case 7: { consts.k = atof(numStr); break; }
//Если строк больше, чем предполагается, то обработка завершается
      default: {
        delete [] line;
        delete [] numStr;
        fclose(fp);
        printf("Error with reading file \"%s\"!\n", fileName);
        exit(EXIT_FAILURE);
       }
    }
    delete [] numStr;
  }
  delete [] line;
  fclose(fp);
}


//Функция для вывода расчётных данных файл "data.dat"
//Они используются для визуализации результатов.
void printDataToFile(const discrMesh &mesh,
                     double *curT_exact, double *curT) {
  static bool first = true;
  FILE * fp = fopen("data.dat", "a+");
  if(fp == nullptr) {
    printf("Erorr with opening file to write data\n");
    exit(EXIT_FAILURE);
  }

  if(first) {
    //"Шапка", чтобы было понятно какая колонка за что отвечает
    fprintf(fp, "#   X      T      T_exact\n");
    first = false;
  }

  for(uns i = 0; i < mesh.Np; i++)
    fprintf(fp, "%lf  %lf   %lf\n", mesh.Nodes[i], curT[i], curT_exact[i]);
  //Разделение между данными на каждом слое в две пустых строки
  fprintf(fp, "\n\n");
  fclose(fp);
}

//Процедура, выводящая погрешности в каждой расчётной точке
//в отдельный файл "errors"
void printErrorsOnTimeLayerToFile(const discrMesh &mesh,
                                  const errors &errs, uns t) {
  static bool first = true;
  FILE * fp = fopen("errors", "a+");
  if(fp == nullptr) {
    printf("Erorr with opening file to write errors\n");
    exit(EXIT_FAILURE);
  }

  if(first) {
    //"Шапка", чтобы было понятно какая колонка за что отвечает
    fprintf(fp, "#   X  AbsErr  RelErr [%%]\n");
    first = false;
  }

  for(uns i = 0; i < mesh.Np; i++) {
    if (i == 0)
      fprintf(fp,"#time: %.3lf\n", mesh.cur_t[t]);
    fprintf(fp,"%.3lf   %.2lf   %.2lf\n", mesh.Nodes[i],errs.abs[i],errs.rel[i]);
  }
  fprintf(fp, "\n\n"); //Разделение между слоями в две пустые строки
  fclose(fp);
}

//Выводит максимальную погрешность на каждом слое в файл "MaxErrors"
void printMaxErrorsToFile(const discrMesh &mesh,
                          const errors &errs, uns t) {
  static bool first = true;
  FILE * fp = fopen("MaxErrors", "a+");
  if(fp == nullptr) {
    printf("Erorr with opening file to write errors\n");
    exit(EXIT_FAILURE);
  }

  if(first) {
    fprintf(fp, "#time  MaxAbs  MaxRel [%%]\n");
    first = false;
  }

  fprintf(fp, "%.3lf   %.2lf   %.2lf\n", mesh.cur_t[t], errs.maxAbs, errs.maxRel);
  fclose(fp);
}
