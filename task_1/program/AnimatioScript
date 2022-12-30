#Установка режима анимации со скоростью 2 кадра в сек.
set terminal gif animate delay 50 

#Выбор имени выходного файла
set output "HeatTransAnim.gif"

set font "Times New Roman,12" #Шрифт для подписей
set xlabel "x" #Подпись оси Ox
set xrange [0:pi] #Диапозон изменения x
set ylabel "T(x)" #Подпись оси Oy
set yrange [0:1] #Диапозон изменеия y = T(x)
set grid 
stats "data.dat" name "File" #Получение информации о файле с данными
#%set nokey #Отключение легенды
set title "Evalution of temperature" #Подпись к анимации

#Цикл для покадрового построения анимации
do for [i=1:int(File_blocks)] { \
plot "data.dat" \
using 1:2 index i with linespoints \
linewidth 2 \
linecolor "red" \
title sprintf("frame = %d", i) }
#Анимация получается следующим образом:
#Выбраются данные из первых двух столбцов файла.
#Первый столбец это коорд. узловой точки, а второй знач. темп.
#Выбранный тип линий такой, что два соседних значения температуры
#в узловых точках соединяются прямой.
#Толщина линии равна 2 единицам.
#Цвет линии - красный.
#Опция title выводит номер кадра

