1. Реализовать LU-разложение матрицы A c выбором ведущего
элемента по по столбцу или по всей матрице. Проверить разложение
сравнением матриц LU и PA (или PAQ), где P — матрица перестановки
строк, а Q — столбцов. Выполнить для системы произвольной
размерности, генерировать случайную матрицу для демонстрации
работы программы.

С использованием LU-разложения найти:
a) Определитель матрицы A;
b) Решение СЛАУ Ax = b, выполнить проверку равенства Ax − b = 0;
c) Матрицу A −1 (выполнить проверку AA −1 и A −1 A);
d) Число обусловленности матрицы A.


2. Модифицировать алгоритм для нахождения ранга вырожденных
матриц: при выборе ведущего элемента только по столбцу надо
приводить матрицу к ступенчатой форме, что потребует обнулять
элементы не только под диагональю; при выборе ведущего элемента
по всей матрице никаких изменений не требуется — в этом случае
матрица приводится к трапецевидной форме). Проверять так же систему с вырожденной матрицей
на совместность и выдавать любое частное решение, если она
совместна.


3. Реализовать QR-разложение матрицы A. Проверить разложение
перемножением матриц Q и R. С его помощью найти решение
невырожденной СЛАУ Ax = b.


4. Реализовать метод Якоби и метод Зейделя решения СЛАУ.
Сравнить на примере СЛАУ с матрицей с диагональным преобладанием
и с положительно определённой матрицей без диагонального
преобладания (генерировать случайные матрицы по размерности).
Дать априорную оценку числа необходимых итераций, сравнить с
апостериорной оценкой.