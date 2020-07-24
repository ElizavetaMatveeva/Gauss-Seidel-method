// Написать программу, решающую СЛАУ методом Гаусса-Зейделя

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#define MIN 0.0001
#define ITER 7
#define CONST 3

void fillMatrix(FILE* f, double** m, int sz1, int sz2); 
// Функция, считывающая значения коэффициентов из файла и заполяющая ими матрицу
void print2(double** m, int sz1); // Функция печати двумерного массива
void print1(double* x, int sz1); // Функция печати одномерного массива
double* create1(int sz); // Функция выделения памяти под массив решения системы
double** create2(int sz1, int sz2); // Функция, создающая двумерный массив
void erase2(double** m, int sz1, int sz2); // Функция, удаляющая из памяти двумерный массив
int check_diagonal(double** m, int* index, int sz2); // Функция, проверяющая, есть ли на диагонали нули
double* get_sum(double** m, int sz1); // Функция, считающая сумму эл-тов каждой строки
int check_convergence(double** m, int* index, double* sum, int sz1); // Функция, проверяющая ДУС
double* iteration(double** m, double* x, double* pred, int sz1); // Итерация
double find_max_diff(double* x, double* pred, int sz1); // Поиск максимальной разности значений после итерации
void copy_int(int* t, int* x, int sz1); // Функция, копирующая массив индексов
void copy_double(double* pred, double* x, int sz1);
// Копирование значений предыдущей итерации для последующего сравнения с новыми значениями
int get_solution(double** m, double* x, int sz1); // Поиск решения системы
double** swap_lines(double** m, int sz1, int* index); // Функция для перестановки строк в соответствии с массивом индексов
void swap(int* a, int* b); // Функция для перестановки в массиве индексов
int* create_int(int sz1); // Создание и инициализация массива индексов
int solve_system(double** m, double* sum, double* x, double* pred, int* index, int sz1); 
// Функция, решающая систему без нулей на диагонали
int solve_system_with_zero(double** m, double* sum, double* x, double* pred, int* index, int sz1);
// Функцция, решающая систему с нулями на диагонали
int solve_divergent_system(double** m, double* x, double* pred, int sz1); // Решение системы, для которой не выполнилось ДУС
int solve_convergent_system(double** m, double* x, double* pred, int sz1); // Решение системы, для которой выполнилось ДУС

void copy_double(double* pred, double* x, int sz1)
{ // Копирование значений предыдущей итерации для последующего сравнения с новыми значениями
	int i;
	for (i=0; i<sz1; i++)
		pred[i]=x[i];
}

int solve_divergent_system(double** m, double* x, double* pred, int sz1) // Решение системы, для которой не выполнилось ДУС
{
	int iter=0, count=0;
	double max=0, max1=max;
	while (iter<ITER)
	{ // Пока не прошло заданное кол-во итераций
		x=iteration(m, x, pred, sz1); 
		print1(x, sz1);
		max=find_max_diff(x, pred, sz1); // Ищем макс. разницу на соседних шагах
		printf("%15E\n", max);
		if (max<max1)
			count++; // Если максимальная разность уменьшается, увеличиваем счётчик
		else
			count=0;
		max1=max;
		iter++;
	}
	if (count<CONST) 
	{ // Если разность не уменьшилась заданное кол-во раз подряд, значит, система не может быть решена
		free(x); 
		x=0; 
	}
	return (x==0) ? 0 : 1;	
}

int solve_convergent_system(double** m, double* x, double* pred, int sz1) // Решение системы, для которой выполнилось ДУС
{
	int b=0;
	while (b==0) // Пока не найдём решение системы, продолжаем выполнять итерации
	{
		x=iteration(m, x, pred, sz1); 
		print1(x, sz1);
		if (find_max_diff(x, pred, sz1)<MIN)
			return b=1; // Если максимальная разность м/у текущими и предыдущими значениями уже очень мала, выходим из функции
	}
}

int* create_int(int sz1) // Создание и инициализация массива индексов
{
	int i;
	int* mas=(int*) malloc (sz1*sizeof(int));
	for (i=0; i<sz1; i++)
		mas[i]=i;
	return mas;
}

void swap(int* a, int* b){ // Функция для перестановки значений в массиве индексов
    int temp=*a;
    *a=*b;
    *b=temp;
}

double** swap_lines(double** m, int sz1, int* index){ // Функция для перестановки строк в соответствии с массивом индексов
	double** m1=(double**) malloc (sz1*sizeof(double*)); // Создаём массив указателей
	int i;	
	for (i=0; i<sz1; i++)
		m1[i]=m[index[i]]; // Строки меняются местами
	free(m);
	return m1;	
}

void copy_int(int* t, int* x, int sz1)
{ // Функция, копирующая массив индексов
	int i;
	for (i=0; i<sz1; i++)
		t[i]=x[i];
}

double find_max_diff(double* x, double* pred, int sz1) // Поиск максимальной разности значений после итерации
{
	double max=0;
	int i;
	for (i=0; i<sz1; i++)
		if(fabs(x[i]-pred[i])>max)
    		max=fabs(x[i]-pred[i]);
    return max;
}

int solve_system(double** m, double* sum, double* x, double* pred, int* index, int sz1)
{  // Функция, решающая систему без нулей на диагонали
	if (check_convergence(m, index, sum, sz1)) // Проверка условия сходимости
		return solve_convergent_system(m, x, pred, sz1); // Если ДУС выполняется, просто решаем систему
	else
		return solve_divergent_system(m, x, pred, sz1); // Если не выполняется, решаем с контролем итераций
}

int transpose_lines(double** m, int* index, int* mix_index, double* sum, int num, int sz1)
{ // Функция для поиска удачной перестановки строк
    int b=0;
	if (num==sz1-1)
	{ // Если перестановка уже выполнилась: 
		if (check_diagonal(m, index, sz1))
		{ 
			if (check_convergence(m, index, sum, sz1)) // Если выполнилось ДУС, значит, перестановка удачная
			{
				printf("Permutation was found successfully\n");
				b=1; // Перестановка удачная
			}
			else // Если ДУС не выполнилось, то запоминаем данную перестановку для случая, если удачная так и не найдётся
				copy_int(mix_index, index, sz1); 
		}
		return b;
   	}
	int i;
    for(i=num; i<sz1; i++)
	{
	   	swap(&index[num], &index[i]);
	    b=transpose_lines(m, index, mix_index, sum, num+1, sz1);
	    if (b==1) return b;
	    swap(&index[num], &index[i]);    
    }
}

int solve_system_with_zero(double** m, double* sum, double* x, double* pred, int* index, int sz1)
{ // Функция, решающая систему с нулями на диагонали
	int* mix_index=create_int(sz1); // Массив индексов для подходящей подстановки
	if (transpose_lines(m, index, mix_index, sum, 0, sz1)==1)
	{ // Если найдена нужная перестановка с выполнением ДУС, то переставляем строки и решаем систему
		m=swap_lines(m, sz1, index);
		print2(m, sz1);
		return solve_convergent_system(m, x, pred, sz1);
	}
	// Если подходящая перестановка с выполнением ДУС не нашлась, то смотрим нули на диагонали
	printf("No permutation found'\n");
	if (check_diagonal(m, mix_index, sz1))
	{ // Если удалось избавиться от нулей на диагонали
		m=swap_lines(m, sz1, mix_index); // Меняем строки в соответствии с перестановкой
		print2(m, sz1);
		free(mix_index);
		return solve_divergent_system(m, x, pred, sz1); // Решаем систему с контролем итераций
	} 
	else // Иначе система не может быть решена этим методом
		return 0;
}

int get_solution(double** m, double* x, int sz1) // Поиск решения системы
{
	int b;
	double* pred=create1(sz1); // Создаётся массив, в котором будут храниться корни предыдущей итерации
	double* sum=create1(sz1); // Создаётся массив сумм
	sum=get_sum(m, sz1);
	int* index=create_int(sz1); // Создаётся массив индексов
	if (check_diagonal(m, index, sz1)) // Если нет нулей на диагонали, вызываем функцию для решения системы без нулей
		b=solve_system(m, sum, x, pred, index, sz1);
	else  // Если на диагонали есть нули, вызываем функцию для решения системы с нулями на диагонали
		b=solve_system_with_zero(m, sum, x, pred, index, sz1);
	free(pred);
	free(sum);
	free(index);
	return b;
}

double* iteration(double** m, double* x, double* pred, int sz1) // Итерация
{
	int i, j;
	double n;
	copy_double(pred, x, sz1); // Записываем значения для последующего сравнения
	for (i=0; i<sz1; i++) 
	{ // Для каждого уравнения системы
		x[i]=1/m[i][i]; // Записали коэффициент
		n=m[i][sz1]; // Записали свободный член
		for (j=0; j<sz1; j++)
		{
			if (i!=j)
				n-=m[i][j]*x[j]; // Вычитаем все остальные неизвестные
		}
		x[i]*=n; 
	}
	return x;
}

double* get_sum(double** m, int sz1) // Функция, считающая сумму эл-тов каждой строки
{
	double* s=create1(sz1);
	int i, j;
	for (i=0; i<sz1; i++)
	{
		for (j=0; j<sz1; j++)
			s[i]+=fabs(m[i][j]); // Считается сумма эл-тов строки
	}
	return s;
}

int check_convergence(double** m, int* index, double* sum, int sz1) // Функция, проверяющая ДУС
{
	int i, count=0;
	double d;
	for (i=0; i<sz1; i++)
	{
		d=sum[index[i]]-fabs(m[index[i]][i]); // Из массива сумм вычитаем диагональный элемент
		if ((fabs(m[index[i]][i]) >= d)) // Если абс. величина диаг. элемента >= абс. величины суммы оставшихся эл-тов
		{
			if ((fabs(m[index[i]][i]) > d)) // проверяем, чтобы был хотя бы 1 эл-т, для которого выполнится строгое неравенство
				count++; // Считаем, чтобы понять, имеется ли такой элемент
		}	 
		else		
			return 0; // Если хотя бы для одной строки условие не выполняется, то не выполняется и ДУС
	}
	return (count > 0);
}

int check_diagonal(double** m, int* index, int sz1) // Функция, проверяющая, есть ли на диагонали нули
{
	int i;
	for (i=0; i<sz1; i++)
		if (m[index[i]][i]==0)
			return 0; // Если хотя бы один диагональный элемент равен нулю, сразу выходим из функции со значением 0
	return 1; 
}

void fillMatrix(FILE* f, double** m, int sz1, int sz2)
{ // Функция, считывающая значения коэффициентов из файла и заполняющая ими матрицу
	int i, j;
	for (i=0; i<sz1; i++)
		for (j=0; j<sz2; j++)
			fscanf(f, "%lf", &m[i][j]);
}

double* create1(int sz) // Функция выделения памяти под массив решения системы
{
	double* m=(double*) malloc (sz*sizeof (double));
	int i;
	for (i=0; i<sz; i++)
		m[i]=0;
	return m;
}

double** create2(int sz1, int sz2) // Функция, создающая двумерный массив
{ 
	double** m=(double**) malloc (sz1*sizeof (double*)); // Выделяется память под массив указателей
	int i;
	for (i=0; i<sz1; i++) 
		m[i]=(double*) malloc (sz2*sizeof (double)); // Выделяется память под массив элементов типа double
	return m;
}

void erase2(double** m, int sz1, int sz2) // Функция, удаляющая из памяти двумерный массив
{ 
	int i;
	for (i=0; i<sz1; i++) 
		free(m[i]); // Очищаются массивы элементов
	free(m); // Очищается массив указателей
}

void print1(double* x, int sz1) // Функция печати одномерного массива решений
{
	int i;
	for (i=0; i<sz1; i++)
		printf("x%d=%15E ", i+1, x[i]);
	printf("\n");
}

void print2(double** m, int sz1) // Функция печати двумерного массива
{ 
	int i, j;
	for (i=0; i<sz1; i++){
		for (j=0; j<sz1+1; j++) 
			printf ("%15E", m[i][j]); 
		printf ("\n");
	}
	printf("\n");
}

int main()
{
	int sz1, sz2, b;
	FILE* f=fopen("SLE for Gauss-Seidel method.txt", "r"); // Открываем файл для чтения
	assert(f!=NULL); // Проверяем, открылся ли файл
	fscanf(f, "%i %i", &sz1, &sz2); // Считываем размеры матрицы, хранящиеся в 1ой строке файла
	double** m=create2(sz1, sz2); // Выделяем память под матрицу
	double* x=create1(sz1); // Выделяем память под массив решений
	fillMatrix(f, m, sz1, sz2); // Заполняем матрицу значениями из файла
	fclose(f);
	print2(m, sz1);
	if (get_solution(m, x, sz1))
		print1(x, sz1);
	else 
		printf("System can't be solved using this method");
	free(x);
	erase2(m, sz1, sz2);
	return 0;
}
