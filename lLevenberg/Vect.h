
// Описание вектор-столбца
struct Vect
{
	int size;
	double * V;
};

// Создание вектор-столбца
Vect InitVect(int n)
{
	Vect T;
	T.size = n;
	T.V = new double[n];
	return T;
}

// Вывод на экран
void Print(Vect T)
{
	printf("\nPrint Vector To Console!!!!\n");
	for (int i = 0; i < T.size; i++)
		printf("%30.28f\n", T.V[i]);
}

// Заоплнение нулями
Vect EnterZero(Vect T)
{
	for (int i = 0; i < T.size; i++)
		T.V[i] = 0;

	return T;
}


// Заполнение случайными числами
Vect EnterRandom(Vect T)
{
	srand((unsigned int)time(0));
	for (int i = 0; i < T.size; i++)
		T.V[i] = (double)rand()/(double)RAND_MAX;

	return T;
}

// Распечатывает два веткора для сравнения
void Print(Vect a, Vect b)
{
	printf("\nPrint Two Vectors To Console!!!!\n");
		for (int i = 0; i < a.size; i++)
			printf("%10.8f  |  %10.8f\n", a.V[i], b.V[i]);
}

// считает норму вектора
double Norma(Vect a)
{
	double temp = 0;
	for (int i = 0; i < a.size; i++)
	{
		temp += pow(a.V[i], 2.0);
	}

	temp = pow(temp, 0.5);
	return temp;
}

// считает разность двух векторов
Vect operator -(Vect a, Vect b)
{
	for (int i = 0; i < a.size; i++)
	{
		a.V[i] -= b.V[i];
	}

	return a;
}

// считает сумму двух векторов
Vect operator +(Vect a, Vect b)
{
	for (int i = 0; i < a.size; i++)
	{
		a.V[i] += b.V[i];
	}

	return a;
}

// считает скалярное умножение двух векторов
double operator * (Vect a, Vect b)
{
	double temp = 0;
	for (int i = 0; i < a.size; i++)
		temp += a.V[i] * b.V[i];

	return temp;
}

// считает умножение вектора на скаляр
Vect operator * (double a, Vect b)
{
	for (int i = 0; i < b.size; i++)
		b.V[i] = a * b.V[i];

	return b;
}

// Запись веткора в файл
void VectToFile(char *Name, Vect b, int step = 1)
{
	FILE *fp;
	printf("Name file It Is ________ %s\n", Name);

	if ((fp = fopen(Name,"w"))==NULL)
	{
		printf("Ошибка при открытии файла.\n");
		exit(1);
	}
	for (int i = b.size - 1; i >= 0; i= i - step)
	{
		fprintf(fp, "%7.6f;\n", b.V[i]);
	}

	fclose(fp);
}
