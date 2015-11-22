//#pragma once
#include <time.h>

double Abs(double t)
{
    if (t < 0) return -t;
    else return t;
}


struct Matr // ����� ��� �������
{
	double** M; // ���������� ���� �������
	int row;    // ���������� �����
	int col;    // ���������� ��������
};

// ������� ����������� �������. �� �� ��������. ���������� ������ �������
Matr InitMatr(int n, int m)
{
	Matr temp;
	temp.M = new double*[n];
	temp.row = n;
	temp.col = m;
	for (int i = 0; i < temp.row; i++)
	{
		temp.M[i] = new double[m];
	}

	return temp;
}
void DeleteMatr(Matr m)
{
	
	for (int i = 0; i < m.row; i++)
	{
		delete[] m.M[i];
	}
	delete[] m.M;
}

// �������� ������� �� ������
void Print(Matr A)
{
	printf("\n\nEnter Matrix To Console!!!\n");
	for (int i = 0; i < A.row; i++)
	{
		for (int j = 0; j < A.col; j++)
			printf("%3.8f    ", A.M[i][j]);
		printf("\n");
	}
}

// ���������� ������� ������
Matr EnterZero(Matr A)
{
	for (int i = 0; i < A.row; i++)
	{
		for (int j = 0; j < A.col; j++)
			A.M[i][j] = 0;
	}
	return A;
}

//���������������� �������
void Transponize(Matr A,Matr T)
{ 
	for (int i = 0; i < A.row; i++)
	{
		for (int j = 0; j < A.col; j++)
			T.M[i][j] = A.M[j][i];
	}
}

// ��������� ������� ���������� ������������� �������
Matr EnterRandom(Matr A)
{
	srand((unsigned int)time(0));
	for (int i = 0; i < A.row; i++)
	{
		for (int j = 0; j < A.col; j++)
			A.M[i][j] = (double)rand()/(double)RAND_MAX;
	}
	return A;
}

// �������� ��������� �������
Matr EnterUnit(Matr A)
{
	for (int i = 0; i < A.row; i++)
	{
		for (int j = 0; j < A.col; j++)
			A.M[i][j] = 0;
		A.M[i][i] = 1;
	}
	return A;
}

// ������ ������� � ����
void MatrToFile(char *Name, Matr A, int steprow = 1, int stepcol = 1)
{
	FILE *fp;
	printf("Name file It Is ________ %s\n", Name);

	if ((fp = fopen(Name,"w"))==NULL)
	{
		printf("������ ��� �������� �����.\n");
		exit(1);
	}
	for (int i = A.row - 1; i >= 0; i= i - steprow)
	{
		for (int j = 0; j < A.col; j = j + stepcol)
			fprintf(fp, "%7.6f;", A.M[i][j]);
		fprintf(fp, "\n");
	}

	fclose(fp);
}

// ���������� �������� ���� ������
Matr operator -(Matr A, Matr B)
{
	for (int i = 0; i < A.row; i++)
		for (int j = 0; j < A.col; j++)
			A.M[i][j] -= B.M[i][j];

	return A;
}

// ���������� ����� ���� ������
Matr operator +(Matr A, Matr B)
{
	for (int i = 0; i < A.row; i++)
		for (int j = 0; j < A.col; j++)
			A.M[i][j] += B.M[i][j];

	return A;
}

// ��������� ���� ������
Matr operator *(Matr A, Matr B)
{
	Matr C;
	C = InitMatr(A.row, B.col);
	C = EnterZero(C);

	for (int i = 0; i < A.row; i++)
    	for (int j = 0; j < B.col; j++)
			for (int l = 0; l < A.col; l++ )
			C.M[i][j] += A.M[i][l] * B.M[l][j];

	return C;
}

// ��������� ������� �� ������
Matr operator * (double a, Matr B)
{
	Matr C;
	C = InitMatr(B.row, B.col);
	C = EnterZero(C);

	for (int i = 0; i < B.row; i++)
    	for (int j = 0; j < B.col; j++)
			C.M[i][j] = B.M[i][j] * a;

	return C;
}

// ������� ����� �������
Matr Minor(Matr A, int row, int col)
{
	Matr T;
	T = InitMatr(A.row - 1, A.col - 1);
	T = EnterZero(T);

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
			T.M[i][j] = A.M[i][j];
		for (int j = col; j < A.col - 1; j++)
			T.M[i][j] = A.M[i][j + 1];
	}

	for (int i = row; i < A.row - 1; i++)
	{
		for (int j = 0; j < col; j++)
			T.M[i][j] = A.M[i + 1][j];
		for (int j = col; j < A.col - 1; j++)
			T.M[i][j] = A.M[i + 1][j + 1];
	}

	return T;
}

// ������� ������������ �������
double Determ(Matr A)
{
	double det = 0;

	if (A.row == 2)
	{
		det = A.M[0][0] * A.M[1][1] - A.M[0][1] * A.M[1][0];
		return det;
	}
	else
	{

		for (int i = 0; i < A.row - 1; i += 2)
		{

			Matr m = (Minor(A, 0, i));
			Matr m1 = Minor(A, 0, i + 1);
			det += A.M[0][i] * Determ(m) - A.M[0][i + 1] * Determ(m1);

			DeleteMatr(m);
			DeleteMatr(m1);

		}
		if (A.row%2 == 1)
		{
			Matr m = Minor(A, 0, A.row - 1);
			det += A.M[0][A.row - 1] * Determ(m);
			DeleteMatr(m);
		}

		return det;
	}
}


// �������� �������
Matr InversMatr(Matr A)
{

	Matr T;
	T = InitMatr(A.row, A.col);
	T = EnterUnit(T);

	if (Determ(A) == 0)
	{
		printf("Inverse Matrix not Found!!!!\n");
		return T;
	}

	// ���� �������� ������� ������� ����������
	printf("\n Forward________ \n");

	for (int i = 0; i < A.row - 1; i++)
	{
		double max = Abs(A.M[i][i]);
		int number = i;
		// ��� ������� ������� ���� ������� �������

		for (int j = i + 1; j < A.row; j++)
			if (Abs(A.M[j][i]) > max)
			{
				max = Abs(A.M[j][i]);
				number = j;
			}

		// ������ ��� ������ ������� �� ������� i
		if (number > i)
		{

			// ������ ����� � ��������
			for (int j = i; j < A.col; j++)
			{
				double temp;
				//������ ����� � �������� �������
				temp = A.M[i][j];
				A.M[i][j] = A.M[number][j];
				A.M[number][j] = temp;

			}

			for (int j = 0; j < A.col; j++)
			{
				double temp;
				//������ ����� � ������� ����������
				temp = T.M[i][j];
				T.M[i][j] = T.M[number][j];
				T.M[number][j] = temp;
			}
		}

		// ���������, ���� �� �������� �������...
		if (A.M[i][i] == 0)
		{
			printf("Inverse Matrix not Found!!!!\n");
			return T;
		}

		// ������ ����������...
		for (int j = i+1; j < A.row; j++)
		{
			double temp = A.M[j][i]/A.M[i][i];
			for(int l = i; l < A.col; l++)
			{
				A.M[j][l] = A.M[j][l] - temp * A.M[i][l];
			}
			for (int l = 0; l < A.col; l++)
			{
				T.M[j][l] = T.M[j][l] - temp * T.M[i][l];
			}

		}
	}

	// �������� ����������
	printf("\n\n Inverse ______________\n");
	for (int i = A.row - 1; i >= 0; i --)
	{
		double temp = A.M[i][i];
		for (int j = 0; j < A.col; j++)
		{
			T.M[i][j] = T.M[i][j] / temp;
		}
		for (int j = i; j < A.col; j++)
		{
			A.M[i][j] = A.M[i][j] / temp;
		}

		for (int l = T.row - 1; l > i; l--)
		{
			for (int j = 0; j < A.col; j++)
			{
				T.M[i][j] = T.M[i][j] - A.M[i][l] * T.M[l][j];
			}
		}
		for (int j = i + 1; j < A.row; j++)
		{
			A.M[i][j] = 0;
		}
	}

	return T;
}
