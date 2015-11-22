#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include "Matrix.h"
#include "Vect.h"


using namespace std;

int a = 0;
int b = 1;
double eps = 0.0005;

Vect operator *(Matr A, Vect b)
{
	Vect c;
	c = InitVect(b.size);
	c = EnterZero(c);

	for (int i = 0; i < b.size; i++)
	{
		for (int j = 0; j < b.size; j++)
			c.V[i] += A.M[i][j] * b.V[j];
	}

	return c;
}


double y(double t)
{
    return t + 7.0 / 3;
}

double K(double t, double s, double x)
{
    return (t + x * x + 2) * (t + x * x + 2) / (t + s * s + 2);
}


Vect F(Vect t, Vect s, Vect x, Vect y, double h)
{
    Vect T = InitVect(t.size);
    T = EnterZero(T);

    for (int i = 0; i < t.size; i++)
    {
        double temp = K(t.V[i], s.V[0], x.V[0]) * h / 2;
        for (int j = 1; j < t.size - 1; j++)
            temp += h * K(t.V[i], s.V[j], x.V[j]);
        temp += K(t.V[i], s.V[s.size - 1], x.V[x.size - 1]) * h / 2;
        T.V[i] = temp - y.V[i];
        //cout << T.V[i] << "!!!!!!!!\n";
    }
    return T;
}

double dK_dx(double t, double s, double x)
{
    return 4 * (t + x * x + 2) * x / (t + s * s + 2);
}

void Jacobi(Vect t, Vect s, Vect x, Vect y, double h,Matr J)
{
    for (int i = 0; i < t.size; i++)
    {
        J.M[i][0] = dK_dx(t.V[i], s.V[0], x.V[0]) * h / 2;
        for (int j = 1; j < t.size - 1; j++)
            J.M[i][j] = h * dK_dx(t.V[i], s.V[j], x.V[j]);
        J.M[i][t.size - 1] = dK_dx(t.V[i], s.V[s.size - 1], x.V[x.size - 1]) * h / 2;
    }
}


int main()
{
    cout << "Hello world!" << endl;
    int N = 6;
    double h = ((double)(b - a)) / N;
    double lambda = 0.01;
    N++;
    Vect t = InitVect(N);
    Vect s = InitVect(N);
    Vect f = InitVect(N);
    Vect x = InitVect(N);

    for (int i = 0; i < N; i++)
    {
        t.V[i] = i * h;
        s.V[i] = i * h;
        f.V[i] = y(i * h);
        x.V[i] = y(i * h);
    }
    int step = 0;
    Vect dx = InitVect(N);
    double nev = 0;
	Matr JT = InitMatr(N, N);
	Matr J = InitMatr(N, N);
    Matr I = InitMatr(N, N);
    I = EnterUnit(I);
	auto lambdaxI = lambda * I;

    do
    {
        //Print(x);
		Jacobi(t, s, x, f, h, J);
       // Print(J);
		Transponize(J, JT);
        //Print(JT);
        //Print(JT * J);
        dx = InversMatr(JT * J + lambdaxI) * (JT * F(t, s, x, f, h));
        //Print(JT);

			cout << "_______\n\n";

        //Print(dx);
        x = x - dx;
        nev = Norma(F(t, s, x, f, h));

			cout << "!!!! Nevyazka - " << nev << "\n\n";

        step++;
        //_getch();
    } while (nev > eps);
    Print(x);
    cout << "Count of Iteration --- " << step;
    return 0;
}
