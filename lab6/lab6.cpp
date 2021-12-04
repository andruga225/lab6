// lab6.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

using namespace std;

int variant = 4;
int x0 = 1;
int xn = 2;
const double eps = 10e-8;

double f(double x)
{
	if (variant == 4)
		return 2 * exp(x) - 5 * x;
	return exp(x) + x + 1;
}

void simpsonMethod()
{
	double sum = 0;
	double runge = 1 / 15.;
	double n = 1, count = 0;
	double error = 1;
	double dx = xn - x0;

	double sumPrev = (xn - x0) * (f(x0) + 4 * f((xn + x0) / 2) + f(xn)) / 6.;
	double sumNext = 0;

	double h = dx / n;

	for (int i = 0; i < n*2; ++i)
	{
		double _x0 = x0 + i*h;
		double _xn = x0 + (i + 1) * h;

		sum += (_xn - _x0) * (f(_x0) + 4 * f((_x0 + _xn) / 2) + f(_xn)) / 6.;
	}

	while (abs(error) > eps) {
		h = (dx) / n;
		double pred_sum = sum;

		sum = 0;

		for (int i = 0; i < n; ++i)
		{
			double _x0 = x0 + i * h;
			double _xn = x0 + (i + 1) * h;

			sum += (_xn - _x0) * (f(_x0) + 4 * f((_x0 + _xn) / 2) + f(_xn)) / 6.;
		}
		count += n;

		error = (sum - pred_sum) * runge;
		double k = log((sum - sumPrev) / (sumNext - sumPrev) - 1) / log(0.5);
		sumPrev = sumNext;
		sumNext = sum;

		cout << sum<<endl;
		n *= 2;
	}

}

void threePointGauss()
{
	double a0 = 5 / 9., a1 = 8 / 9.;
	double sum = 0;
	double runge = 1 / 63.;
	double n = 1, count = 0;
	double error = 1;
	double dx = xn - x0;

	double sumPrev = 0;
	double sumNext = 0;

	double h = dx / n;

	for(int i=0;i<n;++i)
	{
		double x = x0 + (2 * i + 1) * h / 2;

		sumPrev += (h / 2) * (a0 * f(x - h * sqrt(0.6) / 2) + a1 * f(x) + a0 * f(x + h * sqrt(0.6) / 2));
	}

	for (int i = 0; i < 2*n; ++i)
	{
		double x = x0 + (2 * i + 1) * h / 2;

		sumNext += (h / 2) * (a0 * f(x - h * sqrt(0.6) / 2) + a1 * f(x) + a0 * f(x + h * sqrt(0.6) / 2));
	}

	while(abs(error)>eps||n<5)
	{
		h = dx / n;

		double predSum = sum;

		sum = 0;
		for (int i = 0; i < n; ++i)
		{
			double x = x0 + (2 * i + 1) * h / 2;

			sum += (h / 2) * (a0 * f(x - h * sqrt(0.6) / 2) + a1 * f(x) + a0 * f(x + h * sqrt(0.6) / 2));

			count += 3;
		}

		error = (sum - predSum) * runge;

		double k = log((sum - sumPrev) / (sumNext - sumPrev) - 1) / log(0.5);

		sumPrev = sumNext;
		sumNext = sum;

		n *= 2;

		cout << sum << endl;
	}
}

int main()
{
	simpsonMethod();
	threePointGauss();
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
