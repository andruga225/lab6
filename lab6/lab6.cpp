// lab6.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include<iomanip>
#include <fstream>

using namespace std;

int variant = 11;
int x0 = 1;
int xn = 2;
const double eps = 10e-8;

ofstream fout("output.txt");

double f(double x)
{
	if (variant == 4)
		return  2 * exp(x) - 5 * x;
	return exp(x) + x + 1;
}

double F()
{
	if (variant == 4)
		return 1.841548540943211;// 2*exp(2.0)-(15.0/2.0)-2*exp(1.0);
	return 7.170774270471606;//exp(2.0) + (5.0 / 2.0) - exp(1.0);
}

double df(double x)
{
	if (variant == 4)
		return 2 * exp(x) - 5;
	return exp(x) + 1;
}

void simpsonMethod()
{
	fout << "Формула Симпсона\n";

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

		sumNext += (_xn - _x0) * (f(_x0) + 4 * f((_x0 + _xn) / 2) + f(_xn)) / 6.;
	}

	fout << "N | h | Interral | Погрешность | Оценка погр. | k\n";

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
		double pogr = sum - pred_sum;
		double k = log((sum - sumPrev) / (sumNext - sumPrev) - 1) / log(0.5);
		sumPrev = sumNext;
		sumNext = sum;

		if (n == 1)
		{
			fout << fixed << setw(5) << setprecision(0) << n << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << sum << "|"<<scientific<<abs(pogr)<<'\n';
		}
		else if (n == 2)
		{
			fout << fixed << setw(5) << setprecision(0) << n << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << sum << "|" << setw(15) << setprecision(12) << scientific<<abs(pogr)<<"|" << error << "|\n";
		}
		else
		{
			fout << fixed << setw(5) << setprecision(0) << n << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << sum << "|" << setw(15) << setprecision(12) << scientific << abs(pogr) << "|" << error << "|" << fixed << setw(8) << setprecision(5) << k << '\n';
		}

		n *= 2;
	}

	fout << "Result: " << fixed<<setprecision(15)<<sum<<'\n';
	fout << "Kobr: " << fixed << setprecision(0) << count << '\n';
}

void threePointGauss()
{
	fout << "Формула Гаусса-3\n";

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

	fout << "N | h | Interral | Погрешность | Оценка погр. | k\n";

	while(abs(error)>eps||n<5)
	{
		h = dx / n;

		double predSum = sum;

		sum = 0;
		for (int i = 0; i < n; ++i)
		{
			double x = x0 + (2 * i + 1) * h / 2.;

			sum += (h / 2) * (a0 * f(x - h * sqrt(0.6) / 2) + a1 * f(x) + a0 * f(x + h * sqrt(0.6) / 2));

			count += 3;
		}

		error = (sum - predSum) * runge;
		double pogr = sum - predSum;
		double k = log((sum - sumPrev) / (sumNext - sumPrev) - 1) / log(0.5);

		sumPrev = sumNext;
		sumNext = sum;

		if (n == 1)
		{
			fout << fixed << setw(5) << setprecision(0) << n << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << sum << "|" << scientific << abs(pogr) << '\n';
		}
		else if (n == 2)
		{
			fout << fixed << setw(5) << setprecision(0) << n << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << sum << "|" << setw(15) << setprecision(12) << scientific << abs(pogr) << "|" << error << "|\n";
		}
		else
		{
			fout << fixed << setw(5) << setprecision(0) << n << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << sum << "|" << setw(15) << setprecision(12) << scientific << abs(pogr) << "|" << error << "|" << fixed << setw(8) << setprecision(5) << k << '\n';
		}
		n *= 2;
	}
	fout << "Result: " << fixed << setprecision(15) << sum << '\n';
	fout <<"Kobr: " << fixed << setprecision(0) << count << '\n';
}



void first_trapeze_metod()
{
	fout << "Метод трапеций" << endl;
	int N = 1;
	double h = (xn - x0) / N;
	double s1, s2=0, s3, sk;
	double k;
	double ocenka_pogr, pogr;
	double alpha = 0.5;
	sk = (f(x0) + f(xn))/2;
	double count = 2;
	s1 = sk * h;
	fout << "N | h | Interral | Погрешность | Оценка погр. | k\n";

	fout <<fixed<<setw(5)<< N << "|"  <<setw(6)<<setprecision(4)<< h << "|"  << setw(12)<<setprecision(9) << s1 << endl;
	do
	{
		N *= 2;
		h /= 2;
		for (int i = 1; i <= N - 1; i += 2)
		{
			sk += f(x0 + h * i);
			count++;
		}
		s3 = s2;
		s2 = s1;
		s1 =sk*h;

		ocenka_pogr = pow(alpha, 2) / (1 - pow(alpha, 2)) * (s1 - s2);
		pogr = F()-s1;
		if (N == 2) {
			fout <<fixed<<setw(5)<<setprecision(0) << N << "|"  << setw(6)<<setprecision(4) << h << "|"  << setw(12)<<setprecision(9) << s1 << "|" << scientific << setw(15)<<setprecision(12) << abs(pogr) <<"|" << ocenka_pogr << endl;
		}
		else
		{
			k = 1 / log(alpha) * log((s1 - s3) / (s2 - s3) - 1);
			fout << fixed << setw(5) << setprecision(0) << N << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << s1 << "|" << scientific << setw(15) << setprecision(12) << abs(pogr) << "|" << ocenka_pogr<<"|" << fixed << setprecision(4) << k << endl;
		}
	} while (abs(s1-s2)/3 > eps); // тут вроде такая константа

	fout << "Результат: " <<fixed<<setprecision(15)<< s1 << endl;
	fout << "Кол-во обращений: " <<fixed<<setprecision(0)<<count<<endl;
}

void second_trapeze_metod()
{
	fout << "Метод трапеций (модифицированной с помощью сплайна)" << endl;
	int N = 1;
	double h = (xn - x0) / N;
	double s1, s2 = 0, s3, sk;
	double k;
	double ocenka_pogr, pogr;
	double alpha = 0.5;

	double spr = df(x0) - df(xn);
	sk = (f(x0) + f(xn)) / 2;
	double count = 2;
	s1 = sk * h + h*h/12*spr;

	fout << "N | h | Interral | Погрешность | Оценка погр. | k\n";

	fout << fixed << setw(5) << N << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << s1 << endl;
	do
	{
		N *= 2;
		h /= 2;
		for (int i = 1; i <= N - 1; i += 2)
		{
			sk += f(x0 + h * i);
			count++;
		}
		s3 = s2;
		s2 = s1;
		s1 = sk * h+spr*h*h/12;

		ocenka_pogr = pow(alpha, 4) / (1 - pow(alpha, 4)) * (s1 - s2);
		pogr = F()-s1;
		if (N == 2) {
			fout << fixed << setw(5) << setprecision(0) << N << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << s1 << "|" << scientific << setw(15) << setprecision(12) << abs(pogr) << "|" << ocenka_pogr << endl;
		}
		else
		{
			k = 1 / log(alpha) * log((s1 - s3) / (s2 - s3) - 1);
			fout << fixed << setw(5) << setprecision(0) << N << "|" << setw(6) << setprecision(4) << h << "|" << setw(12) << setprecision(9) << s1 << "|" << scientific << setw(15) << setprecision(12) << abs(pogr) << "|" << ocenka_pogr << "|" << fixed << setprecision(4) << k << endl;
		}
	} while (abs(s1 - s2) / 3 > eps); //тут я не уверен насчёт константы

	fout << "Результат: " << fixed << setprecision(15) << s1 << endl;
	fout << "Кол-во обращений: " << fixed << setprecision(0) << count << endl;
}

int main()
{
	setlocale(LC_ALL, "rus");

	/*резюмирую вроде всё сходится с вариантом с трекусов но нужно сделать вывод покрасивше*/
	fout <<fixed<<setprecision(15)<< "J = " << F()<<endl;
	first_trapeze_metod();
	second_trapeze_metod();
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
