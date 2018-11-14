
#include "stdafx.h"
#include "iostream"
#include "fstream"
#include "string"
#include "string.h"
#include "iomanip"
#include <sstream>
#include <ctime>
//#include "interpolation.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;
// N-количесетво строк в функции светимости, K-количество строк в изохроне 
//vis_mag_char, fi_char,mass_char,abs_mag_char - вид зв величина звезд скопления, количество звезд (значение функции светимости), масса, абсолютная зв величина - для считывания из файла
//константы: EBV - погл. r - расстояние до скопления ALPHA - доля двойных в скоплении, step-шаг с которым строится гистограмма для подсчета массы
//vis_mag: 0- видимая зв величина, начиная с которой рассматривается гиста, fin-предельная зв величина, 1,2 - вид зв величины для интерполяции, чтобы найти значение функции светимости на середине отрезка vis_mag0 и vis_mag0+step
//fi: 1,2 значение функции в точках vis_mag1 и vis_mag2, fi- интерполированное значение функции в середине отрезка vis_mag0 и vis_mag0+step
//abs_mag: 0 - абсолютная зв величина из vis_mag0, 1,2 - абс зв величины для интерполяции, чтобы найти значение массы для abs_mag0
//mass: 1,2 - масса в точках abs_mag1 и abs_mag2, mass- интерполированное значение массы в точке abs_mag0
//MASS_DOUBLE: _i- масса ТОЛЬКО ДВОЙНЫХ звезд из интервала vis_mag0 и vis_mag0+step, MASS_DOUBLE - масса ВСЕХ ДВОЙНЫХ звезд от начальной до конечной звездной величины
//MASS_SINGLE: _i- масса ТОЛЬКО ОДИНОЧНЫХ звезд из интервала vis_mag0 и vis_mag0+step, MASS_SINGLE - масса ВСЕХ ОДИНОЧНЫХ звезд от начальной до конечной звездной величины
//MASS_CLUSTER_NEW: масса скопления с приближением, что есть двойные
//MASS_CLUSTER_OLD: масса скопления с приближением, что все одиночные
//q-отношение массы главного компонента, к вторичному
double line_interp(double x1, double y1, double x2, double y2, double x)
{
	return y1 + (y2 - y1)*(x - x1) / (x2 - x1);
}
double lf(int N, int f, fstream& fin, double vis_mag0, double step, double ALPHA)
{
	char vis_mag_char[10], fi_char[10],a[10];
	double  vis_mag1, vis_mag2, fi1, fi2;
	double check;
	for (int i = 0; i <= N+1; i++) // нахождение значение функции светимости на середине интервала*длину интервала - количество звезд двойных
	{
		fin.seekg(i * 80);
		fin.read(a, 10);
		vis_mag1 = atof(a);

		fin.seekg(f, ios::cur);
		fin.read(fi_char, 10);
		fi1 = atof(fi_char);

		fin.seekg((i + 1) * 80);
		fin.read(vis_mag_char, 10);
		vis_mag2 = atof(vis_mag_char);

		fin.seekg(f, ios::cur);
		fin.read(fi_char, 10);
		fi2 = atof(fi_char);

		if ((vis_mag0 + step / 2) >= vis_mag1 && (vis_mag0 + step / 2) <= vis_mag2)
		{
			check= line_interp(vis_mag1, fi1, vis_mag2, fi2, (vis_mag0 + step / 2))*step*ALPHA;
			return check;
			break;
		}
	}
}
void to_mf(int K, fstream& fini, double abs_mag0, double fi, double& mass, double& MASS_DOUBLE)
{
	char mass_char[8], abs_mag_char[8];
	double mass1, mass2, abs_mag1, abs_mag2;
	for (int j = 0; j <= K - 2; j++)  // нахождение массы, соответствующей абс звездной величине серидины интервала
	{
		fini.seekg(j * 14);
		fini.read(mass_char, 5);
		mass1 = atof(mass_char);

		fini.seekg(1, ios::cur);
		fini.read(abs_mag_char, 6);
		abs_mag1 = atof(abs_mag_char);

		fini.seekg((j + 1) * 14);
		fini.read(mass_char, 5);
		mass2 = atof(mass_char);

		fini.seekg(1, ios::cur);
		fini.read(abs_mag_char, 6);
		abs_mag2 = atof(abs_mag_char);
		if (abs_mag0 <= abs_mag2 && abs_mag0 >= abs_mag1)//<> когда зв величины по возрастанию
		{
			mass = line_interp(abs_mag1, mass1, abs_mag2, mass2, abs_mag0);
			MASS_DOUBLE = MASS_DOUBLE + fi*mass;
			break;
		}
	}
}
double monte_karlo()
{
		char G_char[4], q_char[4];
		double q, F, G1, G2, G, q1, q2, qmax = 1, Fmax = 1, fi = 10,s;
		fstream finq;
		finq.open("func_gr.txt", ios::in);
		srand(time(NULL));
		for (int i = 0; i <= 1; i++)
		{
			s = rand();
			q = (rand() % 1000) / 1000.;
			F = (rand() % 1000) / 1000.;
			for (int j = 0; j <= 40; j++)
			{
				finq.seekg(j * 11);
				finq.read(q_char, 4);
				q1 = atof(q_char);

				finq.seekg(1, ios::cur);
				finq.read(G_char, 3);
				G1 = atof(G_char);

				finq.seekg((j + 1) * 11);
				finq.read(q_char, 4);
				q2 = atof(q_char);

				finq.seekg(1, ios::cur);
				finq.read(G_char, 3);
				G2 = atof(G_char);
				if (q <= q2 && q >= q1)
				{
					G = line_interp(q1, G1, q2, G2, q);
					break;
				}
			}
			if (F > G) i = i - 1;
			else return q;
		}
}
double lummass(double M)
{
	return pow(10,( -0.705*(log10(M))*(log10(M)) + 4.655*(log10(M)) - 0.025));
	//return pow(10, (4.040*(log10(M)) - 0.002));
}
int main()
{ 
	int N = 96, K = 58, Nsteps=30;//потом с клавиатуры
	double mod_dist = 10.2,EBV = 0.15, ALPHA = 0., step;// константы для скопления - потом с клавиатуры
	char vis_mag_char[10]; 
	double vis_mag0, vis_mag_fin, fi, abs_mag0, mass;
	double MASS_DOUBLE_i, MASS_DOUBLE = 0, MASS_SINGLE_i, MASS_SINGLE = 0, MASS_CLUSTER_NEW, MASS_CLUSTER_OLD_i, MASS_CLUSTER_OLD = 0, NEW_divided_by_OLD;
	double q, musor, mass_s, mass_d1, mass_d2, lf0, lf1, mass_d1_st,fx[2], ratio[31], mistake;
	fstream fin, fini, fout, fout2;
	fin.open("density_NGC7142.txt", ios::in);
	fini.open("mass-absmag.txt", ios::in); 
	fout.open("number_new.txt", ios_base::out);
	fout2.open("alpha_same_NGC7142.txt", ios_base::out);
	fout2.setf(ios::fixed);
	fout2.setf(ios::fixed);
	fout2 << "ALPHA" << '\t' << "MASS_CLUSTER_NEW" << '\t' << "MASS_CLUSTER_NEW_low" << '\t' << "MASS_CLUSTER_NEW_up" << '\t' << '\t' << "NEW_divided_by_OLD" << endl;
	for (ALPHA; ALPHA <= 0.9; ALPHA = ALPHA + 0.1)
	{
		cout <<endl<< ALPHA << endl;
		fout << endl << ALPHA << endl;
		ratio[0] = 0;
		for (int k=1; k<=30;k++)
		{ 
			cout << k;
			fin.seekg(0);
			fin.read(vis_mag_char, 10);
			vis_mag0 = atof(vis_mag_char); // начальная звездная величина

			fin.seekg((N - 1) * 80);
			fin.read(vis_mag_char, 10);
			vis_mag_fin = atof(vis_mag_char);// конечная звездная величина

			step = (vis_mag_fin - vis_mag0) / Nsteps;

			MASS_DOUBLE = 0;  MASS_SINGLE = 0; MASS_CLUSTER_OLD = 0,mistake=0;
			while (vis_mag0 +step/2<= vis_mag_fin)
		{
			//для только одиночных
			fi = round (lf(N, 7, fin, vis_mag0, step, 1));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную одиночной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_CLUSTER_OLD);
			
			//для доли двойных
			fi = round (lf(N,7, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную двойной
			to_mf(K, fini, abs_mag0, fi, mass_s, musor);
			lf0 = lummass(mass_s);
			for (int i=1; i<= fi ; i++)
			{
				//q = monte_karlo();
				q = 1;
				for (double i = 0.001; i <= 10; i = i + 0.001)
				{
					if (log10(lummass(0.001))*log(10) + log(1 + pow(M_E, log(10)*(-0.705*log10(q)*log10(q) + 4.655*log10(q) - 1.41*log10(0.001)*log10(q)))) - log(lf0) > 0) { mass_d1_st = 0.001; break;}
					mass_d1_st = i;
					fx[0] = log10(lummass(mass_d1_st))*log(10) + log(1 + pow(M_E, log(10)*(-0.705*log10(q)*log10(q) + 4.655*log10(q) - 1.41*log10(mass_d1_st)*log10(q)))) - log(lf0);
					mass_d1_st = i + 0.001;
					fx[1] = log10(lummass(mass_d1_st))*log(10) + log(1 + pow(M_E, log(10)*(-0.705*log10(q)*log10(q) + 4.655*log10(q) - 1.41*log10(mass_d1_st)*log10(q)))) - log(lf0);
					if (fx[0] * fx[1] < 0) break;
				}
				do 
				{
					mass_d1 = mass_d1_st;
					mass_d1_st = mass_d1 - (log10(lummass(mass_d1))*log(10) + log(1 + pow(M_E, log(10)*(-0.705*log10(q)*log10(q) + 4.655*log10(q) - 1.41*log10(mass_d1)*log10(q)))) - log(lf0)) / (log(10)*(4.655 - 1.41*log10(mass_d1)) + log(10)*(-1.41*log10(q))*(pow(M_E, log(10)*(-0.705*log10(q)*log10(q) + 4.655*log10(q) - 1.41*log10(mass_d1)*log10(q)))) / (1 + pow(M_E, log(10)*(-0.705*log10(q)*log10(q) + 4.655*log10(q) - 1.41*log10(mass_d1)*log10(q)))));
				} while (abs(mass_d1-mass_d1_st)>=0.0001);
				
				mass_d1 = mass_d1_st;
				mass_d2 = q*mass_d1;
				MASS_DOUBLE = MASS_DOUBLE + mass_d1+mass_d2;
				}
						
			//для доли одиночных 
			fi= round(lf(N, 7, fin, vis_mag0, step, 1))- round(lf(N, 7, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную одиночной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_SINGLE);
			fout << vis_mag0 << '\t' << MASS_DOUBLE << '\t' << MASS_SINGLE<<endl;
		vis_mag0 = vis_mag0 + step;
		} 

			MASS_CLUSTER_NEW = MASS_DOUBLE + MASS_SINGLE;
			ratio[k] = MASS_CLUSTER_NEW / MASS_CLUSTER_OLD;
			ratio[0] = ratio[0] + ratio[k];
			
		}
		NEW_divided_by_OLD = ratio[0] / 30.;
		for (int l=1; l<=30; l++)
		{
			mistake = mistake+(ratio[l]-NEW_divided_by_OLD)*(ratio[l] - NEW_divided_by_OLD)/29.;
		}
		if (mistake>=0 && mistake<=1) fout2 << ALPHA << '\t' << MASS_CLUSTER_NEW << '\t' << MASS_CLUSTER_OLD << '\t' << '\t' << NEW_divided_by_OLD << '\t' << pow(mistake,0.5)<< endl;
		else ALPHA = ALPHA - 0.1;
	}
	fin.close();
	fini.close();
	fout.close();
	//system("pause");
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
#include "stdafx.h"
#include "iostream"
#include "fstream"
#include "string"
#include "string.h"
#include "iomanip"
#include <sstream>
//#include "interpolation.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;
// N-количесетво строк в функции светимости, K-количество строк в изохроне, Nsteps- число разбиений на интервале вохможных звездных величин
//vis_mag_char - зв величина - для считывания из файла
//константы: EBV - погл. r - расстояние до скопления ALPHA - доля двойных в скоплении, step-шаг с которым строится гистограмма для подсчета массы
//vis_mag: 0- видимая зв величина, начиная с которой рассматривается гиста, fin-предельная зв величина
//fi: fi- интерполированное значение функции в середине отрезка vis_mag0 и vis_mag0+step
//abs_mag: 0 - абсолютная зв величина из vis_mag0
//mass: mass- интерполированное значение массы в точке abs_mag0
//MASS_DOUBLE: _i- масса ТОЛЬКО ДВОЙНЫХ звезд из интервала vis_mag0 и vis_mag0+step, MASS_DOUBLE - масса ВСЕХ ДВОЙНЫХ звезд от начальной до конечной звездной величины
//MASS_SINGLE: _i- масса ТОЛЬКО ОДИНОЧНЫХ звезд из интервала vis_mag0 и vis_mag0+step, MASS_SINGLE - масса ВСЕХ ОДИНОЧНЫХ звезд от начальной до конечной звездной величины
//MASS_CLUSTER_NEW: масса скопления с приближением, что есть двойные
//MASS_CLUSTER_OLD: масса скопления с приближением, что все одиночные
double line_interp(double x1, double y1, double x2, double y2, double x)
{
	return y1 + (y2 - y1)*(x - x1) / (x2 - x1);
}
double lf(int N, int f, fstream& fin, double vis_mag0, double step, double ALPHA)
{
	char vis_mag_char[10], fi_char[10],a[10];
	double  vis_mag1, vis_mag2, fi1, fi2;
	double check;
	for (int i = 0; i <= N+1; i++) // нахождение значение функции светимости на середине интервала*длину интервала - количество звезд двойных
	{
		fin.seekg(i * 80);
		fin.read(a, 10);
		vis_mag1 = atof(a);

		fin.seekg(f, ios::cur);
		fin.read(fi_char, 10);
		fi1 = atof(fi_char);

		fin.seekg((i + 1) * 80);
		fin.read(vis_mag_char, 10);
		vis_mag2 = atof(vis_mag_char);

		fin.seekg(f, ios::cur);
		fin.read(fi_char, 10);
		fi2 = atof(fi_char);

		if ((vis_mag0 + step / 2) >= vis_mag1 && (vis_mag0 + step / 2) <= vis_mag2)
		{
			check= line_interp(vis_mag1, fi1, vis_mag2, fi2, (vis_mag0 + step / 2))*step*ALPHA;
			return check;
			break;
		}
	}
}
void to_mf(int K, fstream& fini, double abs_mag0, double fi, double& mass, double& MASS_DOUBLE, double& MASS_DOUBLE_i, int k)
{
	char mass_char[8], abs_mag_char[8];
	double mass1, mass2, abs_mag1, abs_mag2;
	for (int j = 0; j <= K - 2; j++)  // нахождение массы, соответствующей абс звездной величине серидины интервала
	{
		fini.seekg(j * 14);
		fini.read(mass_char, 5);
		mass1 = atof(mass_char);

		fini.seekg(1, ios::cur);
		fini.read(abs_mag_char, 6);
		abs_mag1 = atof(abs_mag_char);

		fini.seekg((j + 1) * 14);
		fini.read(mass_char, 5);
		mass2 = atof(mass_char);

		fini.seekg(1, ios::cur);
		fini.read(abs_mag_char, 6);
		abs_mag2 = atof(abs_mag_char);
		if (abs_mag0 <= abs_mag2 && abs_mag0 >= abs_mag1)//<> когда зв величины по возрастанию
		{
			mass = line_interp(abs_mag1, mass1, abs_mag2, mass2, abs_mag0);
			MASS_DOUBLE_i = fi*mass * k;
			MASS_DOUBLE = MASS_DOUBLE + MASS_DOUBLE_i;
			break;
		}
	}
}

int main()
{
	int N = 44, K = 58, Nsteps=30;//потом с клавиатуры
	double  mod_dist = 10.2, EBV = 0.2, ALPHA = 0., step;// константы для скопления - потом с клавиатуры
	char vis_mag_char[10];
	double vis_mag0, vis_mag_fin, fi, abs_mag0, mass;
	double MASS_DOUBLE_i, MASS_DOUBLE = 0, MASS_SINGLE_i, MASS_SINGLE = 0, MASS_CLUSTER_NEW, MASS_CLUSTER_OLD_i, MASS_CLUSTER_OLD = 0, NEW_divided_by_OLD;
	double MASS_DOUBLE_i_up, MASS_DOUBLE_up = 0, MASS_SINGLE_i_up, MASS_SINGLE_up = 0, MASS_CLUSTER_NEW_up, MASS_DOUBLE_i_low, MASS_DOUBLE_low = 0, MASS_SINGLE_i_low, MASS_SINGLE_low = 0, MASS_CLUSTER_NEW_low;
		fstream fin, fini, fout, fout2;
	fin.open("density.txt", ios::in);
	fini.open("mass-absmag.txt", ios::in);
	fout.open("number.txt", ios_base::out);
	fout2.open("alpha.txt", ios_base::out);
	fout2.setf(ios::fixed);
	fout2.setf(ios::fixed);
	fout2 << "ALPHA" << '\t' << "MASS_CLUSTER_NEW" << '\t' << "MASS_CLUSTER_NEW_low" << '\t' << "MASS_CLUSTER_NEW_up" << '\t' << '\t' << "NEW_divided_by_OLD" << endl;
	for (ALPHA; ALPHA <= 0.9; ALPHA = ALPHA + 0.1)
	{
		fin.seekg(0);
		fin.read(vis_mag_char, 10);
		vis_mag0 = atof(vis_mag_char); // начальная звездная величина

		fin.seekg((N - 1) * 80);
		fin.read(vis_mag_char, 10);
		vis_mag_fin = atof(vis_mag_char);// конечная звездная величина

		step = (vis_mag_fin - vis_mag0) / Nsteps;

		MASS_DOUBLE = 0;  MASS_SINGLE = 0; MASS_CLUSTER_OLD = 0;	MASS_DOUBLE_low = 0;  MASS_SINGLE_low = 0;	MASS_DOUBLE_up = 0;  MASS_SINGLE_up = 0;
		//масса
		while (vis_mag0 +step/2 <= vis_mag_fin)
		{
			//для только одиночных
			fi = round (lf(N, 7, fin, vis_mag0, step, 1));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную одиночной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_CLUSTER_OLD, MASS_CLUSTER_OLD_i, 1);
			fout << MASS_CLUSTER_OLD << endl;
			//для доли двойных
			fi = round (lf(N, 7, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) + 2.5 * log10(2.) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную двойной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_DOUBLE, MASS_DOUBLE_i, 2);
			//fout << abs_mag0 << '\t' << fi << '\t' << mass << '\t' << MASS_DOUBLE_i << endl;

			//для доли одиночных 
			fi= round(lf(N, 7, fin, vis_mag0, step, 1))- round(lf(N, 7, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную одиночной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_SINGLE, MASS_SINGLE_i, 1);

			vis_mag0 = vis_mag0 + step;
		} 
		
		fin.seekg(0);
		fin.read(vis_mag_char, 10);
		vis_mag0 = atof(vis_mag_char); // начальная звездная величина

		//нижняя масса
		while (vis_mag0 + step <= vis_mag_fin)
		{
			//для доли двойных
			fi = round(lf(N, 24, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) + 2.5 * log10(2.) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную двойной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_DOUBLE_low, MASS_DOUBLE_i_low, 2);

			//для доли одиночных 
			fi = round(lf(N, 24, fin, vis_mag0, step, 1)) - round(lf(N, 24, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную одиночной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_SINGLE_low, MASS_SINGLE_i_low, 1);
			vis_mag0 = vis_mag0 + step;
		}
		fin.seekg(0);
		fin.read(vis_mag_char, 10);
		vis_mag0 = atof(vis_mag_char); // начальная звездная величина

		// верхняя масса
		while (vis_mag0 + step <= vis_mag_fin)
		{
			//для доли двойных
			fi = round(lf(N, 41, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) + 2.5 * log10(2.) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную двойной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_DOUBLE_up, MASS_DOUBLE_i_up, 2);

			//для доли одиночных 
			fi = round(lf(N, 41, fin, vis_mag0, step, 1)) - round(lf(N, 41, fin, vis_mag0, step, ALPHA));
			abs_mag0 = (vis_mag0 + step / 2) - 2.43*0.37*EBV - mod_dist;//перевод из видимой одной звезды в абсолютную одиночной
			to_mf(K, fini, abs_mag0, fi, mass, MASS_SINGLE_up, MASS_SINGLE_i_up, 1);
			vis_mag0 = vis_mag0 + step;
		}
		
		MASS_CLUSTER_NEW = MASS_DOUBLE + MASS_SINGLE;
		MASS_CLUSTER_NEW_low = MASS_DOUBLE_low + MASS_SINGLE_low;
		MASS_CLUSTER_NEW_up = MASS_DOUBLE_up + MASS_SINGLE_up;
		NEW_divided_by_OLD = MASS_CLUSTER_NEW / MASS_CLUSTER_OLD;
		//fout << ALPHA << '\t' << MASS_CLUSTER_NEW << '\t' << MASS_CLUSTER_OLD << '\t' << NEW_divided_by_OLD << endl << endl << endl << endl;
		fout2 << ALPHA << '\t' << MASS_CLUSTER_NEW << '\t' << MASS_CLUSTER_NEW_low << '\t' << MASS_CLUSTER_NEW_up << '\t' << '\t' << NEW_divided_by_OLD <<endl;
	}
	 
	fin.close();
	fini.close();
	fout.close();
	system("pause");
	return 0;
}
*/