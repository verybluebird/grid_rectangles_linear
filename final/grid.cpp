﻿#include "Grid.h"
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

void Grid::input()
{
	fstream file("domain_rz.txt");
	if (file.is_open()) {
		file >> Rw;
		file >> Rb;
		file >> NL;
		H.resize(NL);
		K.resize(NL);
		Phi.resize(NL);
		S2.resize(NL);
		for (int i = 0; i < NL; i++)
			file >> H[i];

		for (int i = 0; i < NL; i++)
		{
			file >> K[i];
			file >> Phi[i];
			file >> S2[i];
		}
		file.close();
	}
	else cout << "File domain_rz.txt is not opened!\n\n" << endl;

	file.open("mesh_rz.txt");
	if (file.is_open()) {
		file >> nr;
		file >> kr;
		nz.resize(NL);
		for (int i = 0; i < NL; i++)
			file >> nz[i];

		file >> M;
		file.close();
	}
	else cout << "File mesh_rz.txt is not opened!\n\n" << endl;

	file.open("zone_perf_rz.txt");
	if (file.is_open()) {
		file >> Nzp;
		Pu.resize(Nzp);
		Pd.resize(Nzp);
		Tetta.resize(Nzp);
		for (int i = 0; i < Nzp; i++)
		{
			file >> Pu[i];
			file >> Pd[i];
			file >> Tetta[i];
		}
		file.close();
	}
	else cout << "File zone_perf_rz.txt is not opened!\n\n" << endl;

	file.open("phaseprop.txt");
	if (file.is_open()) {
		file >> Nph;
		Mu.resize(Nph);

		for (int i = 0; i < Nph; i++)
			file >> Mu[i];

		file.close();
	}
	else cout << "File phaseprop.txt is not opened!\n\n" << endl;

	file.open("plast.txt");
	if (file.is_open())
	{
		file >> Plast;
		file.close();
	}
	else cout << "File plast.txt is not opened!\n\n" << endl;

}

void Grid::add_if_not_exist_and_sort(double L)
{
	if (find(z_coord.begin(), z_coord.end(), L) == z_coord.end())
	{
		// узла нет
		z_coord.push_back(L);
		sort(z_coord.begin(), z_coord.end());
	}
	else
		sort(z_coord.begin(), z_coord.end());
}


void Grid::nodes()
{
	double l = Rb - Rw;
	double h;
	if (kr < 1)
		h = l * (1 - kr) / (1 - pow(kr, nr)); // длина первого справа прямоугольника в сетке
	else if (kr == 1)
		h = l / nr;
	else cout << "kr should be <=1\n";
	int Nz_uz = 1; // кол-во узлов по z
	for (int i = 0; i < NL; i++)
		Nz_uz += nz[i];
	
	double S = 0; // начальная координата

	r_coord.resize(nr + 1); // координаты сетки по r
	r_coord[0] = Rw;
	S = r_coord[0];
	int k = 1;
	for (int i = nr - 1; i >= 0; i--)
	{
		S += h * pow(kr, i);
		r_coord[k] = S;
		k++;
	}
	
	
	z_coord.resize(Nz_uz); // координаты сетки по z
	S = z0;
	z_coord[0] = z0;
	int c = 1;
	for (int i = 0; i < NL; i++)
	{
		double delta = H[i] / nz[i];
		for (int j = 0; j < nz[i]; j++)
		{
			S += delta;
			z_coord[c] = S;
			c++;
		}
	}
	// добавляем узлы для перфорации
	for (int i = 0; i < Nzp; i++)
	{
		add_if_not_exist_and_sort(Pu[i]);
		add_if_not_exist_and_sort(Pd[i]);
	}

	// удаляем повторы в r
	sort(r_coord.begin(), r_coord.end());
	r_coord.erase(unique(r_coord.begin(), r_coord.end()), r_coord.end());
	//делаем вложенную сетку, если надо 
	if (M == 1) {
		nested_grid(r_coord);
		nested_grid(z_coord);
	}
	if (M == 2) {
		nested_grid(r_coord);
		nested_grid(z_coord);
		nested_grid(r_coord);
		nested_grid(z_coord);
	}

	Nuz = z_coord.size() * r_coord.size();
	// вывод в файл
	ofstream file("node.txt");
	if (file.is_open())
	{
		file << Nuz << endl;
		file.precision(16);
		for (int i = 0; i < z_coord.size(); i++)
			for (int j = 0; j < r_coord.size(); j++)
				file << r_coord[j] << " " << z_coord[i] << endl;
		file.close();
	}
	else cout << "File node.txt is not opened!\n\n" << endl;
}

void Grid::elems()
{
	Nel = (z_coord.size() - 1) * (r_coord.size() - 1);
	global_numbers.resize(z_coord.size());
	int c = 0;
	for (int i = 0; i < z_coord.size(); i++)
	{
		global_numbers[i].resize(r_coord.size());
		for (int j = 0; j < r_coord.size(); j++, c++)
			global_numbers[i][j] = c;
	}
	// вывод номеров прямоугольников
	ofstream file("elem.txt ");
	if (file.is_open())
	{
		file << Nel << endl;
		for (int i = 0; i < z_coord.size() - 1; i++)
		{
			for (int j = 0; j < r_coord.size() - 1; j++)
				file << global_numbers[i][j] << " " << global_numbers[i][j + 1] << " " << global_numbers[i + 1][j] << " " << global_numbers[i + 1][j + 1] << endl;
		}
		file.close();
	}
	else cout << "File elem.txt is not opened!\n\n" << endl;
}

void Grid::material()
{
	ofstream out;
	out.open("mat.txt");
	int k = 0;
	double boarder = H[0];
	for (int i = 0; i < z_coord.size() - 1; i++)
	{
		if (z_coord[i] >= boarder) // если граница слоя пройдена, переходим к следующему слою
		{
			boarder += H[k];
			k++;
		}
		for (int j = 0; j < r_coord.size() - 1 && k<NL; j++)
			out << K[k] << ' ' << Phi[k] << ' ' << S2[k] << endl;
	}
	out.close();
}

void Grid::gr_bc1()
{
	Nbc1 = z_coord.size();
	ofstream out;
	out.open("bc1.txt");
	out << Nbc1 << endl;
	for (int i = 1; i <= z_coord.size(); i++)
		out << i * r_coord.size() - 1 << ' ' << Plast << endl;
	out.close();
}

void Grid::gr_bc2()
{
	Nbc2 = 0;
	vector<int> num_zp_elems (Nzp, 0); // кол-во КЭ в ЗП
	vector<int> first_zp_elem_num(Nzp, 0); // номер первого КЭ в ЗП
	int i = 0;
	for (int k = 0; k < Nzp; k++)
	{
		for (i = 0; z_coord[i] < Pu[k]; i++);
			
				first_zp_elem_num[k] = ((r_coord.size()) - 1) * i;
		num_zp_elems[k]++;
		i++;
		for (; z_coord[i] < Pd[k]; i++) num_zp_elems[k]++;
		Nbc2 += num_zp_elems[k];
	}
	ofstream out;
	out.open("bc2.txt");
	out << Nbc2 << endl;
	for (int k = 0; k < Nzp; k++)
		for (int j = 0; j < num_zp_elems[k]; j++)
			out << first_zp_elem_num[k] + ((r_coord.size()) - 1) * j << " " << 0 << " " << 2 << " " << Tetta[k] << endl;
	out.close();
}

void Grid::nested_grid(vector<double>& coord)
{
	vector<double> coord2((coord.size() - 1) * 2 + 1);
	coord2[0] = coord[0];
	int k = 1;
	for (int i = 1; i < coord.size(); i++)
	{
		coord2[k] = (coord[i-1] + coord[i]) / 2;
		k++;
		coord2[k] = coord[i];
		k++;
	}
	coord = coord2;
	coord2.clear();
}

void Grid::print_profile()
{
	vector<vector<int>> c_list; // connections list
	c_list.resize(Nuz);
	// формирование списка связности
	for (int i = 1; i <= r_coord.size(); i++)
		c_list[i].push_back(i-1);
	for (int j = 1; j < z_coord.size(); j++)
	{
		c_list[j * r_coord.size()].push_back((j - 1) * r_coord.size());
		c_list[j * r_coord.size()].push_back((j - 1) * r_coord.size() + 1);

		c_list[(j + 1) * r_coord.size() - 1].push_back(j * r_coord.size() - 2);
		c_list[(j + 1) * r_coord.size() - 1].push_back(j * r_coord.size() - 1);
		c_list[(j + 1) * r_coord.size() - 1].push_back((j + 1) * r_coord.size());

		for (int i = j * r_coord.size()+1; i < (j + 1) * r_coord.size() - 2; i++)
		{
			c_list[i - 1].push_back(i - r_coord.size() - 1);
			c_list[i - 1].push_back(i - r_coord.size());
			c_list[i - 1].push_back(i - r_coord.size()+1);
			c_list[i - 1].push_back(i - 1);
		}
	}

	ofstream out;
	out.open("profile.txt"); // разреженный строчно-столбцовый формат (массивы ia, ja)
	out << "0 ";
	int counter = 0;
	for (int i = 0; i < c_list.size(); i++)
	{
		counter += c_list[i].size();
		out << counter << ' ';
	}
	out << endl;
	for (int i = 0; i < c_list.size(); i++)
		for (int j = 0; j < c_list[i].size(); j++)
			out << c_list[i][j] << ' ';
	out << endl;
	out.close();
}