#pragma once
#include <vector>
using namespace std;
class Grid
{
private:
	/*��������� ��������� �������*/
	double Rw; //������ ��������
	double Rb; //���������� �� ��������� ������� (����)
	int NL; //���������� �����
	vector<double> H; // ������� �����
	vector<double> K; //����������� ������������� �����
	vector<double> Phi; // ���������� �����
	vector<double> S2; //����������������� �����
	int nr; //���-�� ��������� �� r
	double kr; //����. �������� �� r
	vector<int> nz; //���-�� ��������� �� z � ������ ����
	int M; // �������� �����
	int Nzp; //���������� ��� ����������
	vector<double> Pu; //������� ������� �.�.
	vector<double> Pd; //������ ������� �.�.
	vector<double> Tetta; //��������
	int Nph; // ���������� ���
	vector<double> Mu; //�������� ���
	double Plast; //��������� ��������
	vector<vector<int>> global_numbers; //���������� ������ �����
	double z0 = 0; //������ ���������� ������ ������ ���� �� z
	vector<double> r_coord; // ���������� �� r
	vector<double> z_coord; // ���������� �� z

public:
	void input();// ����  ������
	int Nuz; //���������� �����
	vector<vector<double>> coords;//{r_i y_i}
	int Nel; //
	vector<vector<int>> numbers; // ���������� ������
	vector<vector<double>> materials;
	int Nbc1; //���-�� ������ ������� �������
	vector<double> bc1;
	int Nbc2;//���-�� ������ ������� �������
	vector<vector<double>> bc2;
	void nodes();
	void elems();
	void material();
	void gr_bc1();
	void gr_bc2();
	void nested_grid(vector<double>& coord); //����������� ���-�� ����� �����
	void add_if_not_exist_and_sort(double L); // �������� ����� ���� � �����, ���� ��� �� ����������
	void print_profile(); // ����������� ������� �������


};