#pragma once
//��ǵ�Marker2D��
//���ڳ�Ա����
class Marker2D
{
public:
	double X[2];//1-1��ǵ�����
	bool IfBoundary;//1-2-1�ڲ���/�߽��
	double n[2];//1-2-2����
	double r;//2-1�ٽ���뾶
	int Nei;//2-2�ٽ�����
	bool Fixed[2];//3-1�Ƿ�̶�
	double Force[2];//3-2�غ�
	double R[2];//4-1����
	double displacement[2];//��ǵ�λ��5-1
	double epsilon[2][2];//��ǵ�Ӧ��5-2
	double sigma[2][2];//��ǵ�Ӧ��5-3
};

//��ǵ�Marker3D��
//���ڳ�Ա����
class Marker3D
{
public:
	double X[3];//1-1��ǵ�����
	bool IfBoundary;//1-2-1�ڲ���/�߽��
	double n[3];//1-2-2����
	double r;//2-1�ٽ���뾶
	int Nei;//2-2�ٽ�����
	bool Fixed[3];//3-1�Ƿ�̶�
	double Force[3];//3-2�غ�
	double R[3];//4-1����
	double displacement[3];//��ǵ�λ��5-1
	double epsilon[3][3];//��ǵ�Ӧ��5-2
	double sigma[3][3];//��ǵ�Ӧ��5-3
};

//����ȫ�ֱ���
//2-3number of Markers
extern int Num_m;
//2-4Number of Neighbors
extern int* N_Nei;
//2-5GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Marker
extern int** GlobalIndex_Nei;
//2-6Vect[m][Nei][2] is the connectivity matrix of the m-th Marker
extern double*** Vect;
extern double E;//0-1����ģ��
extern double nu;//0-2���ɱ�
extern double G;//0-3����ģ��

//���⺯��
//2-1���η���
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, double tolerance);

//2-2��������
void Stress_from_Strain(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//2Dƽ��Ӧ��
void Stress_from_Strain2D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//3D
void Stress_from_Strain3D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);

//2-3Ӧ��������
//2-3-1Ӧ��������-�ڲ���
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, double tolerance);
//2-3-2Ӧ��������-�߽��
void InForce_from_Stress_bdy(double* inforce, double** sigma_m, double* n, int Dim);