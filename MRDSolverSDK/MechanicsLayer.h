#pragma once
//标记点Marker2D类
//类内成员数据
class Marker2D
{
public:
	double X[2];//1-1标记点坐标
	bool IfBoundary;//1-2-1内部点/边界点
	double n[2];//1-2-2方向
	double r;//2-1临近域半径
	int Nei;//2-2临近点数
	bool Fixed[2];//3-1是否固定
	double Force[2];//3-2载荷
	double R[2];//4-1余量
	double displacement[2];//标记点位移5-1
	double epsilon[2][2];//标记点应变5-2
	double sigma[2][2];//标记点应力5-3
};

//标记点Marker3D类
//类内成员数据
class Marker3D
{
public:
	double X[3];//1-1标记点坐标
	bool IfBoundary;//1-2-1内部点/边界点
	double n[3];//1-2-2方向
	double r;//2-1临近域半径
	int Nei;//2-2临近点数
	bool Fixed[3];//3-1是否固定
	double Force[3];//3-2载荷
	double R[3];//4-1余量
	double displacement[3];//标记点位移5-1
	double epsilon[3][3];//标记点应变5-2
	double sigma[3][3];//标记点应力5-3
};

//类外全局变量
//2-3number of Markers
extern int Num_m;
//2-4Number of Neighbors
extern int* N_Nei;
//2-5GlobalIndex_Nei[m][n] is the global index of the n-th neighbor of the m-th Marker
extern int** GlobalIndex_Nei;
//2-6Vect[m][Nei][2] is the connectivity matrix of the m-th Marker
extern double*** Vect;
extern double E;//0-1弹性模量
extern double nu;//0-2泊松比
extern double G;//0-3剪切模量

//类外函数
//2-1几何方程
void Strain_from_Displacement(double** epsilon, double** Vec, double** disp_nei, double* disp_m, int n_nei, int Dim, double tolerance);

//2-2本构方程
void Stress_from_Strain(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//2D平面应力
void Stress_from_Strain2D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);
//3D
void Stress_from_Strain3D(double** sigma, double E, double nu, double G, double** epsilon, int Dim);

//2-3应力算内力
//2-3-1应力算内力-内部点
void InForce_from_Stress_in(double* inforce, double** Vec, double*** sigma_nei, double** sigma_m, int n_nei, int Dim, double tolerance);
//2-3-2应力算内力-边界点
void InForce_from_Stress_bdy(double* inforce, double** sigma_m, double* n, int Dim);