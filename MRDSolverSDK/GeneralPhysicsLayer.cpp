#include "GeneralPhysicsLayer.h"
#include "MathLayer.h"
#include <iostream>

//3-1标量梯度
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance)
{
	double* dphi = new double[n_nei];//3-1-1设临时变量dphi[]
	for (int n = 0; n < n_nei; n++)
		dphi[n] = phi_nei[n] - phi_m;
	LeastSquareSolution(Gradphi, Vec, dphi, n_nei, Dim, tolerance);//3-1-2调用4-3求解梯度
	delete[] dphi;//3-1-3释放临时变量dphi[]

}

//3-2矢量E的梯度 E_nei[n_nei][Dim]
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-2-1设临时变量Gradphi[] phi_nei[] phi_m
	double* Gradphi = new double[Dim];//矢量的每个分量作为标量来分别计算梯度
	double* phi_nei = new double[n_nei];
	double phi_m;
	//3-2-2调用3-1计算每个分量的梯度
	for (int d_attriV = 0; d_attriV < Dim; d_attriV++)//Vector attribute第d_attriV个分量,逐个分量计算
	{
		for (int n = 0; n < n_nei; n++)//Vector变量的第d_attriV个分量看作一个标量phi
		{
			phi_nei[n] = E_nei[n][d_attriV];
			phi_m = E_m[d_attriV];
		}
		Gradient_of_Scalar(Gradphi, Vec, phi_nei, phi_m, n_nei, Dim, tolerance);//Vector变量的第d_attriV个分量的梯度
		for (int dcol = 0; dcol < Dim; dcol++)
			GradE[d_attriV][dcol] = Gradphi[dcol]; //赋值给张量GradE[][]的第d_attriV行
	}
	//3-2-3释放局部变量内存
	delete[] Gradphi;
	delete[] phi_nei;

}

//3-3矢量E的散度
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance)
{
	//3-3-1设临时变量GradE[][]
	double** GradE = new double* [Dim];
	for (int i = 0; i < Dim; i++)
		GradE[i] = new double[Dim];
	//3-3-2调用3-2计算矢量梯度GradE
	Gradient_of_Vector(GradE, Vec, E_nei, E_m, n_nei, Dim, tolerance);
	//3-3-3trace of GradE=DivE
	DivE = 0;
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
		{
			if (j == i)
				DivE = DivE + GradE[i][j];
		}
	//3-3-4释放临时内存GradE[][]
	for (int i = 0; i < Dim; i++)
		delete[] GradE[i];
	delete[] GradE;
}

//3-4张量的散度
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance)
{
	//每一列的散度T=[E0, E1,]
	//3-4-1临时变量Ej_nei[n_nei][Dim], Ej_m[Dim]
	double DivEj=0;
	double** Ej_nei = new double* [n_nei];
	for (int d = 0; d < n_nei; d++)
		Ej_nei[d] = new double[Dim];
	double* Ej_m = new double[Dim];

	for (int j = 0; j < Dim; j++)
	{
		for (int d = 0; d < n_nei; d++)
			for (int i = 0; i < Dim; i++)
				Ej_nei[d][i] = T_nei[d][i][j];//第j列
		for (int i = 0; i < Dim; i++)
			Ej_m[i] = T_m[i][j];

		//3-4-2调用3-3求第j列的散度
		Divergence_of_Vector(DivEj,Vec, Ej_nei, Ej_m, n_nei, Dim, tolerance);

		//3-4-3分量赋值
		GradT[j] = DivEj;
	}
	//3-4-4释放临时变量
	for (int d = 0; d < n_nei; d++)
		delete[] Ej_nei[d];
	delete[] Ej_nei;
	delete[] Ej_m;

}