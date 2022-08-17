#include "MathLayer.h"
#include  <iostream>
#include <cmath>
using namespace std;

//Guass-Seidel迭代法求线性代数方程组,系数矩阵为方阵
void Guass_Seidel_Slover(double* x,double** A, double* b, int n, double tolerance)
{
	/******Step1 设置解的迭代初值-用户已给******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 设置函数体需要的局部变量******/
	double* x0 = new double[n];//进行迭代之前，把当前的根存起来，求误差范数用
	double tol;//误差范数
	int it_sm = 0;//迭代计次器

	//如果长时间不收敛，break一下
	//cout << "您调用了Gusas-Seidel迭代法求解器"  << endl;

	/******Step3 开始迭代******/
	do
	{
		//存一下第it步的解
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//开始第it步的迭代
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				x[i] =x[i] - A[i][j] * x[j] / A[i][i];
				
			}

		}
		//第it次迭代完成

		/******Step4 计算迭代误差******/
		//计算误差范数tol=||x-x0||max范数，所有分量的绝对值必须小于容忍公差
		tol = x[0] - x0[0];//第0分量误差
		tol = abs(tol);//取绝对值
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it_sm = it_sm + 1;//迭代计次器+1

		/******Step5 输出计算log日志it与tol******/
		//cout << "it=" <<it_sm << "  " << "tol=" << tol <<endl;

	} while (tol > tolerance);

	//cout << "it=" << it_sm << "  收敛" << endl;//同为日志文件输出

	/******Step6 释放局部的动态数组内存******/
	delete [] x0;//释放动态数组

}

//Guass-Seidel迭代法求线性代数方程组-缺省预设精度重载函数，迭代输出语句后期删除
void Guass_Seidel_Slover(double* x, double** A, double* b, int n)
{
	/******Step1 设置解的迭代初值-用户已给******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 设置函数体需要的局部变量******/
	double tolerance=1e-7;//
	double* x0 = new double[n];//进行迭代之前，把当前的根存起来，求误差范数用
	double tol;//误差范数
	int it = 0;//迭代计次器
	cout << "您调用时没有预设精度要求，默认预设精度tolerance=" << tolerance << endl;

	//如果长时间不收敛，break一下，算了先不加此功能了。这只是一个函数

	/******Step3 开始迭代******/
	do
	{
		//存一下第it步的解
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//开始第it步的迭代
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x[j] / A[i][i];

			}

		}
		//第it次迭代完成

		/******Step4 计算迭代误差******/
		//计算误差范数tol=||x-x0||1范数，所有分量的绝对值必须小于容忍公差
		tol = x[0] - x0[0];//第0分量误差
		tol = abs(tol);//取绝对值
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//迭代计次器+1

		/******Step5 输出计算log日志it与tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  收敛" << endl;//同为日志文件输出

	/******Step6 释放局部的动态数组内存******/
	delete[] x0;//释放动态数组

}

//Jacobi迭代法求线性代数方程组
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance)
{
	/******Step1 设置解的迭代初值******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 设置函数体需要的局部变量******/
	double* x0 = new double[n];//进行迭代之前，把当前的根存起来，求误差范数用
	double tol;//误差范数
	int it = 0;//迭代计次器

	//如果长时间不收敛，break一下
	cout << "您调用了Jacobi迭代法求解器" << endl;

	/******Step3 开始迭代******/
	do
	{
		//存一下第it步的解
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//开始第it步的迭代
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x0[j] / A[i][i];

			}

		}
		//第it次迭代完成

		/******Step4 计算迭代误差******/
		//计算误差范数tol=||x-x0||max范数，所有分量的绝对值必须小于容忍公差
		tol = x[0] - x0[0];//第0分量误差
		tol = abs(tol);//取绝对值
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol > temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//迭代计次器+1

		/******Step5 输出计算log日志it与tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  收敛" << endl;//同为日志文件输出

	/******Step6 释放局部的动态数组内存******/
	delete[] x0;//释放动态数组

}

//行m>列n,最小二乘解
void LeastSquareSolution(double* x, double** A, double* b, int m, int n)
{
	//临时变量VTV和VTphi
	double** ATA = new double* [n];
	double* ATb = new double[n];
	for (int i = 0; i < n; i++)
		ATA[i] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ATA[i][j] = 0;
			for (int k = 0; k < m; k++)
				ATA[i][j] += A[k][i] * A[k][j];
		}

		ATb[i] = 0;
		for (int k = 0; k < m; k++)
			ATb[i] += A[k][i] * b[k];
	}
	//调用GuassSeidel矩阵求解器求theta
	Guass_Seidel_Slover(x, ATA, ATb, n);
	//释放临时变量VTV和VTphi
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}

//行m>列n,最小二乘解
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance)
{
	//临时变量VTV和VTb
	double** ATA = new double* [n];
	double* ATb = new double[n];
	for (int i = 0; i < n; i++)
		ATA[i] = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ATA[i][j] = 0;
			for (int k = 0; k < m; k++)
				ATA[i][j] += A[k][i] * A[k][j];
		}

		ATb[i] = 0;
		for (int k = 0; k < m; k++)
			ATb[i] += A[k][i] * b[k];
	}
	//调用GuassSeidel矩阵求解器求x
	Guass_Seidel_Slover(x, ATA, ATb, n, tolerance);
	//释放临时变量VTV和VTb
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}
