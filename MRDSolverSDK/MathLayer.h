#pragma once
//方阵系数矩阵的矩阵求解器Ax=b
//x为解向量的指针A为系数矩阵的指针b为常数向量的指针n为向量维数tolrance为预设精度

//4-1Jacobi迭代法声明
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance);

//4-2Guass-Seidel迭代求解器函数声明
void Guass_Seidel_Slover(double* x, double** A, double* b, int n, double tolerance);
//4-2-2重构一下tolerance缺省值为1e-7
void Guass_Seidel_Slover(double* x, double** A, double* b, int n);

//4-3行m>列n,最小二乘解,ATA方阵迭代容差
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance);
//4-3-2重构行m>列n,最小二乘解ATA方阵迭代容差默认1e-7
void LeastSquareSolution(double* x, double** A, double* b, int m, int n);