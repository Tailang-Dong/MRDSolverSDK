#pragma once
//����ϵ������ľ��������Ax=b
//xΪ��������ָ��AΪϵ�������ָ��bΪ����������ָ��nΪ����ά��tolranceΪԤ�辫��

//4-1Jacobi����������
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance);

//4-2Guass-Seidel�����������������
void Guass_Seidel_Slover(double* x, double** A, double* b, int n, double tolerance);
//4-2-2�ع�һ��toleranceȱʡֵΪ1e-7
void Guass_Seidel_Slover(double* x, double** A, double* b, int n);

//4-3��m>��n,��С���˽�,ATA��������ݲ�
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance);
//4-3-2�ع���m>��n,��С���˽�ATA��������ݲ�Ĭ��1e-7
void LeastSquareSolution(double* x, double** A, double* b, int m, int n);