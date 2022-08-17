#include "MathLayer.h"
#include  <iostream>
#include <cmath>
using namespace std;

//Guass-Seidel�����������Դ���������,ϵ������Ϊ����
void Guass_Seidel_Slover(double* x,double** A, double* b, int n, double tolerance)
{
	/******Step1 ���ý�ĵ�����ֵ-�û��Ѹ�******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 ���ú�������Ҫ�ľֲ�����******/
	double* x0 = new double[n];//���е���֮ǰ���ѵ�ǰ�ĸ�����������������
	double tol;//����
	int it_sm = 0;//�����ƴ���

	//�����ʱ�䲻������breakһ��
	//cout << "��������Gusas-Seidel�����������"  << endl;

	/******Step3 ��ʼ����******/
	do
	{
		//��һ�µ�it���Ľ�
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//��ʼ��it���ĵ���
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				x[i] =x[i] - A[i][j] * x[j] / A[i][i];
				
			}

		}
		//��it�ε������

		/******Step4 ����������******/
		//��������tol=||x-x0||max���������з����ľ���ֵ����С�����̹���
		tol = x[0] - x0[0];//��0�������
		tol = abs(tol);//ȡ����ֵ
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it_sm = it_sm + 1;//�����ƴ���+1

		/******Step5 �������log��־it��tol******/
		//cout << "it=" <<it_sm << "  " << "tol=" << tol <<endl;

	} while (tol > tolerance);

	//cout << "it=" << it_sm << "  ����" << endl;//ͬΪ��־�ļ����

	/******Step6 �ͷžֲ��Ķ�̬�����ڴ�******/
	delete [] x0;//�ͷŶ�̬����

}

//Guass-Seidel�����������Դ���������-ȱʡԤ�辫�����غ������������������ɾ��
void Guass_Seidel_Slover(double* x, double** A, double* b, int n)
{
	/******Step1 ���ý�ĵ�����ֵ-�û��Ѹ�******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 ���ú�������Ҫ�ľֲ�����******/
	double tolerance=1e-7;//
	double* x0 = new double[n];//���е���֮ǰ���ѵ�ǰ�ĸ�����������������
	double tol;//����
	int it = 0;//�����ƴ���
	cout << "������ʱû��Ԥ�辫��Ҫ��Ĭ��Ԥ�辫��tolerance=" << tolerance << endl;

	//�����ʱ�䲻������breakһ�£������Ȳ��Ӵ˹����ˡ���ֻ��һ������

	/******Step3 ��ʼ����******/
	do
	{
		//��һ�µ�it���Ľ�
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//��ʼ��it���ĵ���
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x[j] / A[i][i];

			}

		}
		//��it�ε������

		/******Step4 ����������******/
		//��������tol=||x-x0||1���������з����ľ���ֵ����С�����̹���
		tol = x[0] - x0[0];//��0�������
		tol = abs(tol);//ȡ����ֵ
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol < temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//�����ƴ���+1

		/******Step5 �������log��־it��tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  ����" << endl;//ͬΪ��־�ļ����

	/******Step6 �ͷžֲ��Ķ�̬�����ڴ�******/
	delete[] x0;//�ͷŶ�̬����

}

//Jacobi�����������Դ���������
void Jacobi_Solver(double* x, double** A, double* b, int n, double tolerance)
{
	/******Step1 ���ý�ĵ�����ֵ******/
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}

	/******Step2 ���ú�������Ҫ�ľֲ�����******/
	double* x0 = new double[n];//���е���֮ǰ���ѵ�ǰ�ĸ�����������������
	double tol;//����
	int it = 0;//�����ƴ���

	//�����ʱ�䲻������breakһ��
	cout << "��������Jacobi�����������" << endl;

	/******Step3 ��ʼ����******/
	do
	{
		//��һ�µ�it���Ľ�
		for (int k = 0; k < n; k++)
		{
			x0[k] = x[k];
		}

		//��ʼ��it���ĵ���
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
					x[i] = x[i] - A[i][j] * x0[j] / A[i][i];

			}

		}
		//��it�ε������

		/******Step4 ����������******/
		//��������tol=||x-x0||max���������з����ľ���ֵ����С�����̹���
		tol = x[0] - x0[0];//��0�������
		tol = abs(tol);//ȡ����ֵ
		for (int m = 1; m < n; m++)
		{
			double temp_tol = x[m] - x0[m];
			temp_tol = abs(temp_tol);
			if (tol > temp_tol)
			{
				tol = temp_tol;
			}
		}

		it = it + 1;//�����ƴ���+1

		/******Step5 �������log��־it��tol******/
		cout << "it=" << it << "  " << "tol=" << tol << endl;

	} while (tol > tolerance);

	cout << "it=" << it << "  ����" << endl;//ͬΪ��־�ļ����

	/******Step6 �ͷžֲ��Ķ�̬�����ڴ�******/
	delete[] x0;//�ͷŶ�̬����

}

//��m>��n,��С���˽�
void LeastSquareSolution(double* x, double** A, double* b, int m, int n)
{
	//��ʱ����VTV��VTphi
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
	//����GuassSeidel�����������theta
	Guass_Seidel_Slover(x, ATA, ATb, n);
	//�ͷ���ʱ����VTV��VTphi
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}

//��m>��n,��С���˽�
void LeastSquareSolution(double* x, double** A, double* b, int m, int n, double tolerance)
{
	//��ʱ����VTV��VTb
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
	//����GuassSeidel�����������x
	Guass_Seidel_Slover(x, ATA, ATb, n, tolerance);
	//�ͷ���ʱ����VTV��VTb
	for (int i = 0; i < n; i++)
		delete[] ATA[i];
	delete[] ATA;
	delete[] ATb;
}
