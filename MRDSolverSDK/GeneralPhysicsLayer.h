#pragma once
//3-1����marker m�����Ա���phi���ݶ�, �����ݶ�Gradphi[]��ָ�����ΪVec[][]���ٽ�������phi_nei[]�������ǵ�����phi_m���ٽ�����n_nei��ά��Dim��
void Gradient_of_Scalar(double* Gradphi, double** Vec, double* phi_nei, double phi_m, int n_nei, int Dim, double tolerance);

//3-2����marker m������ʸ��E���ݶ�, �����ݶ�GradE[][]��ָ�����ΪVec[][]���ٽ�������phi_nei[][]�������ǵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Gradient_of_Vector(double** GradE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);

//3-3������ɢ������marker m������ʸ��E��ɢ��, ����ɢ��DivE��ָ�����ΪVec[][]���ٽ�������phi_nei[][]�������ǵ�����phi_m[]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Vector(double &DivE, double** Vec, double** E_nei, double* E_m, int n_nei, int Dim, double tolerance);

//3-4������ɢ������marker m����������T��ɢ��, ����ɢ��GradE��ָ�����ΪVec[][]���ٽ�������T_nei[][][]�������ǵ�����T_m[][]���ٽ�����n_nei��ά��Dim��
void Divergence_of_Tensor(double* GradT, double** Vec, double*** T_nei, double** T_m, int n_nei, int Dim, double tolerance);