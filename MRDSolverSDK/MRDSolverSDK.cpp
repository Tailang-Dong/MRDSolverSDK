// MRDSolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
//移动后非逐个更新

#include <iostream>
#include <fstream>
#include "ApplicationLayer.h"
#include "MechanicsLayer.h"
#include "GeneralPhysicsLayer.h"
#include "MathLayer.h"
#include "FileLayer.h"

using namespace std;

int main()
{
	/******0初始化材料信息******/
	G = E / 2 / (1 + nu);

	//模型文件读取选项设置
	Dimension = 2;//维数
	bool UseBinaryModel = 1;//是否使用二进制模型文件
	//求解信息
	double r = 0.06;//0.06
	double eta = 1.0;
	double eta_in =0.64;//边界标记点移动系数0.64
	double eta_b = 0.8;//0.8
	double K_in_2D = E * (3.0 - nu) / (1 - nu * nu) / r / r / 2;//内部移动刚度
	double K_b_2D = E / (1 - nu * nu) / r;//边界标记点移动刚度
	//计算设置
	alpha = eta_in / K_in_2D;//内部marker移动率
	beta = eta_b / K_b_2D;//外部marker移动率
	double ResiTolerance = 10;//余量判断容差##########全局
	double ResiToleranceInternal = 100;//内部点余量判断容差##########全局
	double ResiToleranceBounadry = r*ResiToleranceInternal;//内部点余量判断容差##########全局
	double DiffTolerance = 1.0e-20;//近似微分调用函数时迭代容差#######
	int it_orien = 0;//记录orientation次数##########全局
	int MaxIt_Rrien =110000;//防止遇到死循环
	//输出设置
	bool RecordRelaxation = 0;//松弛动画
	bool OutputLongitudinalSection = 0;//矩形梁的
	bool OutputCrossSection = 0;

	cout << "G=" << G << "\n";
	cout << "E=" << E << "\n";
	cout << "nu=" << nu << "\n";
	/******1初始化模型信息******/
	if (UseBinaryModel)
		mfin.open("model7.b", ios::binary);
	else
		mfin.open("model7.dat");//打开文件
	double L = 2.0;//特例
	
	/****1-1Marker数量****/
	if (UseBinaryModel)
		mfin.read((char*)&Num_m, sizeof Num_m);
	else
		mfin >> Num_m;
	cout << "Num_m=" << Num_m << "\n";
	/****1-2Marker数据****/
	m = new Marker2D[Num_m];//在应用层提前申请好
	if (UseBinaryModel)
		mfin.read((char*)&(m[0]), (long long)Num_m * sizeof(m[0]));
	else
	{
		for (int k = 0; k < Num_m; k++)
		{
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].X[d];//1-1几何坐标
			mfin >> m[k].IfBoundary;//1-2-1几何是否边界
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].n[d];//1-2-2几何边界外法向
			mfin >> m[k].r;//2-1离散连接性-临近域半径
			mfin >> m[k].Nei;//2-2离散连接性-临近点数
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].Fixed[d];////3-1约束-固定?
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].Force[d];//3-2载荷-载荷
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].R[d];//4余量
			for (int d = 0; d < Dimension; d++)
				mfin >> m[k].displacement[d];//5-1结果-位移
			for (int i = 0; i < Dimension; i++)
				for (int j = 0; j < Dimension; j++)
				{
					mfin >> m[k].epsilon[i][j];//5-2结果-应变
				}
			for (int i = 0; i < Dimension; i++)
				for (int j = 0; j < Dimension; j++)
				{
					mfin >> m[k].sigma[i][j];//5-3结果-应力
				}
		}
	}
	
	/****1-3connectivity数据****/
	N_Nei = new int[Num_m];
	if (UseBinaryModel)
		mfin.read((char*)N_Nei, (long long)Num_m * sizeof(N_Nei[0]));
	else
	{
		for (int k = 0; k < Num_m; k++)
			mfin >> N_Nei[k];
	}

	GlobalIndex_Nei = new int* [Num_m];
	for (int i = 0; i < Num_m; i++)
		GlobalIndex_Nei[i] = new int[N_Nei[i]];
	if (UseBinaryModel)
	{
		for (int i = 0; i < Num_m; i++)//**********************
		{
			mfin.read((char*)GlobalIndex_Nei[i], (long long)N_Nei[i] * sizeof(GlobalIndex_Nei[i][0]));
		}
	}
	else
	{
		for (int i = 0; i < Num_m; i++)
		{
			for (int j = 0; j < N_Nei[i]; j++)
				mfin >> GlobalIndex_Nei[i][j];
		}
	}

	mfin.close();//关闭数据文件
	/****1-4connectivity数据方向矩阵组****/
	Vect = new double** [Num_m];
	for (int p = 0; p < Num_m; p++)
		Vect[p] = new double* [N_Nei[p]];
	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < N_Nei[p]; i++)
			Vect[p][i] = new double[Dimension];
	for (int p = 0; p < Num_m; p++)
	{
		for (int i = 0; i < N_Nei[p]; i++)
		{
			int glodex = GlobalIndex_Nei[p][i];//p的第i个neighbor的全局序号globalindex
			for (int j = 0; j < Dimension; j++)//Dim=2
				Vect[p][i][j] = m[glodex].X[j] - m[p].X[j];
		}
	}
	
	/******1初始化模型信息完成******/
	/*部分文件设置*/
	//收敛记录文件
	Itfout.open("Iteration.dat");
	Itfout << "VARIABLES=" << "\"Iteration\"" << "\"N_Residuals\"" << "\"Rate_Residuals\"" << "\"SpecificR_Max\"" << "\n";
	//松弛记录文件
	if (RecordRelaxation == 1)
	{
		RelaxOut.open("Deformation_Iteration.dat");
		RelaxOut << "TITLE = \"Results of Deformation\"" << "\n";//file header
		RelaxOut << "FILETYPE = FULL" << "\n";
		RelaxOut << "variables= \"x\" \"y\" \"DisplacementX\" \"DisplacementY\" \"sigmaXX\" \"sigmaYY\" \"sigmaXY\" \"ResiX\" \"ResiY\" \"ReX\" \"ReY\"\n";
	}

	do
	{
		int Residuals = 0;
		/******2更新标记点信息******/
			// //2-1申请临时变量
			//申请随时释放的临时变量
			//临时变量-(nei0)临近点个数n_nei#######随时释放不用释放每次重新赋值
		int n_nei = 0;
		//临时变量-(nei1)指向矩阵V###########随时释放
		double** V = NULL;
		//临时变量-(nei2)临近点位移disp_nei#########随时释放
		double** disp_nei = NULL;
		//(nei3)临时临近点应力
		double*** sigma_nei = NULL;
		//申请最后释放的临时变量
		//临时变量-(m1)标记点位移dis_m!!!!!!!!!!!!!!最后释放
		double* disp_m = new double[Dimension];
		//临时变量-(m2)应变epsi!!!!!!!!最后释放
		double** epsi = new double* [Dimension];
		for (int i = 0; i < Dimension; i++)
			epsi[i] = new double[Dimension];
		//(m3)临时变量应力sigm[Dim][Dim]!!!!!!!!!最后释放
		double** sigm = new double* [Dimension];
		for (int i = 0; i < Dimension; i++)
			sigm[i] = new double[Dimension];
		//(m4)边界点外法向!!!!!!!!!最后释放
		double* n = new double[Dimension];
		//(m5)临时标记点内力
		double* inforce = new double[Dimension];

		//2-2-1更新全场应变&应力
		for (int p = 0; p < Num_m; p++)
		{
			/****2-1根据位移算应变****/
			/**2-1-1临时变量赋值**/
			//临时变量-(nei0)临近点个数n_nei
			n_nei = N_Nei[p];
			//临时变量-(nei1)指向矩阵V
			V = new double* [n_nei];
			for (int i = 0; i < n_nei; i++)
				V[i] = new double[Dimension];
			for (int i = 0; i < n_nei; i++)
				for (int j = 0; j < Dimension; j++)
					V[i][j] = Vect[p][i][j];
			//临时变量-(nei2)临近点位移disp_nei
			disp_nei = new double* [n_nei];
			for (int i = 0; i < n_nei; i++)
				disp_nei[i] = new double[Dimension];
			for (int i = 0; i < n_nei; i++)
			{
				int glodex = GlobalIndex_Nei[p][i];
				for (int d = 0; d < Dimension; d++)
					disp_nei[i][d] = m[glodex].displacement[d];
			}
			//临时变量-(m1)标记点位移dis_m
			for (int d = 0; d < Dimension; d++)
				disp_m[d] = m[p].displacement[d];
			//(m2)epsi[][]局外

			/**2-1-2计算应变存储在临时变量中**/
			Strain_from_Displacement(epsi, V, disp_nei, disp_m, n_nei, Dimension, DiffTolerance);

			/**2-1-3把应变传值给标记点**/
			for (int i = 0; i < Dimension; i++)//算完传值给m
				for (int j = 0; j < Dimension; j++)
					m[p].epsilon[i][j] = epsi[i][j];

			/**2-1-4释放需要随时释放的临时变量的内存**/
			for (int i = 0; i < n_nei; i++)
				delete[] V[i];
			delete[] V;//(nei1)释放方向矩阵V
			for (int i = 0; i < n_nei; i++)//(nei2)释放临时临近点位移disp_nei
				delete[] disp_nei[i];
			delete[] disp_nei;

			/****2-1根据位移算应变完成****/

			/****2-2应变算应力****/
			/**2-2-1准备临时变量赋值**/
			//(m2)epsi[][]还在
			//(m3)sigm[][]局外变量

			/**2-2-2本构方程计算应力**/
			Stress_from_Strain(sigm, E, nu, G, epsi, Dimension);

			/**2-2-3临时变量应力传给标记点**/
			for (int i = 0; i < Dimension; i++)//算完传值给m
				for (int j = 0; j < Dimension; j++)
					m[p].sigma[i][j] = sigm[i][j];
			/**2-2-4需要随时释放的内存**/

			/****2-2应变算应力完成****/
		}

		//2-2-2更新全场内力&合力
		for (int p = 0; p < Num_m; p++)
		{
			/****2-3应力算内力****/
			/**2-3-1临时变量赋值**/
			n_nei = N_Nei[p];//(nei0)临近点个数n_nei				
			V = new double* [n_nei];
			for (int i = 0; i < n_nei; i++)//(nei1)指向矩阵V
				V[i] = new double[Dimension];
			for (int i = 0; i < n_nei; i++)
				for (int j = 0; j < Dimension; j++)
					V[i][j] = Vect[p][i][j];
			sigma_nei = new double** [n_nei];//(nei3)临时临近点应力
			for (int n = 0; n < n_nei; n++)
				sigma_nei[n] = new double* [Dimension];
			for (int n = 0; n < n_nei; n++)
				for (int i = 0; i < Dimension; i++)
					sigma_nei[n][i] = new double[Dimension];
			for (int n = 0; n < n_nei; n++)
			{
				int glodex = GlobalIndex_Nei[p][n];
				for (int i = 0; i < Dimension; i++)
					for (int j = 0; j < Dimension; j++)
						sigma_nei[n][i][j] = m[glodex].sigma[i][j];
			}
			for (int i = 0; i < Dimension; i++)//(m3)-sigm[][]
				for (int j = 0; j < Dimension; j++)
					sigm[i][j] = m[p].sigma[i][j];
			for (int i = 0; i < Dimension; i++)//(m4)-n[]
				n[i] = m[p].n[i];

			/**2-3-2计算内力**/
			if ((m[p].IfBoundary) == 0)//内部点
			{
				InForce_from_Stress_in(inforce, V, sigma_nei, sigm, n_nei, Dimension, DiffTolerance);
			}
			else//外部点
			{
				InForce_from_Stress_bdy(inforce, sigm, n, Dimension);
			}
			/**2-3-3释放临时变量**/
			for (int n = 0; n < n_nei; n++)
				for (int i = 0; i < Dimension; i++)
					delete[] sigma_nei[n][i];
			for (int n = 0; n < n_nei; n++)
				delete[] sigma_nei[n];
			delete[] sigma_nei;//(nei3)释放临近点应力				
			for (int i = 0; i < n_nei; i++)
				delete[] V[i];
			delete[] V;//(nei1)释放方向矩阵V在2-1步2-3步用完了
			/****2-3应力算内力完成****/

			/****2-4令余量=内力+外力****/
			for (int i = 0; i < Dimension; i++)
			{
				if (m[p].Fixed[i] == 1)//固定余量为零
					m[p].R[i] = 0;
				else
					m[p].R[i] = inforce[i] + m[p].Force[i];
			}
			/****2-4令余量=内力完成****/
		}

		//2-3最后释放
		delete[] disp_m; //(m1)释放临时标记点位移dis_m
		for (int i = 0; i < Dimension; i++)//(m2)释放临时应变epsi[][]
			delete[] epsi[i];
		delete[] epsi;
		for (int i = 0; i < Dimension; i++)
			delete[] sigm[i];
		delete[] sigm;//(m3)释放标记点应力
		delete[] n;//(m4)释放边界点外法向
		delete[] inforce;//(m5)释放临时内力inforce[Dim]
		/******2更新标记点信息完成******/

		/******2+记录当前云图******/
		if (RecordRelaxation == 1)
		{
			RelaxOut << "zone T=\"Frame1\"\n";//zone header
			RelaxOut << "STRANDID=1\n";
			RelaxOut << "I=" << Num_m << "\n";
			RelaxOut << "ZONETYPE=Ordered\n";
			RelaxOut << "DATAPACKING=POINT\n";
			RelaxOut << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n";
			RelaxOut << "SOLUTIONTIME=" << it_orien << "\n";
			for (int p = 0; p < Num_m; p++)//data
			{
				RelaxOut << m[p].X[0] << " " << m[p].X[1] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << " ";
				if (m[p].IfBoundary == 1)
					RelaxOut << m[p].R[0] / ResiToleranceBounadry << " " << m[p].R[1] / ResiToleranceBounadry << "\n";
				else
					RelaxOut << m[p].R[0] / ResiToleranceInternal << " " << m[p].R[1] / ResiToleranceInternal << "\n";
			}
		}
		/******2+记录当前云图完成******/

		/******3检测余量数&无量纲余量分量最大值******/
		SpecificResidual_Max = 0.0;
		for (int p = 0; p < Num_m; p++)
			for (int i = 0; i < Dimension; i++)
			{
				//3-1检测余量数
				if (m[p].IfBoundary == 1)//区分边界点
					ResiTolerance = ResiToleranceBounadry;
				else
					ResiTolerance = ResiToleranceInternal;
				if (abs(m[p].R[i]) > ResiTolerance)
					Residuals++;
				//3-2SpecificResidual_Max
				if (abs(m[p].R[i]) / ResiTolerance > SpecificResidual_Max)
					SpecificResidual_Max = abs(m[p].R[i]) / ResiTolerance;

			}
		cout << "it_orien=" << it_orien << " Residuals=" << Residuals << " Rate_Resi=" << 100.0 * (double)Residuals / ((double)Num_m * (float)Dimension) <<"%" << "\n";
		Itfout << it_orien << " " << Residuals << " " << 100.0 * (double)Residuals / ((double)Num_m * (double)Dimension) <<" " <<SpecificResidual_Max<< "\n";
		/******4对每个Marker循环Move******/
		if (Residuals != 0)
		{
			for (int p = 0; p < Num_m; p++)
			{
				for (int i = 0; i < Dimension; i++)
				{
					if (m[p].Fixed[i] == 0)//非固定
					{
						if (m[p].IfBoundary == 0)//区分边界点
							m[p].displacement[i] += alpha * m[p].R[i];//内部点deltadisplacement alpha[p]=(eta/(E * (3.0 - nu) / (1 - nu * nu) / m[p].r / m[p].r / 2))
						else
							m[p].displacement[i] += beta * m[p].R[i];//边界点 beta[p]=(eta/(E / (1 - nu * nu) / m[p].r))
					}
					else;//固定点
				}
			}
			it_orien++;
		}
		else
			break;
		/******4对每个Marker循环Move完成******/
		
	} while (it_orien<= MaxIt_Rrien);
	Itfout.close();
	
	/******5输出结果******/
	//松弛记录关闭文件
	RelaxOut.close();
	/****5-1全场结果****/
	fout.open("SimulatedStress.dat");
	fout << "VARIABLES=" << "\"X\"" << "\"Y\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\"ResiX\"" << "\"ResiY\"" << "\"ReX\"" << "\"ReY\"" << "\n";
	for (int p = 0; p < Num_m; p++)
	{
		fout << m[p].X[0] << " " << m[p].X[1] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << " ";
		if(m[p].IfBoundary==1)
			fout << m[p].R[0]/ResiToleranceBounadry << " " << m[p].R[1] / ResiToleranceBounadry << "\n";
		else
			fout << m[p].R[0] / ResiToleranceInternal << " " << m[p].R[1] / ResiToleranceInternal << "\n";
	}
	fout.close();

	/****5-2纵截面Y=0结果****/
	if (OutputLongitudinalSection == 1)
	{
		LongSecfout.open("LongitudinalSection.dat");
		LongSecfout << "VARIABLES=" << "\"X\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\"ResiX\"" << "\"ResiY\"" << "\n";
		for (int p = 0; p < Num_m; p++)
		{
			if (abs(m[p].X[1]) <= 0.02)
				LongSecfout << m[p].X[0] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << "\n";
		}
		LongSecfout.close();
	}
	
	/****5-3横截面X=L/2结果****/
	if (OutputCrossSection == 1)
	{
		CrossSecfout.open("CrossSection.dat");
		CrossSecfout << "VARIABLES=" << "\"Y\"" << "\"dispX\"" << "\"dispY\"" << "\"sigmaxx\"" << "\"sigmayy\"" << "\"sigmaxy\"" << "\"ResiX\"" << "\"ResiY\"" << "\n";
		for (int p = 0; p < Num_m; p++)
		{
			if (abs(m[p].X[0] - L / 2) <= 0.02)
				CrossSecfout << m[p].X[1] << " " << m[p].displacement[0] << " " << m[p].displacement[1] << " " << m[p].sigma[0][0] << " " << m[p].sigma[1][1] << " " << m[p].sigma[0][1] << " " << m[p].R[0] << " " << m[p].R[1] << "\n";
		}
		CrossSecfout.close();
	}
	

	//释放全局变量内存
	delete[] m;
	
	for (int i = 0; i < Num_m; i++)
		delete[] GlobalIndex_Nei[i];
	delete[] GlobalIndex_Nei;

	for (int p = 0; p < Num_m; p++)
		for (int i = 0; i < N_Nei[p]; i++)
			delete[] Vect[p][i];
	for (int p = 0; p < Num_m; p++)
		delete[] Vect[p];
	delete[] Vect;

	delete[] N_Nei;
	LogFout.open("Log.dat");
	cout << "The run time is:" << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
	LogFout << "r=" << r << "\n";
	LogFout << "eta=" << eta << "\n";
	LogFout << "The run time is:" << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
	LogFout.close();
	return 0;
}
