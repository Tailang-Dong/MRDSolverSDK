#include "FileLayer.h"

std::ifstream mfin;//读取模型文件输入流
std::ofstream fout;//结果输出
std::ofstream LongSecfout;//纵截面结果输出Longitudinal section
std::ofstream CrossSecfout;//横截面结果输出Cross Section
std::ofstream Itfout;//迭代至收敛记录
std::ofstream RelaxOut;//松弛记录

std::ofstream LogFout;//日志记录

double SpecificResidual_Max=0;//记录所有标记点的无量纲余量各分量的最大值