#pragma once
#include <fstream>
//5文件层
//提前准备文件I/O流
//文件文件夹地址对象等等

extern std::ifstream mfin;//读取模型文件输入流
extern std::ofstream fout;//结果输出
extern std::ofstream LongSecfout;//纵截面结果输出Longitudinal Section
extern std::ofstream CrossSecfout;//横截面结果输出Cross Section
extern std::ofstream Itfout;//迭代至收敛记录
extern std::ofstream RelaxOut;//松弛记录

extern std::ofstream LogFout;//日志记录

extern double SpecificResidual_Max;//记录所有标记点的无量纲余量各分量的最大值

