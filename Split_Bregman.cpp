//#include<cmath>
//#include<stdlib.h>
//#include<opencv2/opencv.hpp>
//#include<opencv2/highgui.hpp>
//#include<iostream>
//#include<omp.h>
//#include<windows.h>
//#include<fstream>
//using namespace cv;
//using namespace std;
//
///* debug log */
///* 我之前写的Sidden算法有些bug，图像小时结果正确，图像大时结果错误。只因为自己在判定边界条件时用了一个自以为很好的代数骚操作
//*  结果当射线穿过网格交点时会出现重大问题，坑了自己。
//*/
//
///*函数功能及注意事项
//*	Sidden: 求解单角度下的系统子矩阵，每次用Sidden函数都会把把A置零，以方便使用，否则计算错误。注意，调用一次Sidden后A的改变值局部永久存在，到下一次调用Sidden时重置。
//*	 proj : 求解单角度下探测器板上的值矩阵，求解结果放在一个一维的临时容器p_temp里
//*	 Diff : 求解给定图像的X向差分和Y向差分，对照于差分算子D(),统一采用向后差分规范,差分矩阵的目标对应形式实为一维向量
//*	Diff_T: 求解给定图像的x向转置差分和Y向转置差分，这个和Diff大同小异，是差分矩阵次对角线的对称位置相移，结果由向后差分变为了向前差分的负值，只是分块矩阵中每一块的首尾部情况特殊
//*	f_solve:求解目标函数，f = arg{min{|dx|_1+|dy|_1+(mu/2)*pow(|Au-p|_2,2)+(lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)}}
//*  u_update:使用最速下降法进行数据更新，函数结果为求解更新量,调用一次后结果存储在u_temp中
//*	shrink: 更新dx,dy,利用软阈值算法。使用此函数后Dx,Dy依赖于u。
//*/
//
//void Sidden(double** A, const int view, const double da, double point_source_x, double point_source_y, double* bar_x, double* bar_y, const int pNum, const int dNum, const double pSize, const double dSize, double** imgX, double** imgY)
//{
//	for (int i = 0; i < dNum; ++i) { // 置零A
//		for (int j = 0; j < pNum * pNum; ++j) {
//			A[i][j] = 0;
//		}
//	}
//	// 定位单角度下的射线簇
//	double phi{ view * da };
//	double rot_point_source_x{ 0 }; double rot_point_source_y{ 0 };
//	rot_point_source_x = cos(phi) * point_source_x - sin(phi) * point_source_y;
//	rot_point_source_y = sin(phi) * point_source_x + cos(phi) * point_source_y;
//
//	double* rot_bar_x = new double[dNum] {}; double* rot_bar_y = new double [dNum] {};
//	for (int i = 0; i < dNum; ++i) {
//		rot_bar_x[i] = cos(phi) * bar_x[i] - sin(phi) * bar_y[i];
//		rot_bar_y[i] = sin(phi) * bar_x[i] + cos(phi) * bar_y[i];
//	}
//#pragma omp parallel for
//	for (int k = 0; k < dNum; ++k) { // 单角度单射线的考虑，可并行
//		double slope = (rot_bar_y[k] - rot_point_source_y) / (rot_bar_x[k] - rot_point_source_x + 1e-13);
//
//		if (abs(slope) < 1) { // X主序
//			double length_through_pixel = pSize * sqrt(1 + pow(slope, 2)); // 两个layer之间的固定长
//
//			/*接下来对特殊情况进行处理
//			* 1. 首先是射线不满穿网格的判定
//			* 2. 接着是单layer跨越非主序的判定
//			*/
//			double* unmain_sequence_value = new double[pNum + 1]{}; // 图像模型构建是用的像素中心，这里使用像素网格边框
//			for (int i = 0; i < pNum + 1; ++i) {
//				if (i < pNum) {
//					unmain_sequence_value[i] = slope * ((imgX[0][i] - pSize / 2.0) - rot_point_source_x) + rot_point_source_y;
//				}
//				else {
//					unmain_sequence_value[i] = slope * ((imgX[0][i - 1] + pSize / 2.0) - rot_point_source_x) + rot_point_source_y;
//				}
//			}
//			// 先求出unmain_sequence_value对应的Y序号，方便下面程序的使用
//			double* unmain_sequence_number = new double[pNum + 1]{};
//			for (int i = 0; i < pNum + 1; ++i) {
//				unmain_sequence_number[i] = (int)floor(-(unmain_sequence_value[i] - pNum * pSize / 2.0) / pSize);
//			}// 这时候unmain_sequence_number可能有负值，需要留意使用。这里也不能使用int，因为有负值，使用的话会计算错
//
//			for (int i = 0; i < pNum; ++i) { // 判断穿越像素的长度
//				if (unmain_sequence_number[i] == unmain_sequence_number[i + 1]) {
//					// 这时候射线经过每个layer穿过一个像素 // 重构时考虑到了边界条件
//					// 易错误的判断条件：ceil(unmain_sequence_value[i] / pSize) == ceil(unmain_sequence_value[i + 1] / pSize)这种写法不对，对临界情况处理不好
//					// 
//					// 这时候有两种可能：上述线段对应像素在网格内，给予计算赋值；上述线段对应像素不在网格内，赋值0，因为初始化为0，故可不再做处理。
//					if (abs(ceil(unmain_sequence_value[i] / pSize)) <= pNum / 2.0 && abs(floor(unmain_sequence_value[i] / pSize)) <= pNum / 2.0) {
//						// 这里判断条件写成< (unmain_sequence_number[i]>=0) && (unmain_sequence_number[i]<=pNum) >不太好，因为没考虑好边界条件，会多放进来值
//						A[k][((int)unmain_sequence_number[i]) * pNum + i] = length_through_pixel; // 这里数组里面必须加int强制类型转换，否则报错，无奈之举
//						// 第一次这里报bug,因为floor里面的前面那部分少写了个负号
//					}// 这里不加(int)会报错,过不了编译
//				}
//				else { // 这时候射线经过每个layer穿过两个像素，而且两个像素不一定都在网格里面。
//					// 若每个像素都在网格里模型又可分为两种情况：斜率为[0,1)；斜率为(-1,0]；
//					// 特殊的临界情况在这里不存在，因为X主序且前面的if判断了
//					// 这时候还可以用代数方式直接计算，一个规律是（floor前一个+ceil后一个）/2=中间的那个Y序号
//					// 此后利用相似三角形即可计算，相似三角形采取小的与大的比较，而不是对角的那两个三角形
//					int  temp_Y{ 0 }; // 穿过一个layer时一半在网格外，一半在网格内，利用这个Y序号容易比较判断
//					temp_Y = (int)((ceil(unmain_sequence_value[i] / pSize) + floor(unmain_sequence_value[i + 1] / pSize)) / 2.0);
//					// 这里计算的temp_Y是与主序相对应的一个Y坐标，是夹在两个layer的Y之间的网格线
//					if (temp_Y == -pNum / 2) { // 边界情形
//						if (unmain_sequence_value[i] < -pNum * pSize / 2.0) {
//							A[k][(pNum - 1) * pNum + i] = abs(unmain_sequence_value[i + 1] + pNum * pSize / 2.0) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//						else {
//							A[k][(pNum - 1) * pNum + i] = abs(unmain_sequence_value[i] + pNum * pSize / 2.0) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//					}
//					else if (temp_Y == pNum / 2) {
//						if (unmain_sequence_value[i] > pNum * pSize / 2.0) {
//							A[k][i] = abs((pNum * pSize / 2.0) - unmain_sequence_value[i + 1]) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//						else {
//							A[k][i] = abs((pNum * pSize / 2.0) - unmain_sequence_value[i]) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//					}
//					else if (temp_Y<(pNum / 2) && temp_Y>(-pNum / 2)) { // 内部情形
//						// 下面这是重构后的，因为重构后发现根本就不需要写这些判断，因为都一样，那Sidden算法编程更简单了，有意思
//						A[k][((int)unmain_sequence_number[i]) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						A[k][((int)unmain_sequence_number[i + 1]) * pNum + i] = length_through_pixel - A[k][((int)unmain_sequence_number[i]) * pNum + i];
//						//if (abs(unmain_sequence_value[i] - pSize * pNum / 2.0) < -(temp_Y - pNum / 2.0) * pSize) { // 一个layer两段，左段在上面那个像素里面
//						//	A[k][((int)unmain_sequence_number[i]) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	A[k][((int)unmain_sequence_number[i + 1]) * pNum + i] = length_through_pixel - A[k][((int)unmain_sequence_number[i]) * pNum + i];
//						//	//A[k][(-(temp_Y - pNum / 2) - 1) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	//A[k][-(temp_Y - pNum / 2) * pNum + i] = length_through_pixel - A[k][(-(temp_Y - pNum / 2) - 1) * pNum + i];
//						//}// 能跑的代码实际上没必要就不改，这里实际上用unmain_sequence_number更明确
//						//else { // 一个layer两段，左段在下面那个像素里面
//						//	A[k][((int)unmain_sequence_number[i]) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	A[k][((int)unmain_sequence_number[i + 1]) * pNum + i] = length_through_pixel - A[k][((int)unmain_sequence_number[i]) * pNum + i];
//						//	//A[k][-(temp_Y - pNum / 2) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	//A[k][(-(temp_Y - pNum / 2) - 1) * pNum + i] = length_through_pixel - A[k][-(temp_Y - pNum / 2) * pNum + i];
//						//}
//					}
//					// 外部情形就不用考虑了，把该赋值的赋值了，剩下的就是没被射线穿越的像素，而且初始化时就置为了0，所以逻辑恰当
//				}
//			}
//			delete[] unmain_sequence_value; delete[] unmain_sequence_number;
//		}
//
//		else { // Y主序
//			// 同X主序的分析步骤，序号之类的做稍微变动以对应
//			double length_through_pixel = pSize * sqrt(1 + pow(slope, -2)); // 注意这里pow函数里面是-2
//
//			double* unmain_sequence_value = new double[pNum + 1]{}; // 图像模型构建是用的像素中心，这里使用像素网格边框
//			for (int i = 0; i < pNum + 1; ++i) {
//				if (i < pNum) {
//					unmain_sequence_value[i] = ((imgY[i][0] + pSize / 2.0) - rot_point_source_y) / (slope + 1e-13) + rot_point_source_x;
//				}
//				else {
//					unmain_sequence_value[i] = ((imgY[i - 1][0] - pSize / 2.0) - rot_point_source_y) / (slope + 1e-13) + rot_point_source_x;
//				}
//			}
//			double* unmain_sequence_number = new double[pNum + 1]{};
//			for (int i = 0; i < pNum + 1; ++i) {
//				unmain_sequence_number[i] = (int)floor((unmain_sequence_value[i] + pNum * pSize / 2.0) / pSize);
//			}
//
//			for (int i = 0; i < pNum; ++i) {
//				if (unmain_sequence_number[i] == unmain_sequence_number[i + 1]) {
//					if (abs(ceil(unmain_sequence_value[i] / pSize)) <= pNum / 2.0 && abs(floor(unmain_sequence_value[i] / pSize)) <= pNum / 2.0) {
//						A[k][i * pNum + (int)unmain_sequence_number[i]] = length_through_pixel;
//					}
//				}
//				else {
//					int  temp_X{ 0 };
//					temp_X = (int)((ceil(unmain_sequence_value[i] / pSize) + floor(unmain_sequence_value[i + 1] / pSize)) / 2.0);
//					if (temp_X == -pNum / 2) { // 边界情形
//						if (unmain_sequence_value[i] < -pNum * pSize / 2.0) {
//							A[k][i * pNum] = abs(unmain_sequence_value[i + 1] + pNum * pSize / 2.0) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//						else {
//							A[k][i * pNum] = abs(unmain_sequence_value[i] + pNum * pSize / 2.0) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//					}
//					else if (temp_X == pNum / 2) {
//						if (unmain_sequence_value[i] > pNum * pSize / 2.0) {
//							A[k][(i + 1) * pNum - 1] = abs((pNum * pSize / 2.0) - unmain_sequence_value[i + 1]) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//						else {
//							A[k][(i + 1) * pNum - 1] = abs((pNum * pSize / 2.0) - unmain_sequence_value[i]) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						}
//					}
//					else if (temp_X<(pNum / 2) && temp_X>(-pNum / 2)) {
//						A[k][i * pNum + (int)unmain_sequence_number[i]] = abs(abs(unmain_sequence_value[i] + pSize * pNum / 2.0) - (temp_X + pNum / 2.0) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						A[k][i * pNum + (int)unmain_sequence_number[i + 1]] = length_through_pixel - A[k][i * pNum + (int)unmain_sequence_number[i]];
//						//if (abs(unmain_sequence_value[i] + pSize * pNum / 2.0) < (temp_X + pNum / 2.0) * pSize) {
//						//	A[k][i * pNum - 1 + (temp_X + pNum / 2)] = abs(abs(unmain_sequence_value[i] + pSize * pNum / 2.0) - (temp_X + pNum / 2.0) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	A[k][i * pNum + (temp_X + pNum / 2)] = length_through_pixel - A[k][i * pNum - 1 + (temp_X + pNum / 2)];
//						//}
//						//else {
//						//	A[k][i * pNum + (temp_X + pNum / 2)] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (temp_X + pNum / 2.0) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	A[k][i * pNum - 1 + (temp_X + pNum / 2)] = length_through_pixel - A[k][i * pNum + (temp_X + pNum / 2)];
//						//}
//					}
//				}
//			}
//			delete[] unmain_sequence_value; delete[] unmain_sequence_number;
//		}
//	}
//	delete[] rot_bar_x; delete[] rot_bar_y;
//};
//
//void proj(double* p_temp, double** A, const int view, const double da, double point_source_x, double point_source_y, double* bar_x, double* bar_y, const int pNum, const int dNum, const double pSize, const double dSize, double** imgX, double** imgY, double** u)
//{
//	for (int i = 0; i < dNum; ++i) {
//		p_temp[i] = 0; // 这里第一次写忘了置零了
//	}// 但尴尬的是，第一次修改后将dNum写成了pNum，还是报错了，排查了好久才查出来
//	Sidden(A, view, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY);
//	for (int i = 0; i < dNum; ++i) {
//		for (int j = 0; j < pNum * pNum; ++j) {
//			p_temp[i] += A[i][j] * u[j / pNum][j % pNum];
//		}
//	}
//}
//
//void Diff(double** u, double** Dx, double** Dy, int pNum)
//{
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			if (j == 0) Dx[i][j] = 0;
//			else Dx[i][j] = u[i][j] - u[i][j - 1];
//
//			if (i == 0) Dy[i][j] = 0;
//			else Dy[i][j] = u[i][j] - u[i - 1][j];
//		}
//	}
//}
//
//void Diff_T(double** u, double** Dx_T, double** Dy_T, int pNum)
//{
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			if (j == 0) Dx_T[i][j] = -u[i][j + 1];
//			else if (j == (pNum - 1)) Dx_T[i][j] = u[i][j]; // 第一次在这里出现了一个难以排查的bug
//			else Dx_T[i][j] = u[i][j] - u[i][j + 1];
//
//			if (i == 0) Dy_T[i][j] = -u[i + 1][j]; // 逐行调试后发现是这里的Dy_T无法读取内存
//			else if (i == (pNum - 1)) Dy_T[i][j] = u[i][j];
//			else Dy_T[i][j] = u[i][j] - u[i + 1][j];
//			// std::cout << "成功驶出" << i << std::endl; 
//		}
//	}
//}
//
//double f_solve(double** dx, double** dy, double** bx, double** by, double mu, double lamda, double** u, double** Dx, double** Dy, \
//	double* p_temp, double** A, const int view, const double da, double point_source_x, double point_source_y, double* bar_x, double* bar_y, const int pNum, const int dNum, const double pSize, const double dSize, double** imgX, double** imgY, double** img, double** p)
//{
//	// f = arg{min{|dx|_1+|dy|_1+(mu/2)*pow(|Au-p|_2,2)+(lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)}}
//	double f_1{ 0 }, f_2{ 0 }, f_3{ 0 }; // 将f_solve的求解分为三部分
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			f_1 += abs(dx[i][j]) + abs(dy[i][j]); // |dx|_1+|dy|_1
//		}
//	}
//
//	proj(p_temp, A, view, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, u); // 当前解u下的情形
//	for (int i = 0; i < dNum; ++i) {
//		f_2 += pow((p_temp[i] - p[i][view]), 2);
//	}f_2 = mu * 0.5 * f_2; // (mu/2)*pow(|Au-p|_2,2)
//
//	Diff(u, Dx, Dy, pNum);
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			f_3 += pow((dx[i][j] - Dx[i][j] - bx[i][j]), 2) + pow((dy[i][j] - Dy[i][j] - by[i][j]), 2);
//		}
//	}f_3 = lamda * 0.5 * f_3; // (lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)
//
//	return (f_1 + f_2 + f_3);
//}
//
//void u_update(double** u_temp, double** Dx_T, double** Dy_T, double** dx, double** dy, double** bx, double** by, double mu, double lamda, double** u, double** Dx, double** Dy, \
//	double* p_temp, double** A, const int view, const double da, double point_source_x, double point_source_y, double* bar_x, double* bar_y, const int pNum, const int dNum, const double pSize, const double dSize, double** imgX, double** imgY, double** img, double** p)
//{
//	/*求解公式
//	*	d^(k)=mu*A'(Au^(k)-p)+lamda*Dx'(Dx(u^(k))+bx^(k)-dx)+lamda*Dy'(Dy(u^(k))+by^(k)-dy)
//	*	alphi^(k) = (d^(k)'*d^(k))/(d^(k)'*H*d^(k)),其中H为Hessian矩阵
//	*	H = lamda*Dx'Dx + lamda*Dy'Dy+mu*A'A
//	*-->d^(k)'*H*d^(k) = lamda*pow(|Dx(d^(k))|_2,2)+  lamda*pow(|Dy(d^(k))|_2,2) + mu*pow(|Ad^(k)|_2,2)
//	*/
//	proj(p_temp, A, view, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, u);
//	double** d_1 = new double* [pNum];
//	double** d_2 = new double* [pNum];
//	for (int i = 0; i < pNum; ++i) {
//		d_1[i] = new double[pNum] {};
//		d_2[i] = new double[pNum] {};
//	}
//	for (int i = 0; i < dNum; ++i) {
//		p_temp[i] = p_temp[i] - p[i][view];
//	}
//	for (int i = 0; i < pNum * pNum; ++i) {
//		for (int j = 0; j < dNum; ++j) {
//			d_1[i / pNum][i % pNum] += A[j][i] * p_temp[j]; // A'(Au^(k)-p)
//		}
//	}
//	// 这里第一次写的时候有一些些的问题，此后一部分都有点问题
//	Diff(u, Dx, Dy, pNum);
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			d_2[i][j] = Dx[i][j] + bx[i][j] - dx[i][j];// Dx(u^(k))+bx^(k)-dx
//		}
//	}
//
//	Diff_T(d_2, Dx_T, Dy_T, pNum); // Dx'(Dx(u^(k))+bx^(k)-dx)
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			u_temp[i][j] = mu * d_1[i][j] + lamda * Dx_T[i][j]; // mu* A'(Au^(k)-p)+lamda*Dx'(Dx(u ^ (k)) + bx ^ (k)-dx)
//		}
//	}
//
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			d_2[i][j] = Dy[i][j] + by[i][j] - dy[i][j];// Dy(u^(k))+by^(k)-dy
//		}// 这里前面<Diff(u, Dx, Dy, pNum);>调用产生过Dy了
//	}
//	Diff_T(d_2, Dx_T, Dy_T, pNum); // 这里d_2值变化，需要重新生成对应的Dy_T
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			u_temp[i][j] += lamda * Dy_T[i][j]; // lamda*Dy'(Dy(u^(k))+by^(k)-dy)
//		}
//	}
//	// 下面求解最优步长alphi
//	Diff(u_temp, Dx, Dy, pNum); // Dx(d^(k)),Dy(d^(k))
//	proj(p_temp, A, view, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, u_temp);//Ad ^ (k)
//	double alphi{ 0 }, alphi_1{ 0 }, alphi_2{ 0 }, alphi_3{ 0 };
//	// alphi ^ (k) = (d ^ (k)'*d^(k))/(d^(k)' * H * d ^ (k)), 其中H为Hessian矩阵
//	// d^(k)'*H*d^(k) = lamda*pow(|Dx(d^(k))|_2,2)+  lamda*pow(|Dy(d^(k))|_2,2) + mu*pow(|Ad^(k)|_2,2)
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			alphi_1 += pow(u_temp[i][j], 2); // (d ^ (k)'*d^(k))
//			alphi_2 += (pow(Dx[i][j], 2) + pow(Dy[i][j], 2));// pow(|Dx(d^(k))|_2,2),pow(|Dy(d^(k))|_2,2)
//		}
//	}
//	for (int i = 0; i < dNum; ++i) {
//		alphi_3 += pow(p_temp[i], 2); // pow(|Ad^(k)|_2,2)
//	}
//	alphi = alphi_1 / (lamda * alphi_2 + mu * alphi_3);
//
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			u_temp[i][j] = alphi * u_temp[i][j]; // 步长更新量
//		}
//	}
//
//	for (int i = 0; i < pNum; ++i) {
//		delete[] d_1[i];
//		delete[] d_2[i];
//	}delete[] d_1; delete[] d_2;
//}
//
//void shrink(double** dx, double** dy, double** u, double** Dx, double** Dy, double** bx, double** by, int pNum, double lamda)
//{
//	// dx=shrink(Dx(u)+bx,1/lamda),dy=shrink(Dy(u)+by,1/lamda)
//	Diff(u, Dx, Dy, pNum);
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			if (Dx[i][j] + bx[i][j] > 1 / lamda) dx[i][j] = Dx[i][j] + bx[i][j] - 1 / lamda;
//			else if (Dx[i][j] + bx[i][j] < -1 / lamda) dx[i][j] = Dx[i][j] + bx[i][j] + 1 / lamda;
//			else dx[i][j] = 0;
//
//			if (Dy[i][j] + by[i][j] > 1 / lamda) dy[i][j] = Dy[i][j] + by[i][j] - 1 / lamda;
//			else if (Dy[i][j] + by[i][j] < -1 / lamda) dy[i][j] = Dy[i][j] + by[i][j] + 1 / lamda;
//			else dy[i][j] = 0;
//		}
//	}
//}
//
//int main()
//{
//	/*坐标系建立，确定视野圆，再生成图像网格，最后定义图像模体*/
//	constexpr int sod{ 600 }; constexpr int sdd{ 1000 }; constexpr int dNum{ 512 }; constexpr double dSize{ 2 }; constexpr int pNum{ 256 };
//	// 这里第一次写的时候dSize写成0.2结果不报错，写成2就报错
//	constexpr int views{ 60 }; constexpr double pi{ 3.1415926927 };
//	constexpr double da{ 2 * pi / views };
//
//	constexpr double L = (dNum - 1) * dSize / 2;
//	double R = L * sod / sqrt(pow(L, 2) + pow(sdd, 2));
//	double pSize = 2 * R / pNum; // 确定了视野圆和网格大小，接下来要存储一些必要的坐标：定射线、定网格
//
//	double point_source_x{ 0 }, point_source_y{ sod }; // 定射线
//	double* bar_x = new double[dNum];
//	double* bar_y = new double[dNum];
//	for (int i = 0; i < dNum; ++i) {
//		bar_x[i] = -L + i * dSize;
//		bar_y[i] = (double)sod - (double)sdd;
//	}
//	double** imgX = new double* [pNum]; // 定网格
//	double** imgY = new double* [pNum];
//	for (int i = 0; i < pNum; ++i) {
//		imgX[i] = new  double[pNum];
//		imgY[i] = new double[pNum];
//	}
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			imgX[i][j] = -R + pSize / 2 + pSize * j; // 注意这里乘的是j
//			imgY[i][j] = R - pSize / 2 - pSize * i;  // 数组的序号和网格坐标系是不一致的
//		}
//	}
//
//	double** img = new double* [pNum]; // 定义图像模体
//	for (int i = 0; i < pNum; ++i) {
//		img[i] = new double[pNum];
//	}
//	int x_1 = 0, y_1 = 0, r_1 = 200, density_1 = 50; // 这里设置为double后面使用Mat会出问题，主要是我对Mat不熟悉
//	int x_2 = 100, y_2 = 100, r_2 = 50, density_2 = 100;
//	for (int i = 0; i < pNum; ++i) {// 灰度容器的每个值对应一个像素中心点
//		for (int j = 0; j < pNum; ++j) {
//			if ((pow(imgX[i][j] - x_1, 2) + pow(imgY[i][j] - y_1, 2)) < pow(r_1, 2)) {
//				img[i][j] = density_1;
//			}
//			else img[i][j] = 0; // 这里逻辑要处理好，不仅要图像都显示出来，还要图层合理 
//								// 这块相当于数组初始化了，图层一
//
//			if ((pow(imgX[i][j] - x_2, 2) + pow(imgY[i][j] - y_2, 2)) < pow(r_2, 2)) {
//				img[i][j] = density_2;
//			}
//			if ((pow(imgX[i][j], 2) + pow(imgY[i][j], 2)) > pow(R - pSize / 2, 2)) // 视野圆外的置零
//			{
//				img[i][j] = 0;
//			}
//		}
//	}
//
//	/*借助OpenCv来显示图像*/
//	Mat BB = Mat(pNum, pNum, CV_8UC1);
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			BB.data[i * pNum + j] = (int)img[i][j];
//		}
//	}
//	imwrite("C:\\Users\\team8\\source\\repos\\cplusplus/Source.jpg", BB);
//	//imshow("source image", BB);
//	//waitKey();
//
//
//	/*利用Sidden算法计算系统矩阵
//	* 1. 存储空间大小设定为单角度下的系数矩阵所需，要不然全角度数量级为n^4，吃不消
//	* 2. 构造函数Sidden，通过传递关键参数view来确定主序及矩阵生成
//	*/
//	double** A = new double* [dNum];
//	for (int i = 0; i < dNum; ++i) {
//		A[i] = new double[pNum * pNum]{}; // 系统矩阵及初始化
//	} // 这里设计一个求Sidden算法的函数
//
//	/*Split_Bregman算法
//	* 1. 和ADMM方法没啥区别，不过是改变了乘子更新方式
//	* 2. 交替更新时选择最速下降法求解而不用ppt里面的解析解，那个没学过，保守起见现在不用
//	* 3. 初始图像直接选择全零矩阵，先不使用OS-SART迭代出来的图像做基础图像，看看速度
//	* 4. 只使用ADMM方法更新试试，一切先从最简单的组合开始
//	* 5. 使用各向异性模型
//	* 6. 这种情况可能陷入我们不想要的局部最优解，本质上初始图像选择OS-SART迭代后的数据时结果也可能是局部最优（解不稀疏时），但是我们想要的
//	*/
//
//	/*约束模型
//	* f = arg{min{|dx|_1+|dy|_1+(mu/2)*pow(|Au-p|_2,2)+(lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)}}
//	* 其中 ：f是目标函数，Dx是差分算子,bx是Split_Bregman方法的特色（替代lamda的渐大更新）
//	*/
//
//	/*求解公式
//	* 1. dx,dy给定初始值0，并且每次迭代完u后按软阈值方式更新dx=shrink(Dx(u)+bx,1/lamda),dy=shrink(Dy(u)+by,1/lamuda)
//	* 2. 最速下降法更新模型u^(K+1) = u^(k) - alpha^(k)*d^(k)，其中dk为梯度方向，alphi为最优步长
//	*	d^(k)=mu*A'(Au^(k)-p)+lamda*Dx'(Dx(u^(k))+bx^(k)-dx)+lamda*Dy'(Dy(u^(k))+by^(k)-dy)
//	*	alphi^(k) = (d^(k)'*d^(k))/(d^(k)'*H*d^(k)),其中H为Hessian矩阵
//	*	H = lamda*Dx'Dx + lamda*Dy'Dy+mu*A'A
//	*-->d^(k)'*H*d^(k) = lamda*pow(|Dx(d^(k))|_2,2)+  lamda*pow(|Dy(d^(k))|_2,2) + mu*pow(|Ad^(k)|_2,2)
//	* 3. Split_Bregman方法参数bx的更新（不走回头路）：bx^(k+1) = bx^(k) + Dx(u)^(k+1) - dx^(k+1)
//	* 4. 罚参数的更新： mu = 2*mu,初始值设置为0.05，最大值设置为pow(2,13)
//	*/
//
//	// 约束模型的求解
//	double** u = new double* [pNum]; double** dx = new double* [pNum]; double** dy = new double* [pNum]; double** bx = new double* [pNum]; double** by = new double* [pNum];
//	for (int i = 0; i < pNum; ++i) {
//		u[i] = new double[pNum] {}; dx[i] = new double[pNum] {}; dy[i] = new double[pNum] {}; bx[i] = new double[pNum] {}; by[i] = new double[pNum] {};
//	}
//
//	constexpr int maxIte_out{ 1000 }; constexpr int maxIte_inn{ 100 };
//	double mu = 0.05; double lamda = pow(2, 5); // 与迭代有关的参数设定
//
//	double** p = new double* [dNum];
//	for (int i = 0; i < dNum; ++i) {
//		p[i] = new double[views] {};
//	} // 这里设计一个实现单角度正投影的函数
//	double* p_temp = new double[dNum] {}; // 单角度下投影矩阵的临时存储容器
//	for (int i = 0; i < views; ++i) { // 这里不能用openmp，因为A和p_temp定义在外面只有一份
//		proj(p_temp, A, i, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img);
//		for (int j = 0; j < dNum; ++j) {
//			p[j][i] = p_temp[j];
//		}
//	} // 计算得到了全角度真实正投影数据，以做重建质量评价
//
//	/*借助OpenCv来显示图像*/
//	// 这里投影值的范围你不知道，所以你需要压缩显示
//	double proj_max{ p[0][0] };
//	for (int i = 0; i < dNum; ++i) {
//		for (int j = 0; j < views; ++j) {
//			if (proj_max < p[i][j]) {
//				proj_max = p[i][j];
//			}
//			else continue;
//		}
//	}
//	Mat	CC = Mat(dNum, views, CV_8UC1);
//	for (int i = 0; i < dNum; ++i) {
//		for (int j = 0; j < views; ++j) {
//
//			CC.data[i * views + j] = (int)((p[i][j]) / proj_max * 255);
//		}
//	}
//	imwrite("C:\\Users\\team8\\source\\repos\\cplusplus/Projection.jpg", CC);
//	//imshow("projection image", CC);
//	//waitKey();
//
//
//	double** Dx = new double* [pNum]; double** Dy = new double* [pNum];// 这里设计一个函数来实现差分功能
//	double** Dx_T = new double* [pNum]; double** Dy_T = new double* [pNum]; // 这里设计一个函数实现差分算子转置功能
//	for (int i = 0; i < pNum; ++i) {
//		Dx[i] = new double[pNum] {};
//		Dy[i] = new double[pNum] {};
//		Dx_T[i] = new double[pNum] {};
//		Dy_T[i] = new double[pNum] {}; // 第一次我这里把Dy_T写成了Dx_T，结果排查了半天，真的搞人
//	}
//
//	double f_current{ 0 }, f_previous{ 0 }; // 这里设计一个函数来计算目标函数，一个表示迭代后的解，一个表示迭代前的解
//
//	double** u_temp = new double* [pNum]; // 定义一个中间量来存储临时结果图像或者沿梯度更新的数据
//	for (int i = 0; i < pNum; ++i) {
//		u_temp[i] = new double[pNum] {};
//	} // 这里设计一个函数来计算沿梯度更新量
//	double error{ 0 }, error_1{ 0 }, error_2{ 0 }; // 用于判断程序终止条件
//
//	// 迭代开始
//	ULONGLONG t1, t2;
//	t1 = GetTickCount64();
//	for (int ite_out = 0; ite_out < maxIte_out; ++ite_out) { // 设定总共迭代多少次
//		// 先生成随机数,通过随机交换数组中的两个数实现随机效果
//		int ran_views[views]{ 0 };
//		for (int i = 0; i < views; ++i) {
//			ran_views[i] = i;
//		}
//		for (int i = 0; i < views; ++i) {
//			int temp = 0;
//			int j = rand() % views;
//			temp = ran_views[i];
//			ran_views[i] = ran_views[j];
//			ran_views[j] = temp;
//		}
//		for (int k = 0; k < views; ++k) { // 逐角度迭代，每个角度使用Splite_Bregman方法100次
//			// 这里选角度要随机，速度会大大提高
//			k = ran_views[k];
//			for (int ite = 0; ite < maxIte_inn; ++ite) {
//				f_current = f_solve(dx, dy, bx, by, mu, lamda, u, Dx, Dy, \
//					p_temp, A, views, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//
//				u_update(u_temp, Dx_T, Dy_T, dx, dy, bx, by, mu, lamda, u, Dx, Dy, \
//					p_temp, A, k, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						u_temp[i][j] = u[i][j] - u_temp[i][j]; // 此时u_temp代表沿梯度更新后的图像
//						if (u_temp[i][j] < 0) u_temp[i][j] = 0; // 非负约束
//					}
//				}
//				f_previous = f_current;
//				f_current = f_solve(dx, dy, bx, by, mu, lamda, u_temp, Dx, Dy, \
//					p_temp, A, views, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//
//				if (f_current >= f_previous) { // 加一个人工判定，使每次循环尽量有效
//					for (int t = 0; t < 5; t++) { // 若不满足，则至少人工优化五次
//						for (int i = 0; i < pNum; ++i) {
//							for (int j = 0; j < pNum; ++j) {
//								u_temp[i][j] = (u[i][j] + u_temp[i][j]) / 2.0;
//							}
//						}
//						f_current = f_solve(dx, dy, bx, by, mu, lamda, u_temp, Dx, Dy, \
//							p_temp, A, views, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//						if (f_current < f_previous) break;
//					}
//				}
//
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						error_1 += pow((u_temp[i][j] - u[i][j]), 2) / 255; // 防止累加超了double类型准许长度
//						error_2 += pow(u[i][j], 2) / 255;
//					}
//				}error = error_1 / error_2;
//
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						u[i][j] = u_temp[i][j];
//					}
//				}
//				/***************至此u更新完毕*,接下来更新dx,dy，bx,by***********/
//					// dx=shrink(Dx(u)+bx,1/lamda),dy=shrink(Dy(u)+by,1/lamda)
//					// bx^(k+1) = bx^(k) + Dx(u)^(k+1) - dx^(k+1)
//					// 需要写个关于dx,dy的软阈值函数
//				shrink(dx, dy, u, Dx, Dy, bx, by, pNum, lamda);
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						bx[i][j] = bx[i][j] + Dx[i][j] - dx[i][j];
//						by[i][j] = by[i][j] + Dy[i][j] - dy[i][j]; // 这里直接利用shrink函数改变的Dx,Dy的值
//					}
//				}
//				/***************至此dx,dy,bx,by更新完毕*,接下来更新罚参数mu***********/
//				//mu = pow(2, 13);
//				if (mu < pow(2, 13)) mu = 2 * mu;
//
//				std::cout << "ite: " << ite << std::endl;
//			}
//			std::cout << "views: " << k << std::endl;
//
//			//Mat GG = Mat(pNum, pNum, CV_8UC1);
//			//for (int i = 0; i < pNum; ++i) 	{
//			//	for (int j = 0; j < pNum; ++j) 
//			//	{
//			//		GG.data[i * pNum + j] = (int)u[i][j];
//			//	}
//			//}
//			//imshow("u_ image", GG);
//			//waitKey();
//		}
//
//		std::cout << "ite_out: " << ite_out << std::endl;
//		if (error < 1e-5) break;
//	}
//	t2 = GetTickCount64();
//	cout << "Split_Bregman法迭代耗时： " << (t2 - t1) / (1000 * 60.0) << "min" << endl;
//
//	ofstream fout("Split_Bregman.txt", ios_base::out | ios_base::app);
//	fout << "Split_Bregman法迭代耗时： " << (t2 - t1) / (1000 * 60.0) << "min" << endl;
//	fout.close();
//
//	// 显示图像时矩阵是double类型，需要注意
//
//		/*借助OpenCv来显示图像*/
//	Mat EE = Mat(pNum, pNum, CV_8UC1);
//	// 若是double类型img设置为cv_64FC1会出现类型转换错误，但我不了解Mat,所以不清楚这是为什么
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) // 之前这里把j写成i了，真是能粘贴的我绝不再写第二遍！！！
//		{
//			EE.data[i * pNum + j] = (int)u[i][j];
//		}
//	}
//	imwrite("C:\\Users\\team8\\source\\repos\\cplusplus/reconstruction.jpg", EE);
//	//imshow("reconstruction image", EE);
//	//waitKey();
//
//	Mat DD = Mat(pNum, pNum, CV_8UC1);
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			DD.data[i * pNum + j] = abs((int)(img[i][j] - u[i][j]));
//		}
//	}
//	imwrite("C:\\Users\\team8\\source\\repos\\cplusplus/error.jpg", DD);
//	//imshow("error image", DD);
//	//waitKey();
//
//	return 0;
//}
//
///*关于函数式编程
//* 1. 以前总觉得何必要函数式编程？不清楚它的意义在哪？每次传递一堆参数很烦人，现在通过一遍遍debug明白点了
//* 2. 函数式编程容易debug，尤其是代码一长；函数式编程容易扩展，容易copy调用等等，额外的还有减少代码量吧，但我觉得这不算啥，毕竟可以ctrl+c
//* 3. 通过此次debug我深深的明白了函数式编程的重要性，更明白了开发时测试的重要性。函数式编程一定是和测试在一起的，一个个模块的测试保证了这种方法的开发高效性
//*/
///*关于编程
//* 1. 真的，为了测试正确，你最后还要手算一下小规模运算。那就说白了，编程到最后还是要数学好，编程只不过是个设定机械化的过程
//* 2. 实现这个机械化过程需要学习语言，更需要经验，这是个工程问题
//* 3. 内家功法的学习要重视，本末不能倒置啊
//* 4. 写程序多用代数方法去判断/计算，这样不仅快而且体现数学
// */
//
// /*测试Sidden算法的极度简化方法
// * 1. 网格什么的宽度都设定为1
// * 2. 四个像素，四条射线，两个角度（0和45）
// */
//
// /* 关于编程风格
// * 1. 朴实无华才是王道，计算机只是执行，如果没有特别性能要求结果一样都一样
// * 2. 各种范式尽量给自己找个能避免因误写而造就bug的写法，比如两层for循环，里面的j我经常肌肉记忆顺带写成i，这种写循环方法肯定不好
// *	  这也就和我喜欢简单计算或者判断都加括号一样，我是自己看着没毛病，但是一旦手滑了，程序理解错了就出现bug了，还难以排查
// * 3. 我觉得好的代码注释比代码多或者量级上不会有巨大差距
// */