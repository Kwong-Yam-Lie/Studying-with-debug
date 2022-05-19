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
///* ��֮ǰд��Sidden�㷨��Щbug��ͼ��Сʱ�����ȷ��ͼ���ʱ�������ֻ��Ϊ�Լ����ж��߽�����ʱ����һ������Ϊ�ܺõĴ���ɧ����
//*  ��������ߴ������񽻵�ʱ������ش����⣬�����Լ���
//*/
//
///*�������ܼ�ע������
//*	Sidden: ��ⵥ�Ƕ��µ�ϵͳ�Ӿ���ÿ����Sidden��������Ѱ�A���㣬�Է���ʹ�ã�����������ע�⣬����һ��Sidden��A�ĸı�ֵ�ֲ����ô��ڣ�����һ�ε���Siddenʱ���á�
//*	 proj : ��ⵥ�Ƕ���̽�������ϵ�ֵ�������������һ��һά����ʱ����p_temp��
//*	 Diff : ������ͼ���X���ֺ�Y���֣������ڲ������D(),ͳһ��������ֹ淶,��־����Ŀ���Ӧ��ʽʵΪһά����
//*	Diff_T: ������ͼ���x��ת�ò�ֺ�Y��ת�ò�֣������Diff��ͬС�죬�ǲ�־���ζԽ��ߵĶԳ�λ�����ƣ����������ֱ�Ϊ����ǰ��ֵĸ�ֵ��ֻ�Ƿֿ������ÿһ�����β���������
//*	f_solve:���Ŀ�꺯����f = arg{min{|dx|_1+|dy|_1+(mu/2)*pow(|Au-p|_2,2)+(lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)}}
//*  u_update:ʹ�������½����������ݸ��£��������Ϊ��������,����һ�κ����洢��u_temp��
//*	shrink: ����dx,dy,��������ֵ�㷨��ʹ�ô˺�����Dx,Dy������u��
//*/
//
//void Sidden(double** A, const int view, const double da, double point_source_x, double point_source_y, double* bar_x, double* bar_y, const int pNum, const int dNum, const double pSize, const double dSize, double** imgX, double** imgY)
//{
//	for (int i = 0; i < dNum; ++i) { // ����A
//		for (int j = 0; j < pNum * pNum; ++j) {
//			A[i][j] = 0;
//		}
//	}
//	// ��λ���Ƕ��µ����ߴ�
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
//	for (int k = 0; k < dNum; ++k) { // ���Ƕȵ����ߵĿ��ǣ��ɲ���
//		double slope = (rot_bar_y[k] - rot_point_source_y) / (rot_bar_x[k] - rot_point_source_x + 1e-13);
//
//		if (abs(slope) < 1) { // X����
//			double length_through_pixel = pSize * sqrt(1 + pow(slope, 2)); // ����layer֮��Ĺ̶���
//
//			/*������������������д���
//			* 1. ���������߲�����������ж�
//			* 2. �����ǵ�layer��Խ��������ж�
//			*/
//			double* unmain_sequence_value = new double[pNum + 1]{}; // ͼ��ģ�͹������õ��������ģ�����ʹ����������߿�
//			for (int i = 0; i < pNum + 1; ++i) {
//				if (i < pNum) {
//					unmain_sequence_value[i] = slope * ((imgX[0][i] - pSize / 2.0) - rot_point_source_x) + rot_point_source_y;
//				}
//				else {
//					unmain_sequence_value[i] = slope * ((imgX[0][i - 1] + pSize / 2.0) - rot_point_source_x) + rot_point_source_y;
//				}
//			}
//			// �����unmain_sequence_value��Ӧ��Y��ţ�������������ʹ��
//			double* unmain_sequence_number = new double[pNum + 1]{};
//			for (int i = 0; i < pNum + 1; ++i) {
//				unmain_sequence_number[i] = (int)floor(-(unmain_sequence_value[i] - pNum * pSize / 2.0) / pSize);
//			}// ��ʱ��unmain_sequence_number�����и�ֵ����Ҫ����ʹ�á�����Ҳ����ʹ��int����Ϊ�и�ֵ��ʹ�õĻ�������
//
//			for (int i = 0; i < pNum; ++i) { // �жϴ�Խ���صĳ���
//				if (unmain_sequence_number[i] == unmain_sequence_number[i + 1]) {
//					// ��ʱ�����߾���ÿ��layer����һ������ // �ع�ʱ���ǵ��˱߽�����
//					// �״�����ж�������ceil(unmain_sequence_value[i] / pSize) == ceil(unmain_sequence_value[i + 1] / pSize)����д�����ԣ����ٽ����������
//					// 
//					// ��ʱ�������ֿ��ܣ������߶ζ�Ӧ�����������ڣ�������㸳ֵ�������߶ζ�Ӧ���ز��������ڣ���ֵ0����Ϊ��ʼ��Ϊ0���ʿɲ���������
//					if (abs(ceil(unmain_sequence_value[i] / pSize)) <= pNum / 2.0 && abs(floor(unmain_sequence_value[i] / pSize)) <= pNum / 2.0) {
//						// �����ж�����д��< (unmain_sequence_number[i]>=0) && (unmain_sequence_number[i]<=pNum) >��̫�ã���Ϊû���Ǻñ߽����������Ž���ֵ
//						A[k][((int)unmain_sequence_number[i]) * pNum + i] = length_through_pixel; // ����������������intǿ������ת�������򱨴�����֮��
//						// ��һ�����ﱨbug,��Ϊfloor�����ǰ���ǲ�����д�˸�����
//					}// ���ﲻ��(int)�ᱨ��,�����˱���
//				}
//				else { // ��ʱ�����߾���ÿ��layer�����������أ������������ز�һ�������������档
//					// ��ÿ�����ض���������ģ���ֿɷ�Ϊ���������б��Ϊ[0,1)��б��Ϊ(-1,0]��
//					// ������ٽ���������ﲻ���ڣ���ΪX������ǰ���if�ж���
//					// ��ʱ�򻹿����ô�����ʽֱ�Ӽ��㣬һ�������ǣ�floorǰһ��+ceil��һ����/2=�м���Ǹ�Y���
//					// �˺��������������μ��ɼ��㣬���������β�ȡС�����ıȽϣ������ǶԽǵ�������������
//					int  temp_Y{ 0 }; // ����һ��layerʱһ���������⣬һ���������ڣ��������Y������ױȽ��ж�
//					temp_Y = (int)((ceil(unmain_sequence_value[i] / pSize) + floor(unmain_sequence_value[i + 1] / pSize)) / 2.0);
//					// ��������temp_Y�����������Ӧ��һ��Y���꣬�Ǽ�������layer��Y֮���������
//					if (temp_Y == -pNum / 2) { // �߽�����
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
//					else if (temp_Y<(pNum / 2) && temp_Y>(-pNum / 2)) { // �ڲ�����
//						// ���������ع���ģ���Ϊ�ع����ָ����Ͳ���Ҫд��Щ�жϣ���Ϊ��һ������Sidden�㷨��̸����ˣ�����˼
//						A[k][((int)unmain_sequence_number[i]) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						A[k][((int)unmain_sequence_number[i + 1]) * pNum + i] = length_through_pixel - A[k][((int)unmain_sequence_number[i]) * pNum + i];
//						//if (abs(unmain_sequence_value[i] - pSize * pNum / 2.0) < -(temp_Y - pNum / 2.0) * pSize) { // һ��layer���Σ�����������Ǹ���������
//						//	A[k][((int)unmain_sequence_number[i]) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	A[k][((int)unmain_sequence_number[i + 1]) * pNum + i] = length_through_pixel - A[k][((int)unmain_sequence_number[i]) * pNum + i];
//						//	//A[k][(-(temp_Y - pNum / 2) - 1) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	//A[k][-(temp_Y - pNum / 2) * pNum + i] = length_through_pixel - A[k][(-(temp_Y - pNum / 2) - 1) * pNum + i];
//						//}// ���ܵĴ���ʵ����û��Ҫ�Ͳ��ģ�����ʵ������unmain_sequence_number����ȷ
//						//else { // һ��layer���Σ�����������Ǹ���������
//						//	A[k][((int)unmain_sequence_number[i]) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	A[k][((int)unmain_sequence_number[i + 1]) * pNum + i] = length_through_pixel - A[k][((int)unmain_sequence_number[i]) * pNum + i];
//						//	//A[k][-(temp_Y - pNum / 2) * pNum + i] = abs(abs(unmain_sequence_value[i] - pSize * pNum / 2.0) - (-(temp_Y - pNum / 2.0)) * pSize) * length_through_pixel / abs(unmain_sequence_value[i + 1] - unmain_sequence_value[i]);
//						//	//A[k][(-(temp_Y - pNum / 2) - 1) * pNum + i] = length_through_pixel - A[k][-(temp_Y - pNum / 2) * pNum + i];
//						//}
//					}
//					// �ⲿ���ξͲ��ÿ����ˣ��Ѹø�ֵ�ĸ�ֵ�ˣ�ʣ�µľ���û�����ߴ�Խ�����أ����ҳ�ʼ��ʱ����Ϊ��0�������߼�ǡ��
//				}
//			}
//			delete[] unmain_sequence_value; delete[] unmain_sequence_number;
//		}
//
//		else { // Y����
//			// ͬX����ķ������裬���֮�������΢�䶯�Զ�Ӧ
//			double length_through_pixel = pSize * sqrt(1 + pow(slope, -2)); // ע������pow����������-2
//
//			double* unmain_sequence_value = new double[pNum + 1]{}; // ͼ��ģ�͹������õ��������ģ�����ʹ����������߿�
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
//					if (temp_X == -pNum / 2) { // �߽�����
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
//		p_temp[i] = 0; // �����һ��д����������
//	}// �����ε��ǣ���һ���޸ĺ�dNumд����pNum�����Ǳ����ˣ��Ų��˺þòŲ����
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
//			else if (j == (pNum - 1)) Dx_T[i][j] = u[i][j]; // ��һ�������������һ�������Ų��bug
//			else Dx_T[i][j] = u[i][j] - u[i][j + 1];
//
//			if (i == 0) Dy_T[i][j] = -u[i + 1][j]; // ���е��Ժ����������Dy_T�޷���ȡ�ڴ�
//			else if (i == (pNum - 1)) Dy_T[i][j] = u[i][j];
//			else Dy_T[i][j] = u[i][j] - u[i + 1][j];
//			// std::cout << "�ɹ�ʻ��" << i << std::endl; 
//		}
//	}
//}
//
//double f_solve(double** dx, double** dy, double** bx, double** by, double mu, double lamda, double** u, double** Dx, double** Dy, \
//	double* p_temp, double** A, const int view, const double da, double point_source_x, double point_source_y, double* bar_x, double* bar_y, const int pNum, const int dNum, const double pSize, const double dSize, double** imgX, double** imgY, double** img, double** p)
//{
//	// f = arg{min{|dx|_1+|dy|_1+(mu/2)*pow(|Au-p|_2,2)+(lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)}}
//	double f_1{ 0 }, f_2{ 0 }, f_3{ 0 }; // ��f_solve������Ϊ������
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			f_1 += abs(dx[i][j]) + abs(dy[i][j]); // |dx|_1+|dy|_1
//		}
//	}
//
//	proj(p_temp, A, view, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, u); // ��ǰ��u�µ�����
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
//	/*��⹫ʽ
//	*	d^(k)=mu*A'(Au^(k)-p)+lamda*Dx'(Dx(u^(k))+bx^(k)-dx)+lamda*Dy'(Dy(u^(k))+by^(k)-dy)
//	*	alphi^(k) = (d^(k)'*d^(k))/(d^(k)'*H*d^(k)),����HΪHessian����
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
//	// �����һ��д��ʱ����һЩЩ�����⣬�˺�һ���ֶ��е�����
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
//		}// ����ǰ��<Diff(u, Dx, Dy, pNum);>���ò�����Dy��
//	}
//	Diff_T(d_2, Dx_T, Dy_T, pNum); // ����d_2ֵ�仯����Ҫ�������ɶ�Ӧ��Dy_T
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			u_temp[i][j] += lamda * Dy_T[i][j]; // lamda*Dy'(Dy(u^(k))+by^(k)-dy)
//		}
//	}
//	// ����������Ų���alphi
//	Diff(u_temp, Dx, Dy, pNum); // Dx(d^(k)),Dy(d^(k))
//	proj(p_temp, A, view, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, u_temp);//Ad ^ (k)
//	double alphi{ 0 }, alphi_1{ 0 }, alphi_2{ 0 }, alphi_3{ 0 };
//	// alphi ^ (k) = (d ^ (k)'*d^(k))/(d^(k)' * H * d ^ (k)), ����HΪHessian����
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
//			u_temp[i][j] = alphi * u_temp[i][j]; // ����������
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
//	/*����ϵ������ȷ����ҰԲ��������ͼ�����������ͼ��ģ��*/
//	constexpr int sod{ 600 }; constexpr int sdd{ 1000 }; constexpr int dNum{ 512 }; constexpr double dSize{ 2 }; constexpr int pNum{ 256 };
//	// �����һ��д��ʱ��dSizeд��0.2���������д��2�ͱ���
//	constexpr int views{ 60 }; constexpr double pi{ 3.1415926927 };
//	constexpr double da{ 2 * pi / views };
//
//	constexpr double L = (dNum - 1) * dSize / 2;
//	double R = L * sod / sqrt(pow(L, 2) + pow(sdd, 2));
//	double pSize = 2 * R / pNum; // ȷ������ҰԲ�������С��������Ҫ�洢һЩ��Ҫ�����꣺�����ߡ�������
//
//	double point_source_x{ 0 }, point_source_y{ sod }; // ������
//	double* bar_x = new double[dNum];
//	double* bar_y = new double[dNum];
//	for (int i = 0; i < dNum; ++i) {
//		bar_x[i] = -L + i * dSize;
//		bar_y[i] = (double)sod - (double)sdd;
//	}
//	double** imgX = new double* [pNum]; // ������
//	double** imgY = new double* [pNum];
//	for (int i = 0; i < pNum; ++i) {
//		imgX[i] = new  double[pNum];
//		imgY[i] = new double[pNum];
//	}
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) {
//			imgX[i][j] = -R + pSize / 2 + pSize * j; // ע������˵���j
//			imgY[i][j] = R - pSize / 2 - pSize * i;  // �������ź���������ϵ�ǲ�һ�µ�
//		}
//	}
//
//	double** img = new double* [pNum]; // ����ͼ��ģ��
//	for (int i = 0; i < pNum; ++i) {
//		img[i] = new double[pNum];
//	}
//	int x_1 = 0, y_1 = 0, r_1 = 200, density_1 = 50; // ��������Ϊdouble����ʹ��Mat������⣬��Ҫ���Ҷ�Mat����Ϥ
//	int x_2 = 100, y_2 = 100, r_2 = 50, density_2 = 100;
//	for (int i = 0; i < pNum; ++i) {// �Ҷ�������ÿ��ֵ��Ӧһ���������ĵ�
//		for (int j = 0; j < pNum; ++j) {
//			if ((pow(imgX[i][j] - x_1, 2) + pow(imgY[i][j] - y_1, 2)) < pow(r_1, 2)) {
//				img[i][j] = density_1;
//			}
//			else img[i][j] = 0; // �����߼�Ҫ����ã�����Ҫͼ����ʾ��������Ҫͼ����� 
//								// ����൱�������ʼ���ˣ�ͼ��һ
//
//			if ((pow(imgX[i][j] - x_2, 2) + pow(imgY[i][j] - y_2, 2)) < pow(r_2, 2)) {
//				img[i][j] = density_2;
//			}
//			if ((pow(imgX[i][j], 2) + pow(imgY[i][j], 2)) > pow(R - pSize / 2, 2)) // ��ҰԲ�������
//			{
//				img[i][j] = 0;
//			}
//		}
//	}
//
//	/*����OpenCv����ʾͼ��*/
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
//	/*����Sidden�㷨����ϵͳ����
//	* 1. �洢�ռ��С�趨Ϊ���Ƕ��µ�ϵ���������裬Ҫ��Ȼȫ�Ƕ�������Ϊn^4���Բ���
//	* 2. ���캯��Sidden��ͨ�����ݹؼ�����view��ȷ�����򼰾�������
//	*/
//	double** A = new double* [dNum];
//	for (int i = 0; i < dNum; ++i) {
//		A[i] = new double[pNum * pNum]{}; // ϵͳ���󼰳�ʼ��
//	} // �������һ����Sidden�㷨�ĺ���
//
//	/*Split_Bregman�㷨
//	* 1. ��ADMM����ûɶ���𣬲����Ǹı��˳��Ӹ��·�ʽ
//	* 2. �������ʱѡ�������½�����������ppt����Ľ����⣬�Ǹ�ûѧ��������������ڲ���
//	* 3. ��ʼͼ��ֱ��ѡ��ȫ������Ȳ�ʹ��OS-SART����������ͼ��������ͼ�񣬿����ٶ�
//	* 4. ֻʹ��ADMM�����������ԣ�һ���ȴ���򵥵���Ͽ�ʼ
//	* 5. ʹ�ø�������ģ��
//	* 6. ������������������ǲ���Ҫ�ľֲ����Ž⣬�����ϳ�ʼͼ��ѡ��OS-SART�����������ʱ���Ҳ�����Ǿֲ����ţ��ⲻϡ��ʱ��������������Ҫ��
//	*/
//
//	/*Լ��ģ��
//	* f = arg{min{|dx|_1+|dy|_1+(mu/2)*pow(|Au-p|_2,2)+(lamda/2)*pow(|dx-Dx(u)-bx|_2,2)+(lamda/2)*pow(|dy-Dy(u)-by|_2,2)}}
//	* ���� ��f��Ŀ�꺯����Dx�ǲ������,bx��Split_Bregman��������ɫ�����lamda�Ľ�����£�
//	*/
//
//	/*��⹫ʽ
//	* 1. dx,dy������ʼֵ0������ÿ�ε�����u������ֵ��ʽ����dx=shrink(Dx(u)+bx,1/lamda),dy=shrink(Dy(u)+by,1/lamuda)
//	* 2. �����½�������ģ��u^(K+1) = u^(k) - alpha^(k)*d^(k)������dkΪ�ݶȷ���alphiΪ���Ų���
//	*	d^(k)=mu*A'(Au^(k)-p)+lamda*Dx'(Dx(u^(k))+bx^(k)-dx)+lamda*Dy'(Dy(u^(k))+by^(k)-dy)
//	*	alphi^(k) = (d^(k)'*d^(k))/(d^(k)'*H*d^(k)),����HΪHessian����
//	*	H = lamda*Dx'Dx + lamda*Dy'Dy+mu*A'A
//	*-->d^(k)'*H*d^(k) = lamda*pow(|Dx(d^(k))|_2,2)+  lamda*pow(|Dy(d^(k))|_2,2) + mu*pow(|Ad^(k)|_2,2)
//	* 3. Split_Bregman��������bx�ĸ��£����߻�ͷ·����bx^(k+1) = bx^(k) + Dx(u)^(k+1) - dx^(k+1)
//	* 4. �������ĸ��£� mu = 2*mu,��ʼֵ����Ϊ0.05�����ֵ����Ϊpow(2,13)
//	*/
//
//	// Լ��ģ�͵����
//	double** u = new double* [pNum]; double** dx = new double* [pNum]; double** dy = new double* [pNum]; double** bx = new double* [pNum]; double** by = new double* [pNum];
//	for (int i = 0; i < pNum; ++i) {
//		u[i] = new double[pNum] {}; dx[i] = new double[pNum] {}; dy[i] = new double[pNum] {}; bx[i] = new double[pNum] {}; by[i] = new double[pNum] {};
//	}
//
//	constexpr int maxIte_out{ 1000 }; constexpr int maxIte_inn{ 100 };
//	double mu = 0.05; double lamda = pow(2, 5); // ������йصĲ����趨
//
//	double** p = new double* [dNum];
//	for (int i = 0; i < dNum; ++i) {
//		p[i] = new double[views] {};
//	} // �������һ��ʵ�ֵ��Ƕ���ͶӰ�ĺ���
//	double* p_temp = new double[dNum] {}; // ���Ƕ���ͶӰ�������ʱ�洢����
//	for (int i = 0; i < views; ++i) { // ���ﲻ����openmp����ΪA��p_temp����������ֻ��һ��
//		proj(p_temp, A, i, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img);
//		for (int j = 0; j < dNum; ++j) {
//			p[j][i] = p_temp[j];
//		}
//	} // ����õ���ȫ�Ƕ���ʵ��ͶӰ���ݣ������ؽ���������
//
//	/*����OpenCv����ʾͼ��*/
//	// ����ͶӰֵ�ķ�Χ�㲻֪������������Ҫѹ����ʾ
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
//	double** Dx = new double* [pNum]; double** Dy = new double* [pNum];// �������һ��������ʵ�ֲ�ֹ���
//	double** Dx_T = new double* [pNum]; double** Dy_T = new double* [pNum]; // �������һ������ʵ�ֲ������ת�ù���
//	for (int i = 0; i < pNum; ++i) {
//		Dx[i] = new double[pNum] {};
//		Dy[i] = new double[pNum] {};
//		Dx_T[i] = new double[pNum] {};
//		Dy_T[i] = new double[pNum] {}; // ��һ���������Dy_Tд����Dx_T������Ų��˰��죬��ĸ���
//	}
//
//	double f_current{ 0 }, f_previous{ 0 }; // �������һ������������Ŀ�꺯����һ����ʾ������Ľ⣬һ����ʾ����ǰ�Ľ�
//
//	double** u_temp = new double* [pNum]; // ����һ���м������洢��ʱ���ͼ��������ݶȸ��µ�����
//	for (int i = 0; i < pNum; ++i) {
//		u_temp[i] = new double[pNum] {};
//	} // �������һ���������������ݶȸ�����
//	double error{ 0 }, error_1{ 0 }, error_2{ 0 }; // �����жϳ�����ֹ����
//
//	// ������ʼ
//	ULONGLONG t1, t2;
//	t1 = GetTickCount64();
//	for (int ite_out = 0; ite_out < maxIte_out; ++ite_out) { // �趨�ܹ��������ٴ�
//		// �����������,ͨ��������������е�������ʵ�����Ч��
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
//		for (int k = 0; k < views; ++k) { // ��Ƕȵ�����ÿ���Ƕ�ʹ��Splite_Bregman����100��
//			// ����ѡ�Ƕ�Ҫ������ٶȻ������
//			k = ran_views[k];
//			for (int ite = 0; ite < maxIte_inn; ++ite) {
//				f_current = f_solve(dx, dy, bx, by, mu, lamda, u, Dx, Dy, \
//					p_temp, A, views, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//
//				u_update(u_temp, Dx_T, Dy_T, dx, dy, bx, by, mu, lamda, u, Dx, Dy, \
//					p_temp, A, k, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						u_temp[i][j] = u[i][j] - u_temp[i][j]; // ��ʱu_temp�������ݶȸ��º��ͼ��
//						if (u_temp[i][j] < 0) u_temp[i][j] = 0; // �Ǹ�Լ��
//					}
//				}
//				f_previous = f_current;
//				f_current = f_solve(dx, dy, bx, by, mu, lamda, u_temp, Dx, Dy, \
//					p_temp, A, views, da, point_source_x, point_source_y, bar_x, bar_y, pNum, dNum, pSize, dSize, imgX, imgY, img, p);
//
//				if (f_current >= f_previous) { // ��һ���˹��ж���ʹÿ��ѭ��������Ч
//					for (int t = 0; t < 5; t++) { // �������㣬�������˹��Ż����
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
//						error_1 += pow((u_temp[i][j] - u[i][j]), 2) / 255; // ��ֹ�ۼӳ���double����׼����
//						error_2 += pow(u[i][j], 2) / 255;
//					}
//				}error = error_1 / error_2;
//
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						u[i][j] = u_temp[i][j];
//					}
//				}
//				/***************����u�������*,����������dx,dy��bx,by***********/
//					// dx=shrink(Dx(u)+bx,1/lamda),dy=shrink(Dy(u)+by,1/lamda)
//					// bx^(k+1) = bx^(k) + Dx(u)^(k+1) - dx^(k+1)
//					// ��Ҫд������dx,dy������ֵ����
//				shrink(dx, dy, u, Dx, Dy, bx, by, pNum, lamda);
//				for (int i = 0; i < pNum; ++i) {
//					for (int j = 0; j < pNum; ++j) {
//						bx[i][j] = bx[i][j] + Dx[i][j] - dx[i][j];
//						by[i][j] = by[i][j] + Dy[i][j] - dy[i][j]; // ����ֱ������shrink�����ı��Dx,Dy��ֵ
//					}
//				}
//				/***************����dx,dy,bx,by�������*,���������·�����mu***********/
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
//	cout << "Split_Bregman��������ʱ�� " << (t2 - t1) / (1000 * 60.0) << "min" << endl;
//
//	ofstream fout("Split_Bregman.txt", ios_base::out | ios_base::app);
//	fout << "Split_Bregman��������ʱ�� " << (t2 - t1) / (1000 * 60.0) << "min" << endl;
//	fout.close();
//
//	// ��ʾͼ��ʱ������double���ͣ���Ҫע��
//
//		/*����OpenCv����ʾͼ��*/
//	Mat EE = Mat(pNum, pNum, CV_8UC1);
//	// ����double����img����Ϊcv_64FC1���������ת�����󣬵��Ҳ��˽�Mat,���Բ��������Ϊʲô
//	for (int i = 0; i < pNum; ++i) {
//		for (int j = 0; j < pNum; ++j) // ֮ǰ�����jд��i�ˣ�������ճ�����Ҿ�����д�ڶ��飡����
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
///*���ں���ʽ���
//* 1. ��ǰ�ܾ��úα�Ҫ����ʽ��̣�����������������ģ�ÿ�δ���һ�Ѳ����ܷ��ˣ�����ͨ��һ���debug���׵���
//* 2. ����ʽ�������debug�������Ǵ���һ��������ʽ���������չ������copy���õȵȣ�����Ļ��м��ٴ������ɣ����Ҿ����ⲻ��ɶ���Ͼ�����ctrl+c
//* 3. ͨ���˴�debug������������˺���ʽ��̵���Ҫ�ԣ��������˿���ʱ���Ե���Ҫ�ԡ�����ʽ���һ���ǺͲ�����һ��ģ�һ����ģ��Ĳ��Ա�֤�����ַ����Ŀ�����Ч��
//*/
///*���ڱ��
//* 1. ��ģ�Ϊ�˲�����ȷ�������Ҫ����һ��С��ģ���㡣�Ǿ�˵���ˣ���̵������Ҫ��ѧ�ã����ֻ�����Ǹ��趨��е���Ĺ���
//* 2. ʵ�������е��������Ҫѧϰ���ԣ�����Ҫ���飬���Ǹ���������
//* 3. �ڼҹ�����ѧϰҪ���ӣ���ĩ���ܵ��ð�
//* 4. д������ô�������ȥ�ж�/���㣬�������������������ѧ
// */
//
// /*����Sidden�㷨�ļ��ȼ򻯷���
// * 1. ����ʲô�Ŀ�ȶ��趨Ϊ1
// * 2. �ĸ����أ��������ߣ������Ƕȣ�0��45��
// */
//
// /* ���ڱ�̷��
// * 1. ��ʵ�޻����������������ֻ��ִ�У����û���ر�����Ҫ����һ����һ��
// * 2. ���ַ�ʽ�������Լ��Ҹ��ܱ�������д�����bug��д������������forѭ���������j�Ҿ����������˳��д��i������дѭ�������϶�����
// *	  ��Ҳ�ͺ���ϲ���򵥼�������ж϶�������һ���������Լ�����ûë��������һ���ֻ��ˣ����������˾ͳ���bug�ˣ��������Ų�
// * 3. �Ҿ��úõĴ���ע�ͱȴ������������ϲ����о޴���
// */