//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include<stdlib.h>
//#include<windows.h>
//#include<cmath>
//#include<iostream>
//#include<omp.h>
//using namespace std;
//#define num_thread 29 // ��block�߳��� // ע��һ����ֵ�(num_thread,num)->(30,1000)�ٶȱ�(32,1000)��10����
//#define num 3000 // �����С 
////�������ԣ��˷�����ѡȡnum_thread=29������ѣ�num=1000ʱ����Ч�� 1300����num=3000ʱ����Ч�� 7500��
////��ģ������һ���£�������һ���¡���Ƶ��ٺ�ҲҪ��������������������Ʋ���ȥ��������Ч��Ҳ���Ե��Ĳ�����֤�����ɡ�
//
///*shared memory ��ʹ����Ҫ��Ϊ�˼��ٶ�ȫ���ڴ���ظ����ʴ������Դ�Ϊ��Ʊ�׼*/
///*Ϊ�˼ӿ�����ٶ�����һ���־��ȣ�ʹ��float���͡�int�����Ż����ã�double�����ֽ����ӱ�����Ȼfloat���ͻ��������ʧ��������ƽ����������������Խ���*/
///*ȱ�㣺kernel����̫��if���*/
///*"���ȼ�"˼�����ڳ����bug��������ĺ�ʹ������*/
//
///*����ǧ�����գ����ղ����bug:
//* 1. ��blockDim����gridDimʹ�ˣ�
//* 2. ���˹����ֵ�ִ�п��ά�ȴ�С�ȼ���ʵ�ʳ����block��С�ˣ�������������ǰ�dev_C����Ϊ���ִ�п飬ÿ�����С(num_thread,num_thread)
//*	 ���ǣ�ʵ����һ���߳�ִ��һ��ά�ȣ��Ҷ����block��(num_thread,num_thread)������ʹ��blockDim.y��1������num_thread����ɻ�������
//*/
//
//__global__ void matrix_mul_kernel(float* dev_C, float* dev_A, float* dev_B, int size, int weight, int height)
//{
//	/* allocation shared memory */ //(����������ķ����汾��ʵ���Ϻܶ���ⲻ��λ)
//	
//	/*�豸���
//	* ʹ��GPU device 0: Tesla V100S-PCIE-32GB
//	* SM��������80
//	* ÿ���߳̿�Ĺ����ڴ��С��48 KB
//	* ÿ���߳̿������߳�����1024
//	* ÿ��SM������߳�����2048
//	* ÿ��SM������߳�������64
//	*/
//
//	/*����
//	* 1. ��������Ĭ��ʹ�þ�̬�ڴ�Ϊ16k�������򾯸棬���������ⲻ��debug����Ҫ�����Ƿ��Ҫ������
//	*    Ƭ���ڴ�48k��Ҫ��block������á���ÿ�� SM 3��block����1�Ļ��������3*80=240��block��
//	*    ����һ��blockִ������һ��block��������ȥ������block�������ÿ����ˣ������ڴ��ڴ�ʹ���Ͽ��ǡ�
//	*	 ������ǣ�ÿ��SM 3��block��ÿ��block 2048/3=682 -> 32*21=672��thread�����16k�Ĺ����ڴ�
//	*
//	* 2. ÿ��block�Ĺ����ڴ涼����������飬(pow(x,2)+2*x)*4.0/1024=16,���������x=63
//	*	 ��ôÿ��blockӦ��������63��thread����ô���Կ���1��2���໥��Լ�ġ�
//		 Ӧ������ѡȡ���ʵ�num_threads?Ӧ��ѡȡʲôά�ȵ�block��
//	*
//	* 3. ��ʩһ����block����16k���棬ֱ����������ʱ���������x=109 -> x=32*3=96������SM����ִ��һ��block��
//	*	 ��ʩ�����ѼĴ�����16k�ù���������������̫�ã����������ٶȣ�����Ҳ�ѿز�ס��
//	*	 ��ʩ����ѡȡ1ά��block����ά�Ļ�num_threadsҪ���������ܶ�ά��num_threads����Ƴ��򣬵������ֶ�ȫ�ֱ������ʵ��ظ����������
//	*
//	* 4. ���������Grid�϶�άblock;blockDim.y=96,blockDim.x=1;ÿ���̵߳���ִ��1�����ز���������weight�Ρ�
//	*/
//	constexpr int N{ num_thread };
//	__shared__ float temp_A[N * 1]; // 96��֮��������64
//	__shared__ float temp_B[1 * N];
//	__shared__ float temp_C[N][N]; // [size_y][size_x] ��ά����������ģ�
//
//	for (int i = 0; i < N; ++i) { // �Ȱ�temp_C��ʼ����������á�+=�������
//		temp_C[threadIdx.x][i] = 0;
//	}
//	__syncthreads();
//
//	/*compute in each thread*/
//	if (blockIdx.x < (gridDim.x - 1)) {// �б߽�����������������
//		for (int count = 0; count < weight; ++count) {
//			temp_B[threadIdx.x] = dev_B[count * weight + (blockDim.x * blockIdx.x + threadIdx.x)];
//			if ((blockIdx.y * num_thread + threadIdx.x) < height) {
//				temp_A[threadIdx.x] = dev_A[(num_thread * blockIdx.y + threadIdx.x) * weight + count];
//			}// Ҫ��������㸳ֵ����Ҫ��������ݴ���ȥ����
//			else temp_A[threadIdx.x] = 10;
//			__syncthreads();
//
//			for (int col = 0; col < N; ++col) {
//				temp_C[threadIdx.x][col] += temp_A[threadIdx.x] * temp_B[col];
//			}
//			__syncthreads();
//		}
//		__syncthreads();
//		if ((blockIdx.y * num_thread + threadIdx.x) < height) {
//			for (int col = 0; col < N; ++col) {
//				dev_C[(blockIdx.y * num_thread + threadIdx.x) * weight + (blockIdx.x * blockDim.x + col)] = temp_C[threadIdx.x][col];
//			}// ���Ｔʹ����if�ж������ſ��Ҹ�ֵ��Ӹ�1����(num_thread,num)->(3,6)ʱ�������Ҳ���������
//		}	 // �̻߳������ִ�У����ǲ�������ô�ࡢ�Գ����������ô�࣬��������������û���ɰ������ǹ���
//	}
//
//	if (blockIdx.x == (gridDim.x - 1)) {// ������
//		for (int count = 0; count < weight; ++count) {
//			temp_A[threadIdx.x] = dev_A[(num_thread * blockIdx.y + threadIdx.x) * weight + count];
//			if ((blockIdx.x * blockDim.x + threadIdx.x) < weight) {
//				temp_B[threadIdx.x] = dev_B[count * weight + (blockDim.x * blockIdx.x + threadIdx.x)];
//			}
//			else temp_B[threadIdx.x] = 0;
//			__syncthreads();
//
//			for (int col = 0; col < N; ++col) {
//				temp_C[threadIdx.x][col] += temp_A[threadIdx.x] * temp_B[col];
//			}
//			__syncthreads();
//		}
//		__syncthreads();
//		for (int col = 0; col < (weight - blockIdx.x * blockDim.x); ++col) {
//			dev_C[(blockIdx.y * num_thread + threadIdx.x) * weight + (blockIdx.x * blockDim.x + col)] = temp_C[threadIdx.x][col];
//		}
//
//	}
//}
//
//void matrix_mul_withCuda(float* C, float* A, float* B, int size, int weight, int height)
//{
//	/*Device configuration*/
//	cudaSetDevice(0);
//
//	float* dev_C; float* dev_A; float* dev_B;
//	cudaMalloc((void**)&dev_C, size * sizeof(float));
//	cudaMalloc((void**)&dev_A, size * sizeof(float));
//	cudaMalloc((void**)&dev_B, size * sizeof(float));
//
//	cudaMemset(dev_C, 0, size * sizeof(float));
//	cudaMemcpy(dev_A, A, size * sizeof(float), cudaMemcpyHostToDevice);
//	cudaMemcpy(dev_B, B, size * sizeof(float), cudaMemcpyHostToDevice);
//
//	/*kernel : configuration*/
//	dim3 threadsPerBlock(num_thread); // �зֳɶ���block��ÿ��block��thread����һ��
//	// threadsPerblock.xԽ���ظ�����ȫ���ڴ����Խ�٣�����������ô��ã�������Խ��Խ�ã�������kernel
//	dim3 numBlocks((int)ceil((1.0 * height) / threadsPerBlock.x), (int)ceil((1.0 * height) / threadsPerBlock.x));
//	cout << "Block: " << (int)ceil((1.0 * height) / threadsPerBlock.x) << "," << (int)ceil((1.0 * height) / threadsPerBlock.x) << endl;
//	// ����Ŀ�����dev_C������block��ÿ��blockֻ��Ҫ96���̼߳������
//
//	//ʹ��event����ʱ��
//	float time_elapsed{ 0 };
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start); //����Event
//	cudaEventCreate(&stop);
//
//	cudaEventRecord(start, 0); //��¼��ǰʱ��
//	/*kernek : calling function*/
//	matrix_mul_kernel << <numBlocks, threadsPerBlock >> > (dev_C, dev_A, dev_B, size, weight, height);
//	cudaDeviceSynchronize();
//	cudaError_t err = cudaGetLastError();
//	if (err != cudaSuccess) {
//		printf("CUDA ERROR: %s\n", cudaGetErrorString(err));
//	}
//	cudaEventRecord(stop, 0); //��¼��ǰʱ��
//
//	cudaEventSynchronize(start); //Waits for an event to complete.
//	cudaEventSynchronize(stop);  //Waits for an event to complete.Record֮ǰ������
//	cudaEventElapsedTime(&time_elapsed, start, stop); //����ʱ���
//	cout << "CUDA����ʱ�䣨�������ݴ��䣩�� " << time_elapsed << "ms" << endl;
//	cudaEventDestroy(start);
//	cudaEventDestroy(stop);
//
//	/*Device To Host*/
//	cudaMemcpy(C, dev_C, size * sizeof(float), cudaMemcpyDeviceToHost);
//
//	/*Free Memory*/
//	cudaFree(dev_A);	cudaFree(dev_B);	cudaFree(dev_C);
//
//}
//
//int main()
//{
//	constexpr int weight{ num }; constexpr int height{ weight };
//	constexpr int size{ weight * height };
//	float* A = new float[size] {};
//	float* B = new float[size] {};
//	float* C = new float[size] {};
//	srand((unsigned)time(NULL));
//	for (int i = 0; i < size; ++i) {
//		//A[i] = (rand() % 4) / ((float)1); // �������Ծ�������׼
//		//B[i] = (rand() % 4) / ((float)1);
//		A[i] = (rand() % 256) / ((float)255); // ��������ͨ������׼
//		B[i] = (rand() % 256) / ((float)255);
//	}
//	//cout << "*****A********" << endl;
//	//for (int i = 0; i < size; ++i) {
//	//	cout << A[i] << "\t";
//	//	if ((i + 2) % weight == 1) cout << endl;
//	//}
//	//cout << "*****B********" << endl;
//	//for (int i = 0; i < size; ++i) {
//	//	cout << B[i] << "\t";
//	//	if ((i + 2) % weight == 1) cout << endl;
//	//}
//	/*���о���˷�����*/
//	ULONGLONG t1, t2;
//	t1 = GetTickCount64();
//	float* CC = new float[size] {};
//	// nvcc ���� parallel for ��������
//	// Runtime���̣�2000*2000Ҫ75�룬1000*1000Ҫ2�룬���ݴ��˷�������
//	// CPP���̣�2000*2000Ҫ5�룬1000*1000Ҫ0.15��
////#pragma omp parallel for  
//	for (int i = 0; i < height; ++i) { // ����������
//		for (int j = 0; j < weight; ++j) {
//			for (int col = 0; col < weight; ++col) {
//				CC[i * weight + j] += A[i * weight + col] * B[col * weight + j];
//			}
//		}
//	}// 3000*3000�Ĵ�����Ҫʱ�� 318s
//	t2 = GetTickCount64();
//	cout << "���о���˷������ʱ�� " << t2 - t1 << "ms" << endl;
//
//	/*����CUDA���о���˷�����*/
//	matrix_mul_withCuda(C, A, B, size, weight, height);
//
//	double error{ 0 };
//	for (int i = 0; i < height; ++i) {
//		for (int j = 0; j < weight; ++j) {
//			error += abs(CC[i * weight + j] - C[i * weight + j]);
//		}
//	}
//	cout << "sum of absolute error: " << error << endl;
//	cout << "mean error: " << error / size << endl;
//	//cout << "*****CC********" << endl;
//	//for (int i = 0; i < size; ++i) {
//	//	cout << CC[i] << "\t";
//	//	if ((i + 2) % weight == 1) cout << endl;
//	//}
//	//cout << "*****C********" << endl;
//	//for (int i = 0; i < size; ++i) {
//	//	cout << C[i] << "\t";
//	//	if ((i + 2) % weight == 1) cout << endl;
//	//}
//
//	return 0;
//}
///*���Դ˳�������Լ�д�ĺܶ�����߼��Ͽ��Լ򻯴Ӷ����ٴ��������������������Ķ�Ҳ����debug*/
///*�߼������ݵĵ��Գ����ֶΣ� �й��ɸı丳ֵ�����ò����߼��ȣ��۲�Ԥ�ں�����Ƿ�һ��*/
///*�ҷ����ҵ�bug�����������߼�bug���������Ǵ���ֻ��߲�С��û���⼡�����˳дŪ�ɵĴ��󣡹ؼ�����Щ�����Լ������ϣ������øģ���*/