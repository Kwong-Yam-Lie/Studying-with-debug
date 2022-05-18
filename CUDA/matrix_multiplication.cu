//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include<stdlib.h>
//#include<windows.h>
//#include<cmath>
//#include<iostream>
//#include<omp.h>
//using namespace std;
//#define num_thread 29 // 单block线程数 // 注意一个奇怪点(num_thread,num)->(30,1000)速度比(32,1000)快10倍数
//#define num 3000 // 方针大小 
////经过测试，此服务器选取num_thread=29性能最佳，num=1000时加速效果 1300倍；num=3000时加速效果 7500倍
////真的，设计是一回事，调参是一回事。设计的再好也要调调参数，反过来，设计差不多就去调参数，效果也可以调的不错，辩证看待吧。
//
///*shared memory 的使用主要是为了减少对全局内存的重复访问次数，以此为设计标准*/
///*为了加快计算速度牺牲一部分精度，使用float类型。int类型优化不好，double类型字节数加倍。虽然float类型会有误差损失，但是以平均误差来考量还可以接受*/
///*缺点：kernel函数太多if语句*/
///*"极度简化"思想用在程序调bug上是真真的好使！！！*/
//
///*历经千难万险，最终查出了bug:
//* 1. 把blockDim当成gridDim使了；
//* 2. 把人工划分的执行块的维度大小等价于实际程序的block大小了，比如这个程序是把dev_C划分为多个执行块，每个块大小(num_thread,num_thread)
//*	 但是，实际上一个线程执行一个维度，我定义的block是(num_thread,num_thread)，所以使用blockDim.y是1而不是num_thread，造成混淆错误
//*/
//
//__global__ void matrix_mul_kernel(float* dev_C, float* dev_A, float* dev_B, int size, int weight, int height)
//{
//	/* allocation shared memory */ //(下面是最初的分析版本，实际上很多理解不到位)
//	
//	/*设备情况
//	* 使用GPU device 0: Tesla V100S-PCIE-32GB
//	* SM的数量：80
//	* 每个线程块的共享内存大小：48 KB
//	* 每个线程块的最大线程数：1024
//	* 每个SM的最大线程数：2048
//	* 每个SM的最大线程束数：64
//	*/
//
//	/*分析
//	* 1. 单个函数默认使用静态内存为16k，超出则警告，后续有问题不好debug，需要衡量是否必要超出。
//	*    片上内存48k，要被block充分利用。若每个 SM 3个block，则按1的话单次最大3*80=240个block。
//	*    但是一个block执行完另一个block立即跟上去，所以block数量不用考虑了，着重于从内存使用上考虑。
//	*	 结果就是：每个SM 3个block，每个block 2048/3=682 -> 32*21=672个thread，最大16k的共享内存
//	*
//	* 2. 每个block的共享内存都有三个矩阵块，(pow(x,2)+2*x)*4.0/1024=16,解得整数解x=63
//	*	 那么每个block应该有至多63个thread，那么可以看出1和2是相互制约的。
//		 应该怎样选取合适的num_threads?应该选取什么维度的block？
//	*
//	* 3. 措施一：单block超出16k警告，直接拉满，这时解得整数解x=109 -> x=32*3=96，单个SM单次执行一个block。
//	*	 措施二：把寄存器的16k拿过来，但是这样不太好（可能拖慢速度），我也把控不住。
//	*	 措施三：选取1维的block，二维的话num_threads要开方，尽管二维的num_threads好设计程序，但这两种对全局变量访问的重复次数差多了
//	*
//	* 4. 分析结果：Grid上二维block;blockDim.y=96,blockDim.x=1;每个线程单次执行1个像素操作，迭代weight次。
//	*/
//	constexpr int N{ num_thread };
//	__shared__ float temp_A[N * 1]; // 96，之后再试验64
//	__shared__ float temp_B[1 * N];
//	__shared__ float temp_C[N][N]; // [size_y][size_x] 二维是这样定义的？
//
//	for (int i = 0; i < N; ++i) { // 先把temp_C初始化，后面好用“+=”运算符
//		temp_C[threadIdx.x][i] = 0;
//	}
//	__syncthreads();
//
//	/*compute in each thread*/
//	if (blockIdx.x < (gridDim.x - 1)) {// 列边界条件，不含最右列
//		for (int count = 0; count < weight; ++count) {
//			temp_B[threadIdx.x] = dev_B[count * weight + (blockDim.x * blockIdx.x + threadIdx.x)];
//			if ((blockIdx.y * num_thread + threadIdx.x) < height) {
//				temp_A[threadIdx.x] = dev_A[(num_thread * blockIdx.y + threadIdx.x) * weight + count];
//			}// 要不这里计算赋值错误，要不最后数据传回去错误
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
//			}// 这里即使给它if判断条件放宽且赋值多加个1，在(num_thread,num)->(3,6)时最后两行也不会输出。
//		}	 // 线程还会挂起不执行？但是测试了这么多、对程序更改了这么多，但结果就是这个！没理由啊，真是怪了
//	}
//
//	if (blockIdx.x == (gridDim.x - 1)) {// 最右列
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
//	dim3 threadsPerBlock(num_thread); // 行分成多少block，每个block的thread代表一行
//	// threadsPerblock.x越大重复访问全局内存次数越少，因此这里设置大点好，但不是越大越好，分析见kernel
//	dim3 numBlocks((int)ceil((1.0 * height) / threadsPerBlock.x), (int)ceil((1.0 * height) / threadsPerBlock.x));
//	cout << "Block: " << (int)ceil((1.0 * height) / threadsPerBlock.x) << "," << (int)ceil((1.0 * height) / threadsPerBlock.x) << endl;
//	// 根据目标矩阵dev_C来划分block，每个block只需要96个线程即可求解
//
//	//使用event计算时间
//	float time_elapsed{ 0 };
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start); //创建Event
//	cudaEventCreate(&stop);
//
//	cudaEventRecord(start, 0); //记录当前时间
//	/*kernek : calling function*/
//	matrix_mul_kernel << <numBlocks, threadsPerBlock >> > (dev_C, dev_A, dev_B, size, weight, height);
//	cudaDeviceSynchronize();
//	cudaError_t err = cudaGetLastError();
//	if (err != cudaSuccess) {
//		printf("CUDA ERROR: %s\n", cudaGetErrorString(err));
//	}
//	cudaEventRecord(stop, 0); //记录当前时间
//
//	cudaEventSynchronize(start); //Waits for an event to complete.
//	cudaEventSynchronize(stop);  //Waits for an event to complete.Record之前的任务
//	cudaEventElapsedTime(&time_elapsed, start, stop); //计算时间差
//	cout << "CUDA计算时间（不含数据传输）： " << time_elapsed << "ms" << endl;
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
//		//A[i] = (rand() % 4) / ((float)1); // 用来测试绝对误差标准
//		//B[i] = (rand() % 4) / ((float)1);
//		A[i] = (rand() % 256) / ((float)255); // 用来测试通用误差标准
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
//	/*串行矩阵乘法计算*/
//	ULONGLONG t1, t2;
//	t1 = GetTickCount64();
//	float* CC = new float[size] {};
//	// nvcc 编译 parallel for 反而慢了
//	// Runtime工程：2000*2000要75秒，1000*1000要2秒，数据大了反而慢了
//	// CPP工程：2000*2000要5秒，1000*1000要0.15秒
////#pragma omp parallel for  
//	for (int i = 0; i < height; ++i) { // 立方级别倍增
//		for (int j = 0; j < weight; ++j) {
//			for (int col = 0; col < weight; ++col) {
//				CC[i * weight + j] += A[i * weight + col] * B[col * weight + j];
//			}
//		}
//	}// 3000*3000的串行需要时间 318s
//	t2 = GetTickCount64();
//	cout << "串行矩阵乘法计算耗时： " << t2 - t1 << "ms" << endl;
//
//	/*利用CUDA进行矩阵乘法计算*/
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
///*调试此程序感悟：自己写的很多程序逻辑上可以简化从而减少代码量，这样不仅方便阅读也方便debug*/
///*逻辑性内容的调试程序手段： 有规律改变赋值、禁用部分逻辑等，观察预期和想的是否一样*/
///*我发现我的bug基本都不是逻辑bug，几乎都是打错字或者不小心没留意肌肉记忆顺写弄成的错误！关键是这些刻在自己身体上，还不好改！！*/
