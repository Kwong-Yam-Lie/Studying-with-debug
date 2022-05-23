#include<opencv2/opencv.hpp>
#include<opencv2/highgui.hpp>
#include<cstdio>
using namespace cv;

int main()
{
	unsigned int width, height;
	width = 4096;
	height = 2048;
	float* hData = new float[width * height];
	unsigned int size = width * height * sizeof(float);

	FILE* fp;
	fp = fopen("C:\\Users\\15867\\Pictures\\imgData.dat", "rb");
	fread(hData, 1, size, fp);
	fclose(fp);
	
	Mat source = Mat(height, width, CV_8UC1);
	for (unsigned int i = 0; i < width * height; ++i) {
		source.data[i] = (int)(hData[i] * 255);
	}
	namedWindow("source", WINDOW_NORMAL);
	imshow("source", source);
	waitKey();


	return 0;
}
