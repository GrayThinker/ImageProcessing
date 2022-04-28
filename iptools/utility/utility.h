#ifndef UTILITY_H
#define UTILITY_H

#include "../image/image.h"
#include <opencv2/opencv.hpp>
#include <sstream>
#include <math.h>

class utility
{
	public:
		utility();
		virtual ~utility();
		static std::string intToString(int number);
		static int checkValue(int value);
		static int hsiCheckValue(int value, int low, int high);
		static void addGrey(image &src, image &tgt, int value);
		static void binarize(image &src, image &tgt, int threshold);
		static void scale(image &src, image &tgt, float ratio);
		static void dualthres(image &src, image &tgt, int T, int V1, int V2);
		static void reg2dsmooth(image &src, image &tgt, int ws);
		static void sep2dsmooth(image &src, image &tgt, int ws);
		static void colorbirght(image &src, image &tgt, int dr, int dg, int db);
		static void colorvisual(image &src, image &tgt, int t, int v);
		static void histostretch(image &src, image &tgt, int A, int B);
		static void althistostretch(image &src, image &tgt, int A, int B);
		static void histothres(image &src, image &tgt, int p1, char p2, int A, int B);
		static void colorstretch(image &src, image &tgt, char p1, int low, int high);
		static vector<int> rgbtohsi(int R, int G, int B);
		static vector<int> hsitorgb(int H, int S, int I);
		static void hsistretch(image &src, image &tgt, char HSI, int low, int high);
		static void fullhsistretch(image &src, image &tgt, char HSI, int low, int high);
		static void histoplot(image &src, image &tgt);  // to plot histogram - not working
		static void sobel3(image &src, image &tgt);
		static void sobel5(image &src, image &tgt);
		static void binaryedge(image &src, image &tgt, int TH, int TL, bool show_angles);
		static void gausobel(image &src, image &tgt, double sigma);
		static void sobelcv(image &src, image &tgt, int win_size);
		static void cannycv(image &src, image &tgt);
		static void hsvcanny(image &src, image &tgt);

		// new
		static void histostretchcv(image &src, image &tgt);
		static void histoequalcv(image &src, image &tgt);
		static void gausscv(image &src, image &tgt);
		static void houghtrans(image &src, image &tgt);
		static void hesobel(image &src, image &tgt);
		static void hecanny(image &src, image &tgt);
		static void qrdecode(image &src, image &tgt);
};

#endif

