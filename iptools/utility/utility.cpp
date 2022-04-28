#include "utility.h"
#include "../image/image.h"
#include <math.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <chrono>
#include <ctime>

#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/core.hpp"
#include <opencv2/objdetect.hpp>

#define MAXRGB 255
#define MINRGB 0
# define M_PI  3.14159265358979323846

using namespace std;
using namespace cv;

std::string utility::intToString(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

int utility::checkValue(int value)
{
	if (value > MAXRGB)
		return MAXRGB;
	if (value < MINRGB)
		return MINRGB;
	return value;
}

int utility::hsiCheckValue(int value, int low, int high){
	if (value > high){
		return high;
	}
	if (value < low){
		return low;
	}
	return value;
}

/*-----------------------------------------------------------------------**/
void utility::addGrey(image &src, image &tgt, int value)
{
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++)
		for (int j=0; j<src.getNumberOfColumns(); j++)
		{
			tgt.setPixel(i,j,checkValue(src.getPixel(i,j)+value)); 
		}
}

/*-----------------------------------------------------------------------**/
void utility::binarize(image &src, image &tgt, int threshold)
{
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++)
	{
		for (int j=0; j<src.getNumberOfColumns(); j++)
		{
			if (src.getPixel(i,j) < threshold)
				tgt.setPixel(i,j,MINRGB);
			else
				tgt.setPixel(i,j,MAXRGB);
		}
	}
}

/*-----------------------------------------------------------------------**/
void utility::scale(image &src, image &tgt, float ratio)
{
	int rows = (int)((float)src.getNumberOfRows() * ratio);
	int cols  = (int)((float)src.getNumberOfColumns() * ratio);
	tgt.resize(rows, cols);
	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<cols; j++)
		{	
			/* Map the pixel of new image back to original image */
			int i2 = (int)floor((float)i/ratio);
			int j2 = (int)floor((float)j/ratio);
			if (ratio == 2) {
				/* Directly copy the value */
				tgt.setPixel(i,j,checkValue(src.getPixel(i2,j2)));
			}

			if (ratio == 0.5) {
				/* Average the values of four pixels */
				int value = src.getPixel(i2,j2) + src.getPixel(i2,j2+1) + src.getPixel(i2+1,j2) + src.getPixel(i2+1,j2+1);
				tgt.setPixel(i,j,checkValue(value/4));
			}
		}
	}
}


/*-----------------------------------------------------------------------**/
void utility::dualthres(image &src, image &tgt, int T, int V1, int V2)
{
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++)
	{
		for (int j=0; j<src.getNumberOfColumns(); j++)
		{
			if (src.getPixel(i,j) > T)
				tgt.setPixel(i,j,checkValue(src.getPixel(i,j)+V1));
			else if (src.getPixel(i,j) < T)
				tgt.setPixel(i,j,checkValue(src.getPixel(i,j)-V2));
		}
	}	
}


/*-----------------------------------------------------------------------**/
void utility::reg2dsmooth(image &src, image &tgt, int ws){
	cout << "function called\n";
	if (ws % 2 != 1){
		ws += 1;
	}
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	// copy entire image for edge values
	for (int i=0; i<src.getNumberOfRows(); i++)
	{
		for (int j=0; j<src.getNumberOfColumns(); j++){
			tgt.setPixel(i,j, src.getPixel(i,j));
		}
	}
	cout << "copying done\n";

	for(int x = floor(ws/2); x < src.getNumberOfRows() - floor(ws/2); ++x){
		for(int y = floor(ws/2); y < src.getNumberOfColumns() - floor(ws/2); ++y){
			int sum = 0;
			int n = 0;
			for(int i = x - floor(ws/2); i < x + floor(ws/2)+1; ++i){
				for(int j = y - floor(ws/2); j < y + floor(ws/2)+1; ++j){
					sum += src.getPixel(i, j);
					n += 1;
				}
			}
			tgt.setPixel(x, y, round(sum/n));
		}
	}
	cout << "smoothing done\n";
}


void utility::sep2dsmooth(image &src, image &tgt, int ws){
	if (ws % 2 != 1){
		ws += 1;
	}
	image temp;
	temp.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	// copy entire image for edge values
	for (int i=0; i<src.getNumberOfRows(); i++)
	{
		for (int j=0; j<src.getNumberOfColumns(); j++){
			temp.setPixel(i,j, src.getPixel(i,j));
			tgt.setPixel(i,j, src.getPixel(i,j));
		}
	}
	cout << "copying done\n";

	// first pass
	for(int j = 0; j < src.getNumberOfColumns(); ++j){
		for(int x = floor(ws/2); x < src.getNumberOfRows() - floor(ws/2); ++x){
			int sum = 0;
			int n = 0;
			for(int i = x - floor(ws/2); i < x + floor(ws/2)+1; ++i){
				sum += src.getPixel(i, j);
				n += 1;
			}
			temp.setPixel(x, j, round(sum/n));
		}
	}

	// second pass
	for(int i = 0; i < src.getNumberOfRows(); ++i){
		for(int y = floor(ws/2); y < src.getNumberOfColumns() - floor(ws/2); ++y){
			int sum = 0;
			int n = 0;
			for(int j = y - floor(ws/2); j < y + floor(ws/2)+1; ++j){
				sum += temp.getPixel(i, j);
				n += 1;
			}
			tgt.setPixel(i, y, round(sum/n));
		}
	}
}


void utility::colorbirght(image &src, image &tgt, int dr, int dg, int db){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			tgt.setPixel(i, j, channel::RED, checkValue(src.getPixel(i, j, channel::RED) * dr));
			
			tgt.setPixel(i, j, channel::BLUE, checkValue(src.getPixel(i, j, channel::BLUE) * db));
		
			tgt.setPixel(i, j, channel::GREEN, checkValue(src.getPixel(i, j, channel::GREEN) * dg));
		}
	}
}

void utility::colorvisual(image &src, image &tgt, int t, int v){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			if (abs(src.getPixel(i, j) - v) < t){
				tgt.setPixel(i, j, channel::RED, 255);
				tgt.setPixel(i, j, channel::BLUE, 0);
				tgt.setPixel(i, j, channel::GREEN, 0);
			} else {
				tgt.setPixel(i, j, channel::RED, src.getPixel(i, j, channel::RED));
				tgt.setPixel(i, j, channel::BLUE, src.getPixel(i, j, channel::BLUE));
				tgt.setPixel(i, j, channel::GREEN, src.getPixel(i, j, channel::GREEN));
			}
		}
	}
}

// new

void utility::histostretch(image &src, image &tgt, int A, int B){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt.copyImage(src);
	std::vector<int> pixels;
	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();
	
	// get min and max
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			int pixel = src.getPixel(i, j);
			pixels.push_back(pixel);
			min = (pixel < min) ? pixel : min;
			max = (pixel > max) ? pixel : max;
		}
	}

	// avoid divide by zero
	if (max-min == 0){
		max++;
	}	
	std::cout << "histostretch min: " << min << " max: " << max << std::endl;
	// stretch
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			int p_prime = checkValue((src.getPixel(i, j) - min) * ((B - A)/(max - min)) + A);
			tgt.setPixel(i, j, p_prime);
		}
	}
}

void utility::althistostretch(image &src, image &tgt, int A, int B){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	std::vector<int> pixels;
	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();
	
	// get min and max
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			int pixel = src.getPixel(i, j);
			pixels.push_back(pixel);
			min = (pixel < min) ? pixel : min;
			max = (pixel > max) ? pixel : max;
		}
	}

	min *= 1.05;
	max *= 0.95;

	// avoid divie by zero
	if (max-min == 0){
		max++;
	}

	// alt stretch
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			int p_prime;
			if (src.getPixel(i, j) < min){
				p_prime = checkValue(A);
			} else if (src.getPixel(i, j) > max){
				p_prime = checkValue(B);
			} else{
				p_prime = checkValue((src.getPixel(i, j) - min) * ((B - A)/(max - min)) + A);
			}
			tgt.setPixel(i, j, p_prime);
		}
	}	

}

void utility::histothres(image &src, image &tgt, int p1, char p2, int A, int B){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	image temp_bin(tgt);
	image temp_stretch(tgt);
	binarize(src, temp_bin, p1);
	histostretch(src, temp_stretch, A, B);

	int cond;
	if (p2 == 'B'){
		cond = MINRGB;  // background
	} else {
		cond = MAXRGB;  // foreground
	}

	for(int i=0; i < src.getNumberOfRows(); i++){
		for(int j=0; j<src.getNumberOfColumns(); j++){
			if (temp_bin.getPixel(i, j) == cond){
				tgt.setPixel(i, j, temp_stretch.getPixel(i, j));
			} else {
				tgt.setPixel(i, j, src.getPixel(i, j));
			}
		}
	}
}

void utility::colorstretch(image &src, image &tgt, char p1, int low, int high){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt.copyImage(src);

	channel c;
	switch(p1){
		case 'R':
			// std::cout << "channel is red\n";
			c = RED;
			break;
		case 'G':
			// std::cout << "channel is green\n";
			c = GREEN;
			break;
		case 'B':
			// std::cout << "channel is blue\n";
			c = BLUE;
			break;
		default:
			std::cout << "Invalid color channel for stretching\n";
			return;
	}

	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			int pixel = src.getPixel(i, j, c);
			min = (pixel < min) ? pixel : min;
			max = (pixel > max) ? pixel : max;
		}
	}

	// avoid divide by zero
	if (max-min == 0){
		max++;
	}

	std::cout << "colorstretch min: " << min << " max: " << max << std::endl;
	for(int i=0; i < src.getNumberOfRows(); i++){
		for(int j=0; j<src.getNumberOfColumns(); j++){
			int p_prime = checkValue((src.getPixel(i, j, c) - min) * ((high - low)/(max - min)) + low);
			tgt.setPixel(i, j, c, p_prime);
		}
	}
}

vector<int> utility::rgbtohsi(int R, int G, int B){
	double r, g, b;
	if (R+G+B == 0){
		// r = 0;
		// g = 0;
		// b = 0;
		return {0, 0, 0};
	}

	r = 1.0*R/(R+G+B);
	g = 1.0*G/(R+G+B);
	b = 1.0*B/(R+G+B);

	double i = (R + G + B)/(3.0*MAXRGB);
	double s = 1 - (3/(r+g+b)*std::min({r, g, b}));
	double h = std::acos( 0.5*((r-g)+(r-b))/std::pow((std::pow(r-g, 2) + ((r-b)*(g-b))), 0.5));;
	h = (b <= g) ? h : 2*3.14159265 - h;

	int H = std::lround(h * 180/3.14159265);
	int S = std::lround(s * 100);
	int I = std::lround(i * 255);
	return {H, S, I};
}

vector<int> utility::hsitorgb(int H, int S, int I){
	if (H + S + I == 0){
		return {0, 0, 0};
	}

	double r, g, b;
	double h = H * 3.14159265 / 180.0;
	double s = S / 100.0;
	double i = I / 255.0;

	if (H < 120 || H == 360){
		b = i*(1.0 - s);
		r = i*(1.0 + (s*std::cos(h)/std::cos(60 * 3.14159265 / 180.0 - h)));
		g = 3*i - (r+b);
	} else if (H >= 120 && H < 240){
		r = i*(1.0 - s);
		g = i*(1.0 + (s*std::cos(h - (120 * 3.14159265 / 180.0))/std::cos(60 * 3.14159265 / 180.0 - (h - (120 * 3.14159265 / 180.0)))));
		b = 3*i - (r+g);
	} else {
		g = i*(1.0 - s);
		b = i*(1.0 + (s*std::cos(h)/std::cos(60 * 3.14159265 / 180.0 - h)));
		r = 3*i - (g+b);
	}

	int R = std::lround(1.0*r*MAXRGB);
	int G = std::lround(1.0*g*MAXRGB);
	int B = std::lround(1.0*b*MAXRGB);

	return {R, G, B};
}
void utility::hsistretch(image &src, image &tgt, char HSI, int low, int high){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt.copyImage(src);

	vector<int> limits;
	int c;
	switch(HSI){
		case 'H':
			// std::cout << "channel is red\n";
			c = 0;
			limits.push_back(0);
			limits.push_back(360);
			break;
		case 'S':
			// std::cout << "channel is green\n";
			c = 1;
			limits.push_back(0);
			limits.push_back(100);
			break;
		case 'I':
			// std::cout << "channel is blue\n";
			c = 2;
			limits.push_back(0);
			limits.push_back(255);
			break;
		default:
			std::cout << "Invalid color channel for stretching\n";
			return;
	}

	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();
	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			vector<int> hsi = rgbtohsi(src.getPixel(i, j, RED), src.getPixel(i, j, GREEN), src.getPixel(i, j, BLUE));
			int pixel = hsi[c];
			min = (pixel < min) ? pixel : min;
			max = (pixel > max) ? pixel : max;
		}
	}

	// avoid divide by zero
	if (max-min == 0){
		max++;
	}

	std::cout << "hsi stretch min: " << min << " max: " << max << std::endl;
	// std::cout << "min: " << min << " max: " << max << std::endl;
	for(int i=0; i < src.getNumberOfRows(); i++){
		for(int j=0; j<src.getNumberOfColumns(); j++){
			vector<int> hsi = rgbtohsi(src.getPixel(i, j, RED), src.getPixel(i, j, GREEN), src.getPixel(i, j, BLUE));
			hsi[c] = hsiCheckValue((hsi[c] - min) * ((high - low)/(max - min)) + low, limits[0], limits[1]);
			vector<int> rgb = hsitorgb(hsi[0], hsi[1], hsi[2]);
			tgt.setPixel(i, j, RED, rgb[0]);
			tgt.setPixel(i, j, GREEN, rgb[1]);
			tgt.setPixel(i, j, BLUE, rgb[2]);
		}
	}
}
void utility::fullhsistretch(image &src, image &tgt, char RGB, int low, int high){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt.copyImage(src);

	vector<int> hsi = rgbtohsi(src.getPixel(0, 0, RED), src.getPixel(0, 0, GREEN), src.getPixel(0, 0, BLUE));
	int min_h = hsi[0];
	int min_s = hsi[1];
	int min_i = hsi[2];
	int max_h = hsi[0];
	int max_s = hsi[1];
	int max_i = hsi[2];

	for (int i=0; i<src.getNumberOfRows(); i++){
		for (int j=0; j<src.getNumberOfColumns(); j++){
			vector<int> hsi = rgbtohsi(src.getPixel(i, j, RED), src.getPixel(i, j, GREEN), src.getPixel(i, j, BLUE));
			min_h = (hsi[0] < min_h) ? hsi[0] : min_h;
			max_h = (hsi[0] > max_h) ? hsi[0] : max_h;
			min_s = (hsi[1] < min_s) ? hsi[1] : min_s;
			max_s = (hsi[1] > max_s) ? hsi[1] : max_s;
			min_i = (hsi[2] < min_i) ? hsi[2] : min_i;
			max_i = (hsi[2] > max_i) ? hsi[2] : max_i;
		}
	}

	// avoid divide by zero
	// if (max-min == 0){
	// 	max++;
	// }

	// std::cout << "min: " << min << " max: " << max << std::endl;
	for(int i=0; i < src.getNumberOfRows(); i++){
		for(int j=0; j<src.getNumberOfColumns(); j++){
			vector<int> hsi = rgbtohsi(src.getPixel(i, j, RED), src.getPixel(i, j, GREEN), src.getPixel(i, j, BLUE));
			hsi[0] = checkValue((hsi[0] - min_h) * ((high - low)/(max_h - min_h)) + low);
			hsi[1] = checkValue((hsi[1] - min_s) * ((high - low)/(max_s - min_s)) + low);
			hsi[2] = checkValue((hsi[2] - min_i) * ((high - low)/(max_i - min_i)) + low);
			vector<int> rgb = hsitorgb(hsi[0], hsi[1], hsi[2]);
			tgt.setPixel(i, j, RED, rgb[0]);
			tgt.setPixel(i, j, GREEN, rgb[1]);
			tgt.setPixel(i, j, BLUE, rgb[2]);
		}
	}
}

void utility::histoplot(image &src, image &tgt){ // do not use
	int plot_height = 500;
	int bar_thickness = 2;
	int plot_width = bar_thickness * 255;

	tgt.resize(plot_height, plot_width);
	vector<double> freq(256, 0);

	for(int i=0; i < src.getNumberOfRows(); i++){
		for(int j=0; j<src.getNumberOfColumns(); j++){
			freq[checkValue(src.getPixel(i, j))]++;
		}
	}

	// int max = std::numeric_limits<int>::min();
	// for(int i=0; i < freq.size(); ++i){
	// 	max = (freq[i] > max) ? freq[i] : max;
	// }

	auto max = std::max_element(std::begin(freq), std::end(freq));
	// for(int i=0; i<freq.size(); ++i){
	// 	freq[i] = std::round(freq[i]/(*max)*plot_hieght);
	// }

	// for(auto it : freq){
	// 	std::cout << it << " ";
	// }
	// std::cout << std::endl;

	for (int x=0; x<freq.size(); ++x){
		double norm_freq = std::round(freq[x]/(*max) * plot_height);
		for (int i=0; i<norm_freq; ++i){
			for (int j=0; j<bar_thickness; ++j){
				tgt.setPixel(plot_height - i, x + bar_thickness, 255);
			}
		}
	}
}

// new

void utility::sobel3(image &src, image &tgt){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	vector<vector<int>> hor_kernel{
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1}
	};

	vector<vector<int>> ver_kernel{
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1}
	};

	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();

	for(int i=1; i<src.getNumberOfRows()-1; ++i){
		for(int j=1; j<src.getNumberOfColumns()-1; ++j){
			int hor_cur = 0;
			int ver_cur = 0;
			for(int x=0; x<hor_kernel.size(); ++x){
				for(int y=0; y<hor_kernel[x].size(); ++y){
					hor_cur += src.getPixel(i+x-1, j+y-1) * hor_kernel[x][y];
					ver_cur += src.getPixel(i+x-1, j+y-1) * ver_kernel[x][y];
				}
			}
			double g = sqrt(pow(hor_cur, 2) + pow(ver_cur, 2));
			if (g < min) min=g;
			if (g > max) max=g;
			tgt.setPixel(i, j, round(g));
		}
	}

	// stretching
	for(int i=2; i<src.getNumberOfRows()-2; ++i){
		for(int j=2; j<src.getNumberOfColumns()-2; ++j){
			double stretch_g = (tgt.getPixel(i, j)-min) * 255/(max-min);
			tgt.setPixel(i, j, checkValue(round(stretch_g)));
		}
	}	
}

void utility::sobel5(image &src, image &tgt){
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	vector<vector<int>> hor_kernel{
		{-5, -8, -10, -8, -5},
		{-4, -10, -20, -10, -4},
		{0, 0, 0, 0, 0},
		{4, 10, 20, 10, 4},
		{5, 8, 10, 8, 5}
	};

	vector<vector<int>> ver_kernel{
		{-5, -4, 0, 4, 5},
		{-8, -10, 0, 10, 8},
		{-10, -20, 0, 20, 10},
		{-8, -10, 0, 10, 8},
		{-5, -4, 0, 4, 5}
	};

	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();

	for(int i=2; i<src.getNumberOfRows()-2; ++i){
		for(int j=2; j<src.getNumberOfColumns()-2; ++j){
			int hor_cur = 0;
			int ver_cur = 0;
			for(int x=0; x<hor_kernel.size(); ++x){
				for(int y=0; y<hor_kernel[x].size(); ++y){
					hor_cur += src.getPixel(i+x-2, j+y-2) * hor_kernel[x][y];
					ver_cur += src.getPixel(i+x-2, j+y-2) * ver_kernel[x][y];
				}
			}
			double g = sqrt(pow(hor_cur, 2) + pow(ver_cur, 2));
			if (g < min) min=g;
			if (g > max) max=g;
			tgt.setPixel(i, j, round(g));
		}
	}

	// stretching
	for(int i=2; i<src.getNumberOfRows()-2; ++i){
		for(int j=2; j<src.getNumberOfColumns()-2; ++j){
			double stretch_g = (tgt.getPixel(i, j)-min) * 255/(max-min);
			tgt.setPixel(i, j, checkValue(round(stretch_g)));
		}
	}
}

void utility::binaryedge(image &src, image &tgt, int TH, int TL, bool show_angles){

	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	vector<vector<int>> hor_kernel{
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1}
	};

	vector<vector<int>> ver_kernel{
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1}
	};

	vector<vector<double>> angles(src.getNumberOfRows(), vector<double>(src.getNumberOfColumns(), 0));

	int min = std::numeric_limits<int>::max();
	int max = std::numeric_limits<int>::min();

	for(int i=1; i<src.getNumberOfRows()-1; ++i){
		for(int j=1; j<src.getNumberOfColumns()-1; ++j){
			int hor_cur = 0;
			int ver_cur = 0;
			for(int x=0; x<hor_kernel.size(); ++x){
				for(int y=0; y<hor_kernel[x].size(); ++y){
					hor_cur += src.getPixel(i+x-1, j+y-1) * hor_kernel[x][y];
					ver_cur += src.getPixel(i+x-1, j+y-1) * ver_kernel[x][y];
				}
			}

			double g = sqrt(pow(hor_cur, 2) + pow(ver_cur, 2));
			if (g < min) min=g;
			if (g > max) max=g;
			tgt.setPixel(i, j, round(g));
			angles[i][j] = atan2(ver_cur, hor_cur);
		}
	}

	// stretch
	for(int i=1; i<src.getNumberOfRows()-1; ++i){
		for(int j=1; j<src.getNumberOfColumns()-1; ++j){
			double stretch_g = (tgt.getPixel(i, j)-min) * 255/(max-min);
			tgt.setPixel(i, j, checkValue(round(stretch_g)));
		}
	}

	// binarize
	for(int i=1; i<src.getNumberOfRows()-1; ++i){
		for(int j=1; j<src.getNumberOfColumns()-1; ++j){	
			int g = tgt.getPixel(i, j);
			if (g > TH){
				tgt.setPixel(i, j, 255);
			} else if (g < TL){
				tgt.setPixel(i, j, 0);
			} else {
				tgt.setPixel(i, j, 0);
				// check neighbors
				for(int x=-1; x<=1; ++x){
					for(int y=-1; y<=1; ++y){
						if (tgt.getPixel(i+x, j+y) > TH){
								tgt.setPixel(i, j, 255);
						}
					}
				}				
			}
		}		
	}

	// angles
	if (show_angles){
		int angle = 45;  // hardcode angle
		for(int i=1; i<src.getNumberOfRows()-1; ++i){
			for(int j=1; j<src.getNumberOfColumns()-1; ++j){
				int cur = round(angles[i][j] * 180/M_PI);
				if (!((cur >= angle - 10 && cur <= angle + 10) || ((cur >= (angle - 10 - 180)) && (cur <= angle + 10 -180)))){
					tgt.setPixel(i, j, 0);
				}
			}
		}		
	}

}

void utility::gausobel(image &src, image &tgt, double sigma){

	// create kernel
    int radius = (int) ceil(3 * sigma);
    int kernel_size = radius * 2 + 1;

    vector<double> kernel(kernel_size);
    float coeff = 1 / pow(2 * M_PI * pow(sigma, 2), 0.5);
    
    int j = 0;
    for(int i = -radius; i <= radius; i++)
    {
        kernel[j++] = coeff * exp(-pow(i, 2) / (2 * pow(sigma, 2)));
    }

	// run filter
	image temp(src);

	// 1st pass
	for(int i=radius; i<src.getNumberOfRows() - radius; ++i){
		for(int j=radius; j<src.getNumberOfColumns() - radius; ++j){
			double sum = 0;
			for (int k=0; k<kernel_size; ++k){
				sum += src.getPixel(i+k-radius, j)*kernel[k];
			}
			temp.setPixel(i, j, checkValue(round(sum)));
		}
	}

	// 2nd pass
	for(int i=radius; i<src.getNumberOfRows() - radius; ++i){
		for(int j=radius; j<src.getNumberOfColumns() - radius; ++j){
			double sum = 0;
			for (int k=0; k<kernel_size; ++k){
				sum += src.getPixel(i, j+k-radius)*kernel[k];
			}
			temp.setPixel(i, j, checkValue(round(sum)));
		}
	}

	// sobel
	sobel3(temp, tgt);
}

void utility::sobelcv(image &src, image &tgt, int win_size){
	char * TEMP = "temp.pgm";

	if (win_size%2 == 0){
		++win_size; // must be odd
	}
	if (win_size > 31) {
		win_size = 31; // opencv limit is 31
	}

	Mat image, image_gray;
  	Mat grad;

	src.save(TEMP);
	image = imread(TEMP, IMREAD_COLOR);
    cvtColor(image, image_gray, COLOR_BGR2GRAY);
	
	auto start = std::chrono::system_clock::now();
    
	Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;
    Sobel(image_gray, grad_x, CV_64F, 1, 0, win_size, -1, 0, BORDER_DEFAULT);
    Sobel(image_gray, grad_y, CV_64F, 0, 1, win_size, -1, 0, BORDER_DEFAULT);
    
    convertScaleAbs(grad_x, abs_grad_x);
    convertScaleAbs(grad_y, abs_grad_y);
    addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "sobelcv (without) " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	imwrite(TEMP, grad);
	tgt.read(TEMP);
	remove(TEMP);
}

void utility::cannycv(image &src, image &tgt){
	char * TEMP = "temp.pgm";

	Mat image, image_gray;
  	Mat grad;


	src.save(TEMP);
	image = imread(TEMP, IMREAD_COLOR);
    cvtColor(image, image_gray, COLOR_BGR2GRAY);

	auto start = std::chrono::system_clock::now();

    Canny(image_gray, grad, 0, 300);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "cannycv (without) " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";								

	imwrite(TEMP, grad);
	tgt.read(TEMP);
	remove(TEMP);
}

void utility::hsvcanny(image &src, image &tgt){
	char * TEMP = "temp.ppm";

	Mat image, image_hsv, hsv_channels[3];
  	Mat grad, grad_hsv, temp_hsv;

	src.save(TEMP);
	image = imread(TEMP, IMREAD_COLOR);

	auto start = std::chrono::system_clock::now();
    
	cvtColor(image, image_hsv, COLOR_BGR2HSV);
	split(image_hsv, hsv_channels);
    Canny(hsv_channels[1], hsv_channels[1], 0, 300); // v channel
	merge(hsv_channels, 3, grad_hsv);
	cvtColor(grad_hsv, grad, COLOR_HSV2BGR);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "hsvcanny (without) " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	imwrite(TEMP, grad);
	tgt.read(TEMP);
	remove(TEMP);
}

void utility::histostretchcv(image &src, image &tgt){
	Mat src_image, tgt_image;
	char * TEMP = "temp.pgm";
	src.save(TEMP);
	src_image = imread(TEMP, IMREAD_GRAYSCALE);

	// hist stretching
	normalize(src_image, tgt_image, 255, 0, NORM_MINMAX);

	imwrite(TEMP, tgt_image);
	tgt.read(TEMP);
	remove(TEMP);
}

void utility::histoequalcv(image &src, image &tgt){
	Mat src_image, tgt_image;
	char * TEMP = "temp.pgm";
	src.save(TEMP);
	src_image = imread(TEMP, IMREAD_GRAYSCALE);

	// hist equalization
	equalizeHist(src_image, tgt_image);

	imwrite(TEMP, tgt_image);
	tgt.read(TEMP);
	remove(TEMP);	
}

void utility::gausscv(image &src, image &tgt){
	Mat src_image, tgt_image;
	char * TEMP = "temp.ppm";
	src.save(TEMP);
	src_image = imread(TEMP);
	
	GaussianBlur(src_image, tgt_image, Size(9, 9), 0);

	imwrite(TEMP, tgt_image);
	tgt.read(TEMP);
	remove(TEMP);	
}

void utility::houghtrans(image &src, image &tgt){
	image original, blur, canny;
	Mat src_image, tgt_image;
	original.copyImage(src);

	// gaussian blur
	gausscv(src, blur);
	// canny edge 
	hsvcanny(blur, canny);
	
	char * TEMP = "temp.ppm";
    canny.save(TEMP);
	tgt_image = imread(TEMP);
	cvtColor(tgt_image, src_image, COLOR_BGR2GRAY);
	// hough trans
	vector<Vec3f> circles;
    HoughCircles(src_image, circles, HOUGH_GRADIENT, 1,
                 src_image.rows/16,  // change this value to detect circles with different distances to each other
                 100, 30, 2, 100 // change the last two parameters
            // (min_radius & max_radius) to detect larger circles
    );

	// uncomment to draw circles on original image
    original.save(TEMP);
	tgt_image = imread(TEMP);

	std::cout << "circles found: " << circles.size() << std::endl;
    for(int i = 0; i < circles.size(); i++) {
        Vec3i c = circles[i];
        Point center = Point(c[0], c[1]);
        int radius = c[2];

	    // draw circle outline
        circle(tgt_image, center, radius, Scalar(255,0,255), 1, cv::LINE_AA);

    //     // draw circle center
    //     circle(tgt_image, center, 1, Scalar(0,255,0), 3, cv::LINE_AA);
    
    }


	imwrite(TEMP, tgt_image);
	tgt.read(TEMP);
	remove(TEMP);
}

void utility::hesobel(image &src, image &tgt){
	image temp_tgt;
	histoequalcv(src, temp_tgt);
	sobelcv(temp_tgt, tgt, 3);
}

void utility::hecanny(image &src, image &tgt){
	image temp_tgt;
	histoequalcv(src, temp_tgt);
	cannycv(src, tgt);
}

void utility::qrdecode(image &src, image &tgt){
	// Read image
  	Mat src_image, tgt_image;
	char * TEMP = "temp.pgm";
	src.save(TEMP);
	src_image = imread(TEMP, IMREAD_GRAYSCALE); 

  	QRCodeDetector qrDecoder = QRCodeDetector();
	Mat bbox, rectifiedImage;
  	std::string data = qrDecoder.detectAndDecode(src_image, bbox, rectifiedImage);
	std::cout << "qr code: " << data << std::endl;
	tgt.copyImage(src);
}

