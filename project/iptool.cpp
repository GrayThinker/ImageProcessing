/************************************************************
 *															*
 * This sample project include three functions:				*
 * 1. Add intensity for gray-level image.					*
 *    Input: source image, output image name, value			*
 *															*
 * 2. Image thresholding: pixels will become black if the	*
 *    intensity is below the threshold, and white if above	*
 *    or equal the threhold.								*
 *    Input: source image, output image name, threshold		*
 *															*
 * 3. Image scaling: reduction/expansion of 2 for 			*
 *    the width and length. This project uses averaging 	*
 *    technique for reduction and pixel replication			*
 *    technique for expansion.								*
 *    Input: source image, output image name, scale factor	*
 *															*
 ************************************************************/

#include "../iptools/core.h"
#include <strings.h>
#include <string.h>
#include <vector>
#include <chrono>
#include <ctime>

using namespace std;

#define MAXLEN 256

int main (int argc, char** argv)
{
	image src, out;
	FILE *fp;
	char str[MAXLEN];
	char outfile[MAXLEN];
	char *pch;
	if ((fp = fopen(argv[1],"r")) == NULL) {
		fprintf(stderr, "Can't open file: %s\n", argv[1]);
		exit(1);
	}
	cout << "opened file\n";

	
	while(fgets(str,MAXLEN,fp) != NULL) {
		pch = strtok(str, " ");
		src.read(pch);
		out.read(pch);
		pch = strtok(NULL, " ");
		strcpy(outfile, pch);

		cout << "Done making copy\n";
		cout << "file dim: " << src.getNumberOfRows() << " x " << src.getNumberOfColumns() << "\n";
		// # ROI
		pch = strtok(NULL, " ");
		int no_roi = atoi(pch);

		if (no_roi < 1 || no_roi > 3){
			printf("Invalid #ROIs %s\n", pch);
			break;
		}
		vector<image*> rois;
		vector<vector<int>> roi_dims;

		for(int counter=0; counter < no_roi; ++counter){
			pch = strtok(NULL, " ");
			int x = atoi(pch);
			pch = strtok(NULL, " ");
			int y = atoi(pch);
			pch = strtok(NULL, " ");
			int sx = atoi(pch);
			pch = strtok(NULL, " ");
			int sy = atoi(pch);
			// cout << "x:" << x << " y: " << y << " sx: " << sx << " sy: " << sy << "\n";
			if (x + sx > src.getNumberOfColumns() || x + sx - 1 < 0){
				cout << "Invalid ROI x dimension\n";
				break;
			}
			if (y + sy > src.getNumberOfRows() || y + sy - 1 < 0){
				cout << "Invalid ROI y dimension\n";
				break;
			}


			// copy part of src image to roi
			roi_dims.insert(roi_dims.begin(), {x, y, sx, sy});
			image roi, tgt;
			roi.resize(sy, sx);
			for(int i=0; i< sy; ++i){
				for(int j=0; j< sx; ++j){
					roi.setPixel(i, j, RED, src.getPixel(y+i, x+j, RED));
					roi.setPixel(i, j, GREEN, src.getPixel(y+i, x+j, GREEN));
					roi.setPixel(i, j, BLUE, src.getPixel(y+i, x+j, BLUE));
				}
			}


			// functions
			pch = strtok(NULL, " ");
			if (strncasecmp(pch,"add",MAXLEN)==0) {
				/* Add Intensity */
				pch = strtok(NULL, " ");
					utility::addGrey(roi,tgt,atoi(pch));
			}

			else if (strncasecmp(pch,"binarize",MAXLEN)==0) {
				/* Thresholding */
				pch = strtok(NULL, " ");
				utility::binarize(roi,tgt,atoi(pch));
			}

			else if (strncasecmp(pch,"scale",MAXLEN)==0) {
				/* Image scaling */
				pch = strtok(NULL, " ");
				utility::scale(roi,tgt,atof(pch));
			}

			else if (strncasecmp(pch, "dualthres", MAXLEN)==0) {
				/* Dual thresholding */
				char *pch1 = strtok(NULL, " ");
				char *pch2 = strtok(NULL, " ");
				char *pch3 = strtok(NULL, " ");

				utility::dualthres(roi, tgt, atoi(pch1), atoi(pch2), atoi(pch3));
			}

			else if (strncasecmp(pch, "reg2dsmooth", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				auto start = std::chrono::system_clock::now();
    			utility::reg2dsmooth(roi, tgt, atoi(p1));

				// timing
				// auto end = std::chrono::system_clock::now();
				// std::chrono::duration<double> elapsed_seconds = end-start;
				// std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				// std::cout << "reg2d " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
				
			}

			else if (strncasecmp(pch, "sep2dsmooth", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				auto start = std::chrono::system_clock::now();				
				utility::sep2dsmooth(roi, tgt, atoi(p1));

				// timing
				// auto end = std::chrono::system_clock::now();
				// std::chrono::duration<double> elapsed_seconds = end-start;
				// std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				// std::cout << "sep2d " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s\n";	
			}

			else if (strncasecmp(pch, "colorbright", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				char *p2 = strtok(NULL, " ");
				char *p3 = strtok(NULL, " ");
				utility::colorbirght(roi, tgt, atoi(p1), atoi(p2), atoi(p3));
			}

			else if (strncasecmp(pch, "colorvisual", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				char *p2 = strtok(NULL, " ");
				utility::colorvisual(roi, tgt, atoi(p1), atoi(p2));
			}

// new
			else if (strncasecmp(pch, "histostretch", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				char *p2 = strtok(NULL, " ");
				utility::histostretch(roi, tgt, atoi(p1), atoi(p2));

			}

			else if (strncasecmp(pch, "althistostretch", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				char *p2 = strtok(NULL, " ");
				utility::althistostretch(roi, tgt, atoi(p1), atoi(p2));

			}

			else if (strncasecmp(pch, "histothres", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				char *p2 = strtok(NULL, " ");  // B/F
				char *p3 = strtok(NULL, " ");  // A
				char *p4 = strtok(NULL, " ");  // B
				utility::histothres(roi, tgt, atoi(p1), *p2, atoi(p3), atoi(p4));
			}

			else if (strncasecmp(pch, "colorstretch", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");  // RGB
				char *p2 = strtok(NULL, " ");  // low
				char *p3 = strtok(NULL, " ");  // high
				utility::colorstretch(roi, tgt, *p1, atoi(p2), atoi(p3));
			}

			// for testing only
			// else if (strncasecmp(pch, "rgbtohsi", MAXLEN)==0){
			// 	// std::cout << "RGB to HSI cannot be used by parameters\n";
			// 	std::cout << "RGB: 255, 255, 255 HSI: " << utility::rgbtohsi(255, 255, 255)[0] << ", " << utility::rgbtohsi(255, 255, 255)[1] << ", " << utility::rgbtohsi(255, 255, 255)[2] << std::endl;
			// 	std::cout << "RGB: 0, 0, 0 HSI: " << utility::rgbtohsi(0, 0, 0)[0] << ", " << utility::rgbtohsi(0, 0, 0)[1] << ", " << utility::rgbtohsi(0, 0, 0)[2] << std::endl;
			// 	std::cout << "RGB: 255, 0, 0 HSI: " << utility::rgbtohsi(255, 0, 0)[0] << ", " << utility::rgbtohsi(255, 0, 0)[1] << ", " << utility::rgbtohsi(255, 0, 0)[2] << std::endl;
			// 	std::cout << "RGB: 107, 200, 10 HSI: " << utility::rgbtohsi(107, 200, 10)[0] << ", " << utility::rgbtohsi(107, 200, 10)[1] << ", " << utility::rgbtohsi(107, 200, 10)[2] << std::endl;
			// 	std::cout << "RGB: 24, 98, 118 HSI: " << utility::rgbtohsi(24, 98, 118)[0] << ", " << utility::rgbtohsi(24, 98, 118)[1] << ", " << utility::rgbtohsi(24, 98, 118)[2] << std::endl;

			// 	std::cout << "AND BACK\n";

			// 	std::cout << "HSI: 0, 0, 255 RGB: " << utility::hsitorgb(0, 0, 255)[0] << ", " << utility::hsitorgb(0, 0, 255)[1] << ", " << utility::hsitorgb(0, 0, 255)[2] << std::endl;
			// 	std::cout << "HSI: 0, 0, 0 RGB: " << utility::hsitorgb(0, 0, 0)[0] << ", " << utility::hsitorgb(0, 0, 0)[1] << ", " << utility::hsitorgb(0, 0, 0)[2] << std::endl;
			// 	std::cout << "HSI: 0, 100, 85 RGB: " << utility::hsitorgb(0, 100, 85)[0] << ", " << utility::hsitorgb(0, 100, 85)[1] << ", " << utility::hsitorgb(0, 100, 85)[2] << std::endl;
			// 	std::cout << "HSI: 89, 91, 106 RGB: " << utility::hsitorgb(89, 91, 106)[0] << ", " << utility::hsitorgb(89, 91, 106)[1] << ", " << utility::hsitorgb(89, 91, 106)[2] << std::endl;
			// 	std::cout << "HSI: 192, 70, 80 RGB: " << utility::hsitorgb(192, 70, 80)[0] << ", " << utility::hsitorgb(192, 70, 80)[1] << ", " << utility::hsitorgb(192, 70, 80)[2] << std::endl;
			// 	break;
			// }

			else if (strncasecmp(pch, "hsistretch", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");  // HSI
				char *p2 = strtok(NULL, " ");  // low
				char *p3 = strtok(NULL, " ");  // high
				utility::hsistretch(roi, tgt, *p1, atoi(p2), atoi(p3));
			}

			else if (strncasecmp(pch, "fullhsistretch", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");  // HSI
				char *p2 = strtok(NULL, " ");  // low
				char *p3 = strtok(NULL, " ");  // high
				utility::fullhsistretch(roi, tgt, *p1, atoi(p2), atoi(p3));
			}			

			else if (strncasecmp(pch, "histoplot", MAXLEN)==0){
				utility::histoplot(roi, tgt);
			}

			else if(strncasecmp(pch, "sobel3", MAXLEN)==0 || strncasecmp(pch, "sobel3\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::sobel3(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "sobel3 " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "sobel5", MAXLEN)==0 || strncasecmp(pch, "sobel5\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::sobel5(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "sobel5 " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "binaryedge", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				char *p2 = strtok(NULL, " ");
				auto start = std::chrono::system_clock::now();
				utility::binaryedge(roi, tgt, atoi(p1), atoi(p2), false);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "binaryedge " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
				
				// To show angles
				utility::binaryedge(roi, tgt, atoi(p1), atoi(p2), true);
			}

			else if(strncasecmp(pch, "gausobel", MAXLEN)==0){
				char *p1 = strtok(NULL, " ");
				auto start = std::chrono::system_clock::now();
				utility::gausobel(roi, tgt, atof(p1));
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "gausobel " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";				
			}

			else if(strncasecmp(pch, "sobelcv", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				char *p1 = strtok(NULL, " ");
				utility::sobelcv(roi, tgt, atoi(p1));
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "sobelcv " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";								
			}

			else if(strncasecmp(pch, "cannycv", strlen("cannycv"))==0 || strncasecmp(pch, "cannycv\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::cannycv(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "cannycv " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";								
			}

			else if(strncasecmp(pch, "hsvcanny", strlen("hsvcanny"))==0 || strncasecmp(pch, "hsvcanny\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::hsvcanny(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "hsvcanny " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";												
			}
			
			// new

			else if(strncasecmp(pch, "histostretchcv", strlen("histostretchcv"))==0 || strncasecmp(pch, "histostretchcv\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::histostretchcv(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "histostretchcv " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "histoequalcv", strlen("histoequalcv"))==0 || strncasecmp(pch, "histoequalcv\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::histoequalcv(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "histoequalcv " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}
			
			else if(strncasecmp(pch, "gausscv", strlen("gausscv"))==0 || strncasecmp(pch, "gausscv\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::gausscv(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "gausscv " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "houghtrans", strlen("houghtrans"))==0 || strncasecmp(pch, "houghtrans\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::houghtrans(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "houghtrans " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "hesobel", strlen("hesobel"))==0 || strncasecmp(pch, "hesobel\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::hesobel(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "hesobel " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "hecanny", strlen("hecanny"))==0 || strncasecmp(pch, "hecanny\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::hecanny(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "hecanny " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else if(strncasecmp(pch, "qrdecode", strlen("qrdecode"))==0 || strncasecmp(pch, "qrdecode\n", MAXLEN)==0){
				auto start = std::chrono::system_clock::now();
				utility::qrdecode(roi, tgt);
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "qrdecode " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
			}

			else {
				printf("No such function: '%s'\n", pch);
				continue;
			}
			
			cout << "combining roi\n";
			rois.insert(rois.begin(), new image(tgt));
			if (counter == no_roi-1){
				for(int n=0; n<rois.size(); ++n){
					for(int i=0; i<roi_dims[n][3]; ++i){  // rows
						for(int j=0; j<roi_dims[n][2]; ++j){  // cols
							out.setPixel(i+roi_dims[n][1], j+roi_dims[n][0], RED, rois[n]->getPixel(i, j, RED));
							out.setPixel(i+roi_dims[n][1], j+roi_dims[n][0], GREEN, rois[n]->getPixel(i, j, GREEN));
							out.setPixel(i+roi_dims[n][1], j+roi_dims[n][0], BLUE, rois[n]->getPixel(i, j, BLUE));
						}
					}
				}
			}
			cout << "rois combined\n\n";
		}
		out.save(outfile);
	}
	fclose(fp);
	return 0;
}
