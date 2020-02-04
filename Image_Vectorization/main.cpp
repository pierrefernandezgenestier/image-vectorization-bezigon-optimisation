#include "main.h"

using namespace cv;
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
	//int nb_etapes = atoi(argv[1]);
	//float alpha = atof(argv[2]);
	//float epsilon = atof(argv[3]);
	//string nom = string(argv[4]);
	//string path = string(argv[5]);
	int nb_etapes = 10;
	double alpha = 15.; //learning rate
	double epsilon = 1.0;
	string name = "test";
	string path = "../data/heart.png";

	const Image<Vec3b> I = imread(path);
	cout << ">>> Initializing : CLIK on edge points. PRESS 0" << endl;
	VectorizationData vd = initialisation(I);
	propagate(&vd, nb_etapes, alpha, epsilon, name);
	cout << endl;

	waitKey();
	return 0;
};

