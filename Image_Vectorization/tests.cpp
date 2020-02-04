#include "Tests.h"

using namespace cv;
using namespace Eigen;


void test_intersect_bez() {
	std::cout << ">>> TEST Bezier intersections" << std::endl;
	Bezier bez1(Point2d(10, 100), Point2d(90, 30), Point2d(40, 140), Point2d(220, 220));
	Bezier bez2(Point2d(5, 150), Point2d(180, 20), Point2d(80, 250), Point2d(210, 190));
	vector<array<double, 2>> vector_ts = intersect(bez1, bez2);
	for (auto ts : vector_ts) {
		cout << "t bez 1 : " << ts[0] << " , t bez 2 : " << ts[1] << endl;
	}
}

void test_compute_angles(VectorizationData vd) {
	std::cout << " >>> TEST Angles" << std::endl;
	for (int j = 0; j < vd.B.Bx.rows(); j++) {
		cout << energy_apt(vd.B, j);
	}
	std::cout << std::endl;
};

void test_rms(VectorizationData vd) {
	cout << ">>> TEST rms" << endl;
	rms(vd, "test_rms");
	cout << endl;
};

//
//
//void test_compute_energy_angles(VectorizationData vd) {
//	Energy* E = new Energy();
//	std::cout << "TEST The Angles Energy is " << E->energy_angles(*vd.B) << std::endl;
//	std::cout << std::endl;
//
//};
//
//
//void test_compute_energy_data(VectorizationData vd) {
//	Energy* E = new Energy();
//	std::cout << "TEST The Energy data is " << E->energy_data(vd) << std::endl;
//	std::cout << std::endl;
//}
//
//void test_compute_bezier_handles(VectorizationData vd) {
//	Energy* E = new Energy();
//	std::cout << "TEST The Bezier Handles Energy is " << E->energy_bezier_handles(*vd.B) << std::endl;
//	std::cout << std::endl;
//};
//
//void test_espilon(VectorizationData vd){
//	Propagation * p = new Propagation(vd);
//	cout << ">>> testing epsilon" << endl;
//	for (int i = 0; i <10; i++){
//		p->propagate(2, 0.2, pow(10,-i), " eps" + to_string(i));
//	}
//
//}
//
//void test_alpha(VectorizationData vd){
//	Propagation * p = new Propagation(vd);
//	cout << ">>> testing alpha" << endl;
//	for (int i = 2; i <5; i++){
//		p->propagate(3, pow(10,i), 1, " alpha 10^^" + to_string(i));
//	}
//
//}
