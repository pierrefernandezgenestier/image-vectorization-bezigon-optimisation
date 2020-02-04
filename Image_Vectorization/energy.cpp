#include "Energy.h"

using namespace cv;
using namespace std;
using namespace Eigen;

double energy_spt(Bezigon B, int j) {
	double energy = 0.;
	vector<array<double, 2>> intersections;
	Bezier bez = B.get_bezier(j);
	intersections = self_intersect(bez);
	if (intersections.size() != 0) energy += B.get_length_gl(j, intersections[0][0], intersections[0][1]);

	Bezier bez_prev = B.get_bezier(j - 1);
	intersections = intersect(bez_prev, bez);
	if (intersections.size() == 1) {
		energy += B.get_length_gl(j - 1, intersections[0][0], 1);
		energy += B.get_length_gl(j, 0, intersections[0][1]);
	}
	if (intersections.size() == 2) {
		energy += B.get_length_gl(j - 1, intersections[0][0], intersections[1][0]);
		energy += B.get_length_gl(j, intersections[0][1], intersections[1][1]);
	}

	Bezier bez_next = B.get_bezier(j + 1);
	intersections = intersect(bez, bez_next);
	if (intersections.size() == 1) {
		energy += B.get_length_gl(j, intersections[0][0], 1);
		energy += B.get_length_gl(j + 1, 0, intersections[0][1]);
	}
	if (intersections.size() == 2) {
		energy += B.get_length_gl(j, intersections[0][0], intersections[1][0]);
		energy += B.get_length_gl(j + 1, intersections[0][1], intersections[1][1]);
	}

	return energy;
}

double energy_spt(Bezigon B) {
	double energy = 0.0;
	for (int j = 0; j < B.Bx.rows(); j++) {
		energy = energy + energy_spt(B, j);
	}
	return energy;
}

double energy_apt(Bezigon B, int j) {
	array<Point2d, 2> tangents_j = B.get_tangent(j);
	Point2d a_j = tangents_j[0];
	Point2d b_j = tangents_j[1];
	double angle = acos(a_j.dot(b_j) / (norm(a_j) * norm(b_j)));
	return (angle != 0 && !isnan(angle)) ? angle : 0.;
}

double energy_apt(Bezigon B) {
	double energy = 0.0;
	for (int j = 0; j < B.Bx.rows(); j++) {
		energy = energy + energy_apt(B, j);
	}
	return energy;
}

//double energy_angles(Bezigon B, int point) {
//	// 3 control points are influenced when we optimize 5 points.
//	return energy_angles_ctrl_pt(B, point) + energy_angles_ctrl_pt(B, point + 1) + energy_angles_ctrl_pt(B, point + 2);
//};

double energy_hpt(Bezigon B, int j) {
	array<Point2d, 2> tangents_j = B.get_tangent(j);
	return (norm(tangents_j[0]) > 0 && norm(tangents_j[1]) > 0)
		? (1 / norm(tangents_j[0]) + 1 / norm(tangents_j[1]))
		: 0.;
}

double energy_hpt(Bezigon B) {
	double energy = 0.0;
	for (int j = 0; j < B.Bx.rows(); j++) {
		energy = energy + energy_hpt(B, j);
	}
	return energy;
};

//double energy_bezier_handles(Bezigon B, int point) {
//	// 3 control points are influenced when we optimize 5 points.
//	return energy_bezier_handles_ctrl_pt(B, point) + energy_bezier_handles_ctrl_pt(B, point + 1) + energy_bezier_handles_ctrl_pt(B, point + 2);
//};

double energy_lpt(Bezigon B, int j) {
	return B.get_length_gl(j);
}

double energy_lpt(Bezigon B) {
	return B.get_length();
}

double energy_data(VectorizationData vd) {
	double lo = vd.B.lo;
	double energy = 0.0;
	Image<Vec3b> rasterized = get_rasterized(vd);
	for (size_t xx = 0; xx < vd.I.width(); xx++) {
		for (size_t yy = 0; yy < vd.I.height(); yy++) {
			Vec3b ras = rasterized(xx, yy);
			Vec3b im = vd.I(xx, yy);
			energy += norm(ras - im);
		}
	};
	return energy / lo;
};

//double energy_partial(VectorizationData vd, int j) {
//	double energy_prior = energy_bezier_handles(vd.B, j);
//	//energy_angles(vd.B, j) ++energy_curve_length(vd.B,j);
//	double energy = energy_data(vd);
//	//cout << "angles"<<lambda_angles * energy_angles(vd.B, j)<<endl;
//	//cout << "hangles "<<lambda_handles * energy_bezier_handles(vd.B, j)<<endl;
//	//cout << "data"<<energy<<endl;
//	//return energy_curve_length(vd.B,j);
//	//return energy;
//	//return energy_bezier_handles(vd.B,j);
//	//return  energy_angles(vd.B, j);
//	return energy_prior + energy;
//
//};

array<double, 2> energy_to_minimize(VectorizationData vd, int j, array<double, 10> vals_inp) {
	MatrixXd Bx_prime = vd.B.Bx;
	MatrixXd By_prime = vd.B.By;
	Bezigon B_prime = Bezigon(Bx_prime, By_prime);
	B_prime.update(vals_inp, j);
	VectorizationData vd_prime = VectorizationData(B_prime, vd.I);
	double e_data = energy_data(vd_prime);
	//double e_data = 0;
	double e_prior = lambda_spt * (energy_spt(vd_prime.B, j - 1) + energy_spt(vd_prime.B, j) + (energy_spt(vd_prime.B, j + 1)))
		+ lambda_apt * (energy_apt(vd_prime.B, j) + energy_apt(vd_prime.B, j + 1) + energy_apt(vd_prime.B, j + 2))
		+ lambda_hpt * (energy_hpt(vd_prime.B, j) + energy_hpt(vd_prime.B, j + 1) + energy_hpt(vd_prime.B, j + 2))
		+ lambda_lpt * (energy_lpt(vd_prime.B, j) + energy_lpt(vd_prime.B, j + 1));
	return { lambda_data * e_data, e_prior };
}

//double energy_to_minimize_prior(VectorizationData vd, int j, array<double, 10> vals_inp) {
//	MatrixXd Bx_prime = vd.B.Bx;
//	MatrixXd By_prime = vd.B.By;
//	Bezigon B_prime = Bezigon(Bx_prime, By_prime);
//	B_prime.update(vals_inp, j);
//	VectorizationData vd_prime = VectorizationData(B_prime, vd.I);
//	double E = lambda_spt * (energy_spt(vd_prime.B, j))
//		+ lambda_apt * (energy_apt(vd_prime.B, j) + energy_apt(vd_prime.B, j + 1) + energy_apt(vd_prime.B, j + 2))
//		+ lambda_hpt * (energy_hpt(vd_prime.B, j) + energy_hpt(vd_prime.B, j + 1))
//		+ lambda_lpt * (energy_lpt(vd_prime.B, j) + energy_lpt(vd_prime.B, j + 1));
//	return E;
//}

