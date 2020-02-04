#include "Utils.h"

using namespace cv;
using namespace std;

Point2d barycenter(double t, Point2d A, Point2d B) {
	Point2d p;
	p.x = (1 - t) * (float)A.x + t * (float)B.x;
	p.y = (1 - t) * (float)A.y + t * (float)B.y;
	return p;
};

double distance(Point2d A, Point2d B) {
	return norm(Mat(A), Mat(B));
}

double angle(Point2d A, Point2d B, Point2d C) {
	if (A == B || B == C) return 0;
	return acos((C - B).dot(A - B) / (distance(A, B) * distance(C, B)));
};

double det(Point2d A, Point2d B) {
	return A.x * B.y - A.x * B.y;
}

double overlap(array<Point2d, 2> bbox1, array<Point2d, 2> bbox2) {
	Point2d corner_min1 = bbox1[0];
	Point2d corner_max1 = bbox1[1];
	Point2d corner_min2 = bbox2[0];
	Point2d corner_max2 = bbox2[1];
	return MAX(0, (MIN(corner_max1.x, corner_max2.x) - MAX(corner_min1.x, corner_min2.x)))
		* MAX(0, (MIN(corner_max1.y, corner_max2.y) - MAX(corner_min1.y, corner_min2.y)));
}

/***********************
* Credits to ufoymtest *
***********************/

/********************************************************
*                                                       *
* This function determines if a double is small enough  *
* to be zero. The purpose of the subroutine is to try   *
* to overcome precision problems in math routines.      *
*                                                       *
********************************************************/

static int isZero(double x) {
	return x > -EQN_EPS && x < EQN_EPS;
}


int solveLinear(double c1, double c0, double& s0) {
	if (isZero(c1)) return 0;
	s0 = -c0 / c1;
	return 1;
}

/********************************************************
*                                                       *
* This function determines the roots of a quadric       *
* equation.                                             *
* It takes as parameters a pointer to the three         *
* coefficient of the quadric equation (the c[2] is the  *
* coefficient of x2 and so on) and a pointer to the     *
* two element array in which the roots are to be        *
* placed.                                               *
* It outputs the number of roots found.                 *
*                                                       *
********************************************************/

int solveQuadric(double c2, double c1, double c0, double& s0, double& s1)
{
	double p, q, D;
	if (isZero(c2)) return solveLinear(c1, c0, s0);

	p = c1 / (2.0 * c2);
	q = c0 / c2;
	D = p * p - q;

	if (isZero(D)) {
		s0 = s1 = -p;
		return 1;
	}
	if (D < 0.0) return 0;

	else {
		double sqrt_D = sqrt(D);
		s0 = sqrt_D - p;
		s1 = -sqrt_D - p;
		return 2;
	}
}

/********************************************************
*                                                       *
* This function determines the roots of a cubic         *
* equation.                                             *
* It takes as parameters a pointer to the four          *
* coefficient of the cubic equation (the c[3] is the    *
* coefficient of x3 and so on) and a pointer to the     *
* three element array in which the roots are to be      *
* placed.                                               *
* It outputs the number of roots found                  *
*                                                       *
********************************************************/

int solveCubic(double c3, double c2, double c1, double c0, double& s0, double& s1, double& s2) {
	int i, num;
	double  sub,
		A, B, C,
		sq_A, p, q,
		cb_p, D;

	if (isZero(c3)) return solveQuadric(c2, c1, c0, s0, s1);

	A = c2 / c3;
	B = c1 / c3;
	C = c0 / c3;

	sq_A = A * A;
	p = 1.0 / 3.0 * (-1.0 / 3.0 * sq_A + B);
	q = 1.0 / 2.0 * (2.0 / 27.0 * A * sq_A - 1.0 / 3.0 * A * B + C);

	cb_p = p * p * p;
	D = q * q + cb_p;

	if (isZero(D)) {
		if (isZero(q)) {
			s0 = 0.0;
			num = 1;
		}
		else {
			double u = cbrt(-q);
			s0 = 2.0 * u;
			s1 = -u;
			num = 2;
		}
	}
	else if (D < 0.0) {
		double phi = 1.0 / 3.0 * acos(-q / sqrt(-cb_p));
		double t = 2.0 * sqrt(-p);
		s0 = t * cos(phi);
		s1 = -t * cos(phi + M_PI / 3.0);
		s2 = -t * cos(phi - M_PI / 3.0);
		num = 3;
	}
	else {
		double sqrt_D = sqrt(D);
		double u = cbrt(sqrt_D + fabs(q));
		if (q > 0.0) s0 = -u + p / u;
		else s0 = u - p / u;
		num = 1;
	}

	sub = 1.0 / 3.0 * A;
	s0 -= sub;
	s1 -= sub;
	s2 -= sub;
	return num;
}

/********************************************************
*                                                       *
* end of credits										*
*                                                       *
********************************************************/

vector<double> get_roots(double c2, double c1, double c0) {
	vector<double> roots;
	double s0, s1;
	int num = solveQuadric(c2, c1, c0, s0, s1);
	if (num == 0) return roots;
	else {
		roots.push_back(s0);
		roots.push_back(s1);
	}
	return roots;
}

vector<double> get_roots(double c3, double c2, double c1, double c0) {
	vector<double> roots;
	double s0, s1, s2;
	int num = solveCubic(c3, c2, c1, c0, s0, s1, s2);
	if (num == 1) roots.push_back(s0);
	else {
		roots.push_back(s0);
		roots.push_back(s1);
		roots.push_back(s2);
	}
	return roots;
}