#pragma once

#ifndef ENTITIES_H_INCLUDED
#define ENTITIES_H_INCLUDED

#include <iostream>
#include <ostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>
#include <map> 

#include "Image.h"
#include "Utils.h"

class Bezier
{
public:
	Bezier();
	Bezier(cv::Point2d p0, cv::Point2d p1, cv::Point2d p2, cv::Point2d p3);
	cv::Point2d control_points[4];
	cv::Point2d cubic_interpolation(double t);
	std::array<Bezier, 2> subdivide(double t0);
	Bezier subdivide(double t0, double t1);
	std::array<cv::Point2d, 2> get_bounding_box();
};

class Bezigon
{
public:
	Eigen::MatrixXd Bx;
	Eigen::MatrixXd By;
	Vec3b C;
	double lo;
	Bezigon();
	Bezigon(Eigen::MatrixXd _Bx, Eigen::MatrixXd _By, Vec3b _C = Vec3b(0, 0, 0));
	Bezigon(std::vector<cv::Point2d> vector_points, Vec3b _C = Vec3b(0, 0, 0), int method = 2);

	//SET 
	void set_point_x(int j, int i, double coord_x);
	void set_point_y(int j, int i, double coord_y);

	//GET
	Bezier get_bezier(int j);
	cv::Point2d get_pt(int j, int i);
	array<Point2d, 2> get_tangent(int j);
	//COMPUTE
	double get_length(int j, double t0 = 0., double t1 = 1.);
	double get_length_gl(int j, double t0 = 0., double t1 = 1.);
	double get_length();

	//PROPAGATION FUNCTIONS
	void update(std::array<double, 10> vals_inp, int j);
	std::array<double, 10> input_propagation(int j);

	//PLOT
	void plot_curve(Image<Vec3b> I, std::string nom);
};

struct VectorizationData
{
	Bezigon B;
	Image<Vec3b> I;
	VectorizationData(Bezigon _B, Image<Vec3b> _I);
};

int intersect(Bezier bez1, Bezier bez2, std::vector<std::array<double, 2>>* vector_ts, std::array<double, 2> ts, int n_rec);
std::vector<std::array<double, 2>> intersect(Bezier bez1, Bezier bez2);
std::vector<std::array<double, 2>> self_intersect(Bezier bez);

bool to_svg(vector<VectorizationData> vector_vd, string name);

#endif