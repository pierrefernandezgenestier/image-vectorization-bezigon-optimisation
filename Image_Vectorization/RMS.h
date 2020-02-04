#pragma once

#ifndef RMS_H_INCLUDED
#define RMS_H_INCLUDED

#include <iostream>
#include <ostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "Entities.h"
#include "Image.h"
#include "Utils.h"

int phi(double t);
int psi(double t);
int phi_s_k(double t, int s, int k);
int psi_s_k(double t, int s, int k);
Eigen::Matrix2i haar(int s, int kx, int ky, double x, double y);
cv::Point2d normalize(cv::Point2d p, int wh);
Bezigon normalize(Bezigon B, int wh);
vector<Bezier> clip(Bezier bez, int left, int right, int bottom, int top);
array<cv::Point2d, 2> get_KL(Bezier bez);
Bezier transform(Bezier bez, cv::Point Q, int ss, int kx, int ky);
array<double, 3> get_coef(Bezigon B, int jj, int ss, int kx, int ky);
double c_00(Bezier bez);
double color(Bezigon B, int d, int N, map<array<int, 4>, array<double, 3>> coef_map, cv::Point2d p);
Image<Vec3b> get_rasterized_wavelet(VectorizationData vd, int d, int N, map<array<int, 4>, array<double, 3>> coef_map, int wh);
void rms_wavelet(VectorizationData vd);

bool is_interior(cv::Point, Bezigon B);
Image<cv::Vec3b> get_rasterized(VectorizationData vd);
void rms(VectorizationData vd, string name);

bool is_interior(Point p, std::map<int, std::vector<double>> map_intersect);
Image<cv::Vec3b> get_rasterized_fast(VectorizationData vd);

#endif