#pragma once

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <ostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <math.h>
#define EQN_EPS 1e-9

cv::Point2d barycenter(double t, cv::Point2d A, cv::Point2d B);
double distance(cv::Point2d A, cv::Point2d B);
double angle(cv::Point2d A, cv::Point2d B, cv::Point2d C);
double det(cv::Point2d A, cv::Point2d B);
double overlap(std::array<cv::Point2d, 2> bbox1, std::array<cv::Point2d, 2> bbox2);


static int isZero(double x);
int solveLinear(double c1, double c0, double& s0);
int solveQuadric(double c2, double c1, double c0, double& s0, double& s1);
int solveCubic(double c3, double c2, double c1, double c0, double& s0, double& s1, double& s2);
std::vector<double> get_roots(double c3, double c2, double c1, double c0);
std::vector<double> get_roots(double c2, double c1, double c0);

//string convertToString(char* a);

#endif