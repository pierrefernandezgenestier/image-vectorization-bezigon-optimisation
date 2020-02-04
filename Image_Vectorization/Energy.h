#pragma once

#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

#include <iostream>
#include <ostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "RMS.h"
#include "Entities.h"

double energy_spt(Bezigon B, int j);
double energy_spt(Bezigon B);

double energy_apt(Bezigon B, int j);
double energy_apt(Bezigon B);

double energy_hpt(Bezigon B, int j);
double energy_hpt(Bezigon B);

double energy_lpt(Bezigon B, int j);
double energy_lpt(Bezigon B);

double energy_data(VectorizationData vd);

//double energy_partial(VectorizationData vd, int j);
array<double, 2> energy_to_minimize(VectorizationData vd, int j, array<double, 10> vals_inp);
//double energy_to_minimize_prior(VectorizationData vd, int j, array<double, 10> vals_inp);

const double lambda_spt = 1.0;
const double lambda_apt = 0.08;
const double lambda_hpt = 0.8;
const double lambda_lpt = 0.1;
const double lambda_data = 1e-2;

#endif