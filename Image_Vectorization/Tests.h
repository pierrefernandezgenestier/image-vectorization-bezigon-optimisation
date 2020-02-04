#pragma once

#ifndef TESTS_H_INCLUDED
#define TESTS_H_INCLUDED

#include <iostream>
#include <ostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "Entities.h"
#include "Energy.h"
#include "RMS.h"
#include "Propagation.h"

void test_intersect_bez();
void test_compute_angles(VectorizationData vd);
void test_rms(VectorizationData vd);

void test_compute_energy_angles(VectorizationData vd);
void test_compute_energy_data(VectorizationData vd);
void test_compute_bezier_handles(VectorizationData vd);
void test_espilon(VectorizationData vd);
void test_alpha(VectorizationData vd);


#endif






