#pragma once

#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <iostream>
#include <ostream>
#include <fstream>
#include <csignal>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "Initialisation.h"
#include "Entities.h"
#include "Propagation.h"
#include "Tests.h"
#include "Energy.h"

int main(int argc, char* argv[]);

#endif