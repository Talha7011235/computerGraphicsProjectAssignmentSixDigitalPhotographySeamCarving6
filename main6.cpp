// Computer Graphic Project Assignment 6 Digital Photography Seam Carving
// main6.cpp
#define cimg_display 0
#include "CImg.h"
#include "carver.h"
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace cimg_library;

typedef unsigned int uint32_t;

int main(int argc, char *argv[])
{
    Carver carver;
    carver.load(argv[1]);

    carver.resize(atoi(argv[3]), atoi(argv[4]));

    carver.save(argv[2]);
    return 0;
}

