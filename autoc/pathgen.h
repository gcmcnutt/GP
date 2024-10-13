// Define a structure for 3D points
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>

#include "minisim.h"

#define NUM_SEGMENTS_PER_PATH 16

std::vector<std::vector<Path>> generateSmoothPaths(int numPaths, int numPoints, double radius, double height);

#endif