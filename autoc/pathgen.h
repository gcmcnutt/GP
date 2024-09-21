// Define a structure for 3D points
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>
#include <Eigen/Dense>

#include "minisim.h"

#define NUM_SEGMENTS_PER_PATH 16
#define FIELD_SIZE 100.0
#define FIELD_GAP 10.0

std::vector<std::vector<Path>> generateSmoothPaths(int numPaths, int numPoints, double radius, double height);

#endif