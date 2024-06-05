// Define a structure for 3D points
#ifndef POINT3D_H
#define POINT3D_H

#include <vector>

struct Point3D {
    double x, y, z;
};

std::vector<Point3D> generateSmoothPath(int numPoints, double radius);

#endif