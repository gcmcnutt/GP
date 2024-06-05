// Define a structure for 3D points
#ifndef POINT3D_H
#define POINT3D_H

#include <vector>

class Point3D {
  public:
    double x, y, z;

    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}
};

std::vector<Point3D> generateSmoothPath(int numPoints, double radius);

#endif