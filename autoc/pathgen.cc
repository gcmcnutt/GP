#include "pathgen.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a half-sphere
Point3D randomPointInHalfSphere(double radius) {
    Point3D point;
    double theta = ((double) rand() / RAND_MAX) * 2 * M_PI;
    double phi = ((double) rand() / RAND_MAX) * M_PI / 2;
    double r = radius * std::cbrt((double) rand() / RAND_MAX);

    point.x = r * std::sin(phi) * std::cos(theta);
    point.y = r * std::sin(phi) * std::sin(theta);
    point.z = r * std::cos(phi);
    
    return point;
}

// Function to interpolate between points using cubic splines
Point3D cubicInterpolate(const Point3D& p0, const Point3D& p1, const Point3D& p2, const Point3D& p3, double t) {
    Point3D result;
    result.x = 0.5 * ((2 * p1.x) + (-p0.x + p2.x) * t + (2*p0.x - 5*p1.x + 4*p2.x - p3.x) * t*t + (-p0.x + 3*p1.x - 3*p2.x + p3.x) * t*t*t);
    result.y = 0.5 * ((2 * p1.y) + (-p0.y + p2.y) * t + (2*p0.y - 5*p1.y + 4*p2.y - p3.y) * t*t + (-p0.y + 3*p1.y - 3*p2.y + p3.y) * t*t*t);
    result.z = 0.5 * ((2 * p1.z) + (-p0.z + p2.z) * t + (2*p0.z - 5*p1.z + 4*p2.z - p3.z) * t*t + (-p0.z + 3*p1.z - 3*p2.z + p3.z) * t*t*t);
    return result;
}

// Function to generate a smooth random path within a half-sphere
std::vector<Point3D> generateSmoothPath(int numPoints, double radius) {
    std::vector<Point3D> controlPoints;
    std::vector<Point3D> path;

    // Generate random control points
    for (int i = 0; i < numPoints; ++i) {
        controlPoints.push_back(randomPointInHalfSphere(radius));
    }

    // Ensure the path is continuous by looping through control points
    for (size_t i = 1; i < controlPoints.size() - 2; ++i) {
        for (double t = 0; t <= 1; t += 0.1) {
            Point3D interpolatedPoint = cubicInterpolate(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t);
            path.push_back(interpolatedPoint);
        }
    }

    return path;
}

