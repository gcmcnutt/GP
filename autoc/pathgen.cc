#include "pathgen.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a half-sphere
Point3D randomPointInHalfSphere(double radius) {
    double theta = ((double) rand() / RAND_MAX) * 2 * M_PI;
    double phi = ((double) rand() / RAND_MAX) * M_PI / 2;
    double r = radius * std::cbrt((double) rand() / RAND_MAX);

    double x = r * std::sin(phi) * std::cos(theta);
    double y = r * std::sin(phi) * std::sin(theta);
    double z = r * std::cos(phi);
    z = SIM_INITIAL_ALTITUDE; // TODO constrained Z for now

    return Point3D(x, y, z);
}

// Function to interpolate between points using cubic splines
Point3D cubicInterpolate(const Point3D& p0, const Point3D& p1, const Point3D& p2, const Point3D& p3, double t) {
    double x = 0.5 * ((2 * p1.x) + (-p0.x + p2.x) * t + (2*p0.x - 5*p1.x + 4*p2.x - p3.x) * t*t + (-p0.x + 3*p1.x - 3*p2.x + p3.x) * t*t*t);
    double y = 0.5 * ((2 * p1.y) + (-p0.y + p2.y) * t + (2*p0.y - 5*p1.y + 4*p2.y - p3.y) * t*t + (-p0.y + 3*p1.y - 3*p2.y + p3.y) * t*t*t);
    double z = 0.5 * ((2 * p1.z) + (-p0.z + p2.z) * t + (2*p0.z - 5*p1.z + 4*p2.z - p3.z) * t*t + (-p0.z + 3*p1.z - 3*p2.z + p3.z) * t*t*t);
    return Point3D(x, y, z);
}

// Function to generate a smooth random path within a half-sphere
std::vector<Path> generateSmoothPath(int numPoints, double radius) {
    std::vector<Point3D> controlPoints;
    std::vector<Path> path;

    // Initial control point2
    Point3D initialPoint = {0, 0, SIM_INITIAL_ALTITUDE};
    controlPoints.push_back(initialPoint);
    controlPoints.push_back(initialPoint); // XXX i really want this to be the first point
    Point3D initialPoint2 = {SIM_INITIAL_VELOCITY * 0.5, 0, SIM_INITIAL_ALTITUDE};
    controlPoints.push_back(initialPoint2);

    // Generate random control points
    for (int i = 0; i < numPoints - 2; ++i) {
        controlPoints.push_back(randomPointInHalfSphere(radius));
    }

    // Ensure the path is continuous by looping through control points
    double distance = 0;
    Point3D lastPoint = controlPoints[0];
    for (size_t i = 1; i < controlPoints.size() - 2; ++i) {
        for (double t = 0; t <= 1; t += 0.05) {
            Point3D interpolatedPoint = cubicInterpolate(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t);
            double newDistance = std::sqrt(std::pow(interpolatedPoint.x - lastPoint.x, 2) + std::pow(interpolatedPoint.y - lastPoint.y, 2) + std::pow(interpolatedPoint.z - lastPoint.z, 2));
            distance += newDistance;
            Path pathSegment = {interpolatedPoint, distance};
            path.push_back(pathSegment);
            lastPoint = interpolatedPoint;
        }
    }

    controlPoints.clear();

    return path;
}

void Path::toString(char* output) {
    sprintf(output, "Path: (%f, %f, %f), Distance: %f", start.x, start.y, start.z, distanceFromStart);
}

