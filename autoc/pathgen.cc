#include "gp.h"
#include "pathgen.h"
#include "minisim.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a half-sphere
Eigen::Vector3d randomPointInHalfSphere(double radius) {
    double theta = ((double) GPrand() / RAND_MAX) * 2 * M_PI;
    double phi = ((double) GPrand() / RAND_MAX) * M_PI / 2;
    double r = radius * std::cbrt((double) GPrand() / RAND_MAX);

    double x = r * std::sin(phi) * std::cos(theta);
    double y = r * std::sin(phi) * std::sin(theta);
    double z = (SIM_MIN_ELEVATION - r) * std::cos(phi);

    return Eigen::Vector3d(x, y, z);
}

// Function to interpolate between points using cubic splines
Eigen::Vector3d cubicInterpolate(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double t) {
    double x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2*p0[0] - 5*p1[0] + 4*p2[0] - p3[0]) * t*t + (-p0[0] + 3*p1[0] - 3*p2[0] + p3[0]) * t*t*t);
    double y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2*p0[1] - 5*p1[1] + 4*p2[1] - p3[1]) * t*t + (-p0[1] + 3*p1[1] - 3*p2[1] + p3[1]) * t*t*t);
    double z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2*p0[2] - 5*p1[2] + 4*p2[2] - p3[2]) * t*t + (-p0[2] + 3*p1[2] - 3*p2[2] + p3[2]) * t*t*t);
    return Eigen::Vector3d(x, y, z);
}

std::vector<Path> generateSmoothPath(int numPoints, double radius) {
    std::vector<Eigen::Vector3d> controlPoints;
    std::vector<Path> path;

    // Initial control points in forward direction
    Eigen::Vector3d initialPoint = {0, 0, SIM_INITIAL_ALTITUDE};
    controlPoints.push_back(initialPoint);

#define PATHGEN_FIXED_PATH 1
#ifdef PATHGEN_FIXED_PATH
    // Sin
    double x, y, z;
    for (size_t i = 0; i < numPoints; ++i) {
        x = sin(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        z = SIM_INITIAL_ALTITUDE - SIM_PATH_BOUNDS/2 + cos(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        y = i;
        controlPoints.push_back(Eigen::Vector3d(x, y, z)); 
    }
    for (size_t i = 0; i < numPoints; ++i) {
        x = sin(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        y = cos(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        controlPoints.push_back(Eigen::Vector3d(x, y, z)); 
        z = SIM_INITIAL_ALTITUDE - i;
    }
#else
    // Generate random control points
    for (size_t i = 1; i <= numPoints; ++i) {
        controlPoints.push_back(randomPointInHalfSphere(radius));
    }
#endif

    // Ensure the path is continuous by looping through control points
    double odometer = 0;
    double turnmeter = 0;
    Eigen::Vector3d lastPoint;
    Eigen::Vector3d lastDirection;
    bool first = true;

    for (size_t i = 1; i < controlPoints.size() - 3; ++i) {
        for (double t = 0; t <= 1; t += 0.05) {
            Eigen::Vector3d interpolatedPoint = cubicInterpolate(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t);
            Eigen::Vector3d newDirection;

            if (!first) {
                // Compute the new distance
                double newDistance = (interpolatedPoint - lastPoint).norm();

                // Compute new new direction
                newDirection = (interpolatedPoint - lastPoint).normalized();

                // difference to prior direction
                double dVector = lastDirection.dot(newDirection);
                double dAngle = std::acos(std::clamp(dVector / (lastDirection.norm() * newDirection.norm()), -1.0, 1.0));

                // add next segment            
                Path pathSegment = {interpolatedPoint, odometer, turnmeter};
                path.push_back(pathSegment);

                // char output[200];
                // pathSegment.toString(output);
                // cout << output << endl;
                
                odometer += newDistance;
                turnmeter += dAngle;
            } else {
                newDirection = {1.0, 0.0, 0.0};
                first = false;
            }

            lastPoint = interpolatedPoint;
            lastDirection = newDirection;

            // stop around what should be 75 seconds of data
            if (odometer > SIM_TOTAL_TIME * SIM_INITIAL_VELOCITY) {
                goto exitLoop;
            }
        }
    }
exitLoop:

    controlPoints.clear();

    return path;
}

// Function to generate a smooth random paths within a half-sphere
std::vector<std::vector<Path>> generateSmoothPaths(int numPaths, int numPoints, double radius) {
    std::vector<std::vector<Path>> paths;

    for (size_t i = 0; i < numPaths; ++i) {
        paths.push_back(generateSmoothPath(numPoints, radius));
    }

    return paths;
}

void Path::toString(char* output) {
    sprintf(output, "Path: (%f, %f, %f), Odometer: %f, Turnmeter: %f", start[0], start[1], start[2], distanceFromStart, radiansFromStart);
}

