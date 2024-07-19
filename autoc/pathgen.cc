#include "gp.h"
#include "pathgen.h"
#include "minisim.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a cylinder
Eigen::Vector3d randomPointInCylinder(double radius, double height, double base = SIM_INITIAL_ALTITUDE) {
    // Generate random values
    double r = radius * std::cbrt((double) GPrand() / RAND_MAX);
    double theta = ((double) GPrand() / RAND_MAX) * M_PI * 2;
    double z = base - ((double) GPrand() / RAND_MAX) * height;

    // Convert to Cartesian coordinates
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    return Eigen::Vector3d(x, y, z);
}

// Function to interpolate between points using cubic splines
Eigen::Vector3d cubicInterpolate(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double t) {
    double x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2*p0[0] - 5*p1[0] + 4*p2[0] - p3[0]) * t*t + (-p0[0] + 3*p1[0] - 3*p2[0] + p3[0]) * t*t*t);
    double y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2*p0[1] - 5*p1[1] + 4*p2[1] - p3[1]) * t*t + (-p0[1] + 3*p1[1] - 3*p2[1] + p3[1]) * t*t*t);
    double z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2*p0[2] - 5*p1[2] + 4*p2[2] - p3[2]) * t*t + (-p0[2] + 3*p1[2] - 3*p2[2] + p3[2]) * t*t*t);
    return Eigen::Vector3d(x, y, z);
}

std::vector<Path> generateSmoothPath(int numPoints, double radius, double height, double base) {
    std::vector<Eigen::Vector3d> controlPoints;
    std::vector<Path> path;

// #define PATHGEN_FIXED_PATH 1
#ifdef PATHGEN_FIXED_PATH
    controlPoints.push_back({0, 0, base});

    // Sin
    double x, y, z;
    for (size_t i = 0; i < numPoints; ++i) {
        x = -(cos(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2 - SIM_PATH_BOUNDS/2);
        y = sin(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        controlPoints.push_back(Eigen::Vector3d(x, y, z)); 
        z = base - i;
    }
    for (size_t i = 0; i < numPoints; ++i) {
        x = sin(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        z = base - SIM_PATH_BOUNDS/2 + cos(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS/2;
        y = i;
        controlPoints.push_back(Eigen::Vector3d(x, y, z)); 
    }

#else
    // Generate random control points
    for (size_t i = 0; i < numPoints; ++i) {
        controlPoints.push_back(randomPointInCylinder(radius, height));
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
                first = false;
            }

            lastPoint = interpolatedPoint;
            lastDirection = newDirection;

            // stop around what should be TOTAL_TIME seconds of data
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
std::vector<std::vector<Path>> generateSmoothPaths(int numPaths, int numPoints, double radius, double height) {
    std::vector<std::vector<Path>> paths;

    // std::vector<Path> path;
    // double distance = 0.0;
    // for (double i = 0; i <= SIM_PATH_RADIUS_LIMIT; i += 5) {
    //     Eigen::Vector3d interpolatedPoint = {i, i, SIM_INITIAL_ALTITUDE};
    //     Path pathSegment = {interpolatedPoint, distance, 0.0};
    //     path.push_back(pathSegment);
    //     distance += 5;
    // }
    // paths.push_back(path);

    for (size_t i = 0; i < numPaths; ++i) {
        paths.push_back(generateSmoothPath(numPoints, radius, height, SIM_INITIAL_ALTITUDE));
    }

    return paths;
}

void Path::toString(char* output) {
    sprintf(output, "Path: (%f, %f, %f), Odometer: %f, Turnmeter: %f", start[0], start[1], start[2], distanceFromStart, radiansFromStart);
}

