// Define a structure for 3D points
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>
#include <Eigen/Dense>

class Path {
    public:
        Eigen::Vector3d start;
        double distanceFromStart;

        void toString(char* output);
};

std::vector<Path> generateSmoothPath(int numPoints, double radius);

#endif