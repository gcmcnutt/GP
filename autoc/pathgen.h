// Define a structure for 3D points
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>
#include <Eigen/Dense>

#define NUM_SEGMENTS_PER_PATH 16
#define FIELD_SIZE 100.0
#define FIELD_GAP 10.0

class Path {
public:
  Eigen::Vector3d start;
  double distanceFromStart;
  double radiansFromStart;

  void toString(char* output);
};

std::vector<std::vector<Path>> generateSmoothPaths(int numPaths, int numPoints, double radius, double height);

#endif