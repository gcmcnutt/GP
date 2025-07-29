#include "gp.h"
#include "pathgen.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a cylinder
Eigen::Vector3d randomPointInCylinder(double radius, double height, double base) {
  // Generate random values
  double r = radius * std::cbrt((double)GPrand() / RAND_MAX);
  double theta = ((double)GPrand() / RAND_MAX) * M_PI * 2;
  double z = base - ((double)GPrand() / RAND_MAX) * height;

  // Convert to Cartesian coordinates
  double x = r * std::cos(theta);
  double y = r * std::sin(theta);
  return Eigen::Vector3d(x, y, z);
}

// Function to interpolate between points using cubic splines
Eigen::Vector3d cubicInterpolate(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double t) {
  double x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t * t + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t * t * t);
  double y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t * t + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t * t * t);
  double z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t * t + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t * t * t);
  return Eigen::Vector3d(x, y, z);
}

// Function to generate a smooth random paths within a half-sphere
std::vector<std::vector<Path>> generateSmoothPaths(char* method, int numPaths, double radius, double height) {
  std::vector<std::vector<Path>> paths;

  GeneratorMethod* generatorMethod;

  if (strcmp(method, "random") == 0) {
    generatorMethod = new GenerateRandom();
  }
  else if (strcmp(method, "classic") == 0) {
    generatorMethod = new GenerateClassic();
  }
  else if (strcmp(method, "computedPaths") == 0) {
    generatorMethod = new GenerateComputedPaths();
  }
  else if (strcmp(method, "longSequential") == 0) {
    generatorMethod = new GenerateLongSequential();
  }
  else if (strcmp(method, "line") == 0) {
    generatorMethod = new GenerateLine();
  }
  else {
    std::cerr << "Unknown path generation method: " << method << std::endl;
    assert(false);
  }

  for (size_t i = 0; i < numPaths; ++i) {
    paths.push_back(generatorMethod->method(i, radius, height, SIM_INITIAL_ALTITUDE));
  }

  return paths;
}

