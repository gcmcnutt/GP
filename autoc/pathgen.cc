#include "gp.h"
#include "pathgen.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a cylinder
gp_vec3 randomPointInCylinder(gp_scalar radius, gp_scalar height, gp_scalar base) {
  // Generate random values
  gp_scalar r = radius * std::cbrtf(static_cast<gp_scalar>(GPrand()) / static_cast<gp_scalar>(RAND_MAX));
  gp_scalar theta = (static_cast<gp_scalar>(GPrand()) / static_cast<gp_scalar>(RAND_MAX)) * static_cast<gp_scalar>(M_PI * 2.0);
  gp_scalar z = base - (static_cast<gp_scalar>(GPrand()) / static_cast<gp_scalar>(RAND_MAX)) * height;

  // Convert to Cartesian coordinates
  gp_scalar x = r * std::cos(theta);
  gp_scalar y = r * std::sin(theta);
  return gp_vec3(x, y, z);
}

// Function to interpolate between points using cubic splines
gp_vec3 cubicInterpolate(const gp_vec3& p0, const gp_vec3& p1, const gp_vec3& p2, const gp_vec3& p3, gp_scalar t) {
  const gp_scalar half = static_cast<gp_scalar>(0.5f);
  gp_scalar x = half * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t * t + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t * t * t);
  gp_scalar y = half * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t * t + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t * t * t);
  gp_scalar z = half * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t * t + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t * t * t);
  return gp_vec3(x, y, z);
}

// Function to generate a smooth random paths within a half-sphere
std::vector<std::vector<Path>> generateSmoothPaths(char* method, int numPaths, gp_scalar radius, gp_scalar height, unsigned int baseSeed) {
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
  else if (strcmp(method, "aeroStandard") == 0) {
    generatorMethod = new GenerateAeroStandard();
  }
  else if (strcmp(method, "line") == 0) {
    generatorMethod = new GenerateLine();
  }
  else {
    std::cerr << "Unknown path generation method: " << method << std::endl;
    assert(false);
  }

  // For 'random' method: create a NEW seed per path using an mt19937 instance
  // For all other methods (aeroStandard, etc): use the SAME baseSeed for all paths
  bool useNewSeedPerPath = (strcmp(method, "random") == 0);
  std::mt19937 seedGenerator(baseSeed);

  for (int i = 0; i < numPaths; ++i) {
    unsigned int pathSeed = useNewSeedPerPath ? seedGenerator() : baseSeed;
    paths.push_back(generatorMethod->method(i, radius, height, SIM_INITIAL_ALTITUDE, pathSeed));
  }

  delete generatorMethod;
  return paths;
}
