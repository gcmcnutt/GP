#include "autoc/eval/pathgen.h"
#include "autoc/util/rng.h"
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>


// Function to generate a random point within a cylinder
gp_vec3 randomPointInCylinder(gp_scalar radius, gp_scalar height, gp_scalar base) {
  // Generate random values
  gp_scalar r = radius * std::cbrtf(static_cast<gp_scalar>(rng::randDouble()));
  gp_scalar theta = static_cast<gp_scalar>(rng::randDouble()) * static_cast<gp_scalar>(M_PI * 2.0);
  gp_scalar z = base - static_cast<gp_scalar>(rng::randDouble()) * height;

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
std::vector<std::vector<Path>> generateSmoothPaths(const std::string& method, int numPaths, gp_scalar radius, gp_scalar height, unsigned int baseSeed) {
  std::vector<std::vector<Path>> paths;

  GeneratorMethod* generatorMethod;

  if (method == "random") {
    generatorMethod = new GenerateRandom();
  }
  else if (method == "classic") {
    generatorMethod = new GenerateClassic();
  }
  else if (method == "computedPaths") {
    generatorMethod = new GenerateComputedPaths();
  }
  else if (method == "longSequential") {
    generatorMethod = new GenerateLongSequential();
  }
  else if (method == "aeroStandard") {
    generatorMethod = new GenerateAeroStandard();
  }
  else if (method == "line") {
    generatorMethod = new GenerateLine();
  }
  else if (method == "progressiveDistance") {
    generatorMethod = new GenerateProgressiveDistance();
  }
  else {
    std::cerr << "Unknown path generation method: " << method << std::endl;
    assert(false);
  }

  // For 'random' method: create a NEW seed per path using an mt19937 instance
  // For all other methods (aeroStandard, etc): use the SAME baseSeed for all paths
  bool useNewSeedPerPath = (method == "random");
  std::mt19937 seedGenerator(baseSeed);

  for (int i = 0; i < numPaths; ++i) {
    unsigned int pathSeed = useNewSeedPerPath ? seedGenerator() : baseSeed;

    // Generate path in canonical coordinate frame (0,0,0)
    std::vector<Path> canonicalPath = generatorMethod->method(i, radius, height, 0.0f, pathSeed);

    // Apply NED offset for desktop simulation (paths start at SIM_INITIAL_ALTITUDE)
    gp_vec3 offset(0.0f, 0.0f, SIM_INITIAL_ALTITUDE);  // -25m in NED
    std::vector<Path> offsetPath;
    for (const auto& segment : canonicalPath) {
      Path offsetSegment = segment;
      offsetSegment.start += offset;
      offsetPath.push_back(offsetSegment);
    }

    paths.push_back(offsetPath);
  }

  delete generatorMethod;
  return paths;
}
