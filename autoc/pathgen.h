// Path generation for GP training
//
// USAGE PATTERN FOR autoc.cc:
// TODO: Update autoc.cc to pass RandomPathSeedB from config manager as baseSeed parameter
//       Example usage:
//         unsigned int seed = extraConfig.randomPathSeedB;  // Read from .ini (default: 67890)
//         auto paths = generateSmoothPaths("aeroStandard", 6, radius, height, seed);
//
//       For 'aeroStandard': All 6 paths use the SAME seed (paths are deterministic based on pathIndex)
//       For 'random': Each path gets a DIFFERENT seed (generated from baseSeed using mt19937)
//
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>
#include <random>

#include "minisim.h"
#include "gp_types.h"

#define NUM_SEGMENTS_PER_PATH 16

gp_vec3 randomPointInCylinder(gp_scalar radius, gp_scalar height, gp_scalar base = SIM_INITIAL_ALTITUDE);
gp_vec3 cubicInterpolate(const gp_vec3& p0, const gp_vec3& p1, const gp_vec3& p2, const gp_vec3& p3, gp_scalar t);
// Generate paths with controllable random seed
// For 'aeroStandard' and deterministic methods: use ONE seed for all numPaths
// For 'random' method: generate a NEW seed for each path
std::vector<std::vector<Path>> generateSmoothPaths(char* method, int numPaths, gp_scalar radius, gp_scalar height, unsigned int baseSeed);

class GeneratorMethod {
public:
  virtual ~GeneratorMethod() = default;
  // pathIndex: which path variant to generate (0-5 for aeroStandard, etc)
  // seed: random seed for reproducible generation (creates private PRNG instance)
  // NOTE: Callers should use ONE seed for all paths in aeroStandard/computedPaths,
  //       but generate NEW seeds per-path for 'random' method
  virtual std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) = 0;

protected:
  // Generate random point in cylinder using local mt19937 PRNG (does not affect global rand state)
  static gp_vec3 localRandomPointInCylinder(std::mt19937& rng, gp_scalar radius, gp_scalar height, gp_scalar base = SIM_INITIAL_ALTITUDE) {
    std::uniform_real_distribution<gp_scalar> dist(static_cast<gp_scalar>(0.0), static_cast<gp_scalar>(1.0));
    gp_scalar r = radius * std::cbrt(dist(rng));
    gp_scalar theta = dist(rng) * static_cast<gp_scalar>(M_PI * 2.0);
    gp_scalar z = base - dist(rng) * height;
    gp_scalar x = r * std::cos(theta);
    gp_scalar y = r * std::sin(theta);
    return gp_vec3(x, y, z);
  }
};

class GenerateRandom : public GeneratorMethod {
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) override {
    std::vector<gp_vec3> controlPoints;
    std::vector<Path> path;
    int numPoints = NUM_SEGMENTS_PER_PATH;

    // Initialize local mt19937 PRNG (does not affect global rand())
    std::mt19937 rng(seed);

    // Generate random control points
    for (size_t i = 0; i < numPoints; ++i) {
      controlPoints.push_back(localRandomPointInCylinder(rng, radius, height));
    }

    // Ensure the path is continuous by looping through control points
    gp_scalar odometer = 0;
    gp_scalar turnmeter = 0;
    gp_vec3 lastPoint;
    gp_vec3 lastDirection;
    bool first = true;

    for (size_t i = 1; i < controlPoints.size() - 3; ++i) {
      for (gp_scalar t = 0; t <= 1; t += static_cast<gp_scalar>(0.05f)) {
        gp_vec3 interpolatedPoint = cubicInterpolate(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t);
        gp_vec3 newDirection;

        if (!first) {
          // Compute the new distance
          gp_scalar newDistance = (interpolatedPoint - lastPoint).norm();

          // Compute new new direction
          newDirection = (interpolatedPoint - lastPoint).normalized();

          // difference to prior direction
          gp_scalar dVector = lastDirection.dot(newDirection);
          gp_scalar dAngle = std::acos(std::clamp(dVector / (lastDirection.norm() * newDirection.norm()), -1.0f, 1.0f));

          // add next segment with simulation timestamp
          gp_scalar simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), odometer, turnmeter, simTimeMsec);
          path.push_back(pathSegment);

          odometer += newDistance;
          turnmeter += dAngle;
        }
        else {
          first = false;
        }

        lastPoint = interpolatedPoint;
        lastDirection = newDirection;

        // stop around what should be TOTAL_TIME seconds of data
        if (odometer > SIM_TOTAL_TIME_MSEC / 1000.0 * SIM_RABBIT_VELOCITY) {
          goto exitLoop;
        }
      }
    }
  exitLoop:

    return path;
  }
};

class GenerateClassic : public GeneratorMethod {
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) override {
    // Classic path is deterministic, seed is unused
    std::vector<gp_vec3> controlPoints;
    std::vector<Path> path;
    int numPoints = NUM_SEGMENTS_PER_PATH;

    controlPoints.push_back({ 0, 0, base });

    // Sin
    gp_scalar x, y, z;
    for (size_t i = 0; i < numPoints; ++i) {
      x = -(std::cos(static_cast<gp_scalar>(2 * M_PI) * static_cast<gp_scalar>(i) / static_cast<gp_scalar>(numPoints)) * SIM_PATH_BOUNDS / static_cast<gp_scalar>(2.0f) - SIM_PATH_BOUNDS / static_cast<gp_scalar>(2.0f));
      y = std::sin(static_cast<gp_scalar>(2 * M_PI) * static_cast<gp_scalar>(i) / static_cast<gp_scalar>(numPoints)) * SIM_PATH_BOUNDS / static_cast<gp_scalar>(2.0f);
      z = base - static_cast<gp_scalar>(i);
      controlPoints.push_back(gp_vec3(x, y, z));
    }

    gp_scalar lastX = x;
    gp_scalar lastY = y;
    gp_scalar lastZ = z;
    for (size_t i = 0; i < numPoints; ++i) {
      y = lastX + std::sin(static_cast<gp_scalar>(2 * M_PI) * static_cast<gp_scalar>(i) / static_cast<gp_scalar>(numPoints)) * SIM_PATH_BOUNDS / static_cast<gp_scalar>(3.0f);
      z = lastZ - SIM_PATH_BOUNDS / static_cast<gp_scalar>(2.0f) + std::cos(static_cast<gp_scalar>(2 * M_PI) * static_cast<gp_scalar>(i) / static_cast<gp_scalar>(numPoints)) * SIM_PATH_BOUNDS / static_cast<gp_scalar>(2.0f);
      x = lastY - static_cast<gp_scalar>(i);
      controlPoints.push_back(gp_vec3(x, y, z));
    }

    // Ensure the path is continuous by looping through control points
    gp_scalar odometer = 0;
    gp_scalar turnmeter = 0;
    gp_vec3 lastPoint;
    gp_vec3 lastDirection;
    bool first = true;

    for (size_t i = 1; i < controlPoints.size() - 3; ++i) {
      for (gp_scalar t = 0; t <= 1; t += static_cast<gp_scalar>(0.05f)) {
        gp_vec3 interpolatedPoint = cubicInterpolate(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t);
        gp_vec3 newDirection;

        if (!first) {
          // Compute the new distance
          gp_scalar newDistance = (interpolatedPoint - lastPoint).norm();

          // Compute new new direction
          newDirection = (interpolatedPoint - lastPoint).normalized();

          // difference to prior direction
          gp_scalar dVector = lastDirection.dot(newDirection);
          gp_scalar dAngle = std::acos(std::clamp(dVector / (lastDirection.norm() * newDirection.norm()), -1.0f, 1.0f));

          // add next segment with simulation timestamp
          gp_scalar simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), odometer, turnmeter, simTimeMsec);
          path.push_back(pathSegment);

          odometer += newDistance;
          turnmeter += dAngle;
        }
        else {
          first = false;
        }

        lastPoint = interpolatedPoint;
        lastDirection = newDirection;

        // stop around what should be TOTAL_TIME seconds of data
        if (odometer > SIM_TOTAL_TIME_MSEC / 1000.0 * SIM_RABBIT_VELOCITY) {
          goto exitLoop;
        }
      }
    }
  exitLoop:

    return path;
  }
};

enum PathType {
  HorizontalCircle = 0,
  VerticalCircle,
  FigureEight,
  HorizontalSquare,
  VerticalSquare,
  FortyFiveLoop,
  END_MARKER
};

class GenerateComputedPaths : public GeneratorMethod {
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) override {
    // Computed paths are deterministic, seed is unused
    std::vector<Path> path;

    PathType pathType = static_cast<PathType>(pathIndex % static_cast<int>(PathType::END_MARKER));

    switch(pathType) {
      case PathType::HorizontalCircle: {
        gp_scalar cSize = static_cast<gp_scalar>(0.8f) * radius;
        for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
          gp_vec3 interpolatedPoint = { cSize * std::cos(turn), cSize * std::sin(turn), base - SIM_INITIAL_LOCATION_DITHER };
          gp_scalar simTimeMsecLocal = (turn * radius / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), turn * radius, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::VerticalCircle: {
        gp_scalar cSize = static_cast<gp_scalar>(0.8f) * radius;
        for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
          gp_vec3 interpolatedPoint = { cSize * std::cos(turn), 0.0f, base - SIM_INITIAL_LOCATION_DITHER - cSize * std::sin(turn)};
          gp_scalar simTimeMsecLocal = (turn * cSize / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), turn * cSize, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::FigureEight: {
        gp_scalar cSize = static_cast<gp_scalar>(0.4f) * radius;
        for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
          gp_vec3 interpolatedPoint = { (cSize * std::cos(turn)) - cSize, 0.0f, base - SIM_INITIAL_LOCATION_DITHER - (cSize * std::sin(turn))};
          gp_scalar simTimeMsecLocal = ((turn * cSize) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), (turn * cSize), 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }
        for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
          gp_vec3 interpolatedPoint = { (-cSize * std::cos(turn)) + cSize, 0.0f, base - SIM_INITIAL_LOCATION_DITHER - (cSize * std::sin(turn))};
          gp_scalar simTimeMsecLocal = (((turn * cSize) + (static_cast<gp_scalar>(M_PI * 2.0) * cSize)) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), (turn * cSize) + (static_cast<gp_scalar>(M_PI * 2.0) * cSize), 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::HorizontalSquare: {
        gp_scalar sSize = static_cast<gp_scalar>(0.7f) * radius;
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { span, -sSize, base - SIM_INITIAL_LOCATION_DITHER };
          gp_scalar simTimeMsecLocal = (span / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }  
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { sSize, span, base - SIM_INITIAL_LOCATION_DITHER };
          gp_scalar simTimeMsecLocal = ((static_cast<gp_scalar>(2.0f) * sSize + span) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), static_cast<gp_scalar>(2.0f) * sSize + span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        } 
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { -span, sSize, base - SIM_INITIAL_LOCATION_DITHER };
          gp_scalar simTimeMsecLocal = ((static_cast<gp_scalar>(4.0f) * sSize + span) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), static_cast<gp_scalar>(4.0f) * sSize + span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }  
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { -sSize, -span, base - SIM_INITIAL_LOCATION_DITHER };
          gp_scalar simTimeMsecLocal = ((static_cast<gp_scalar>(6.0f) * sSize + span) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), static_cast<gp_scalar>(6.0f) * sSize + span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }      
        break;
      }
      case PathType::VerticalSquare: {
        gp_scalar sSize = static_cast<gp_scalar>(0.7f) * radius;
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { span, 0.0f, base - SIM_INITIAL_LOCATION_DITHER - sSize };
          gp_scalar simTimeMsecLocal = (span / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }  
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { sSize, 0.0f, base - SIM_INITIAL_LOCATION_DITHER + span };
          gp_scalar simTimeMsecLocal = ((static_cast<gp_scalar>(2.0f) * sSize + span) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), static_cast<gp_scalar>(2.0f) * sSize + span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        } 
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { -span, 0.0f, base - SIM_INITIAL_LOCATION_DITHER + sSize };
          gp_scalar simTimeMsecLocal = ((static_cast<gp_scalar>(4.0f) * sSize + span) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), static_cast<gp_scalar>(4.0f) * sSize + span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }  
        for (gp_scalar span = -sSize; span < sSize; span += static_cast<gp_scalar>(0.1f)) {
          gp_vec3 interpolatedPoint = { -sSize, 0.0f, base - SIM_INITIAL_LOCATION_DITHER - span };
          gp_scalar simTimeMsecLocal = ((static_cast<gp_scalar>(6.0f) * sSize + span) / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), static_cast<gp_scalar>(6.0f) * sSize + span, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }      
        break;
      }
      case PathType::FortyFiveLoop: {
        gp_scalar cSize = static_cast<gp_scalar>(0.8f) * radius;
        for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
          gp_vec3 interpolatedPoint = {
            cSize * std::cos(turn),
            cSize * std::sin(turn),
            base - SIM_INITIAL_LOCATION_DITHER - cSize * std::sin(turn)
          };
          gp_scalar simTimeMsecLocal = (turn * cSize / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), turn * cSize, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }
        break;
      }
      default: {
        std::cerr << "Unknown path type: " << pathType << std::endl;
        assert(false);
      }
    }

    return path;
  }
};

class GenerateLine : public GeneratorMethod {
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) override {
    // Line path is deterministic, seed is unused

    std::vector<Path> path;
    gp_scalar distance = 0.0f;
    for (gp_scalar i = -SIM_PATH_RADIUS_LIMIT; i <= SIM_PATH_RADIUS_LIMIT; i += static_cast<gp_scalar>(5.0f)) {
      gp_vec3 interpolatedPoint = { static_cast<gp_scalar>(-30.0f), i, SIM_INITIAL_ALTITUDE };
      gp_scalar simTimeMsecLocal = (distance / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
      Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), distance, 0.0f, simTimeMsecLocal);
      path.push_back(pathSegment);
      distance += static_cast<gp_scalar>(5.0f);
    }
    return path;
  }
};

enum AeroStandardPathType {
  StraightAndLevel = 0,
  SpiralClimb,
  HorizontalFigureEight,
  FortyFiveDegreeAngledLoop,
  HighPerchSplitS,
  SeededRandomB,
  AERO_END_MARKER  // Total: 6 paths
};

class GenerateAeroStandard : public GeneratorMethod {
public:
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) override {
    std::vector<Path> path;
    gp_scalar totalDistance = 0.0f;

    AeroStandardPathType pathType = static_cast<AeroStandardPathType>(pathIndex % static_cast<int>(AeroStandardPathType::AERO_END_MARKER));

    // Entry point: heading south at x=0, y=0, z=-25
    gp_vec3 entryPoint(0.0f, 0.0f, base);

    switch(pathType) {
      case AeroStandardPathType::StraightAndLevel: {
        // Straight and level: racetrack pattern
        // 1. Head south for 20 meters
        addStraightSegment(path, entryPoint, gp_vec3(-1.0f, 0.0f, 0.0f), 20.0f, totalDistance);

        // 2. Right 180 turn with 20m radius
        gp_vec3 turn1Start = path.back().start;
        addHorizontalTurn(path, turn1Start, 20.0f, M_PI, true, totalDistance); // true = right turn (clockwise)

        // 3. Head north for 40 meters
        gp_vec3 northStart = path.back().start;
        addStraightSegment(path, northStart, gp_vec3(1.0f, 0.0f, 0.0f), 40.0f, totalDistance);

        // 4. Right 180 turn with 20m radius
        gp_vec3 turn2Start = path.back().start;
        addHorizontalTurn(path, turn2Start, 20.0f, M_PI, true, totalDistance);

        // 5. Head south back to near origin for 20 meters
        gp_vec3 returnStart = path.back().start;
        addStraightSegment(path, returnStart, gp_vec3(-1.0f, 0.0f, 0.0f), 20.0f, totalDistance);
        break;
      }

      case AeroStandardPathType::SpiralClimb: {
        // 1. Head south for 20 meters at z=-25
        addStraightSegment(path, entryPoint, gp_vec3(-1.0f, 0.0f, 0.0f), 20.0f, totalDistance);

        // 2. LEFT 540 degree turn (1.5 circles) clockwise when viewed from above, with 20m radius, climbing from z=-25 to z=-75
        // NOTE: "left" in spec means the turn goes counterclockwise in standard math (but appears clockwise when looking down +z)
        gp_vec3 spiralStart = path.back().start;
        gp_scalar climbAmount = -75.0f - spiralStart[2]; // should be -50
        addSpiralTurn(path, spiralStart, 20.0f, 540.0f * M_PI / 180.0f, true, climbAmount, totalDistance); // true = clockwise from above

        // 3. Head north for 40 meters at z=-75
        gp_vec3 northStart = path.back().start;
        addStraightSegment(path, northStart, gp_vec3(1.0f, 0.0f, 0.0f), 40.0f, totalDistance);
        break;
      }

      case AeroStandardPathType::HorizontalFigureEight: {
        // Same as longSequential (reuse the figure-8 logic)
        gp_scalar loopRadius = 20.0f;

        // Lead-in: 1s straight south
        const gp_scalar leadSeconds = 1.0f;
        const gp_scalar leadDistance = SIM_RABBIT_VELOCITY * leadSeconds;
        addStraightSegment(path, entryPoint, gp_vec3(-1.0f, 0.0f, 0.0f), leadDistance, totalDistance);

        gp_vec3 loopOrigin = path.back().start;

        // Left horizontal loop
        addHorizontalLoop(path, loopOrigin, loopRadius, false, totalDistance);

        // Right horizontal loop
        addHorizontalLoop(path, loopOrigin, loopRadius, true, totalDistance);
        break;
      }

      case AeroStandardPathType::FortyFiveDegreeAngledLoop: {
        // 45 degree angled loop - vertical loop rotated 45° around x-axis
        // Loop goes UP (into negative Z), not down
        // Entry at bottom of loop at origin (0, 0, base)
        const gp_scalar cos45 = std::sqrt(2.0f) / 2.0f;
        const gp_scalar sin45 = std::sqrt(2.0f) / 2.0f;

        gp_scalar loopRadius = 15.0f;
        // Rotate by -45° to make loop go up: (x, z*sin45, centerAlt - z*cos45)
        // At angle=-π/2, z_plane=-r, want point at (0, 0, base)
        // y = z_plane*sin45 + yOffset = -r*sin45 + yOffset = 0 → yOffset = r*sin45
        // z = centerAlt - z_plane*cos45 = centerAlt + r*cos45 = base → centerAlt = base - r*cos45
        gp_scalar centerAlt = base - loopRadius * cos45;
        gp_scalar yOffset = loopRadius * sin45;

        for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += 0.5f) {  // Reduced from 0.05 to 0.5 rad
          // Start at bottom of loop: turn starts at -90° (or 3π/2)
          // So at turn=0, we're at the bottom (south side)
          gp_scalar angle = turn - static_cast<gp_scalar>(M_PI / 2.0); // Start at bottom

          // Vertical loop in xz plane before rotation
          // NEGATE x_plane so loop continues southward (negative x) from entry
          gp_scalar x_plane = -loopRadius * std::cos(angle);
          gp_scalar z_plane = loopRadius * std::sin(angle);

          // Rotate around x-axis by -45° (to go up, not down) and translate
          gp_vec3 interpolatedPoint = {
            x_plane,
            z_plane * sin45 + yOffset,
            centerAlt - z_plane * cos45
          };

          if (!path.empty()) {
            gp_scalar distance = (interpolatedPoint - path.back().start).norm();
            totalDistance += distance;
          }
          gp_scalar simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0f;
          Path pathSegment = Path(interpolatedPoint, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsecLocal);
          path.push_back(pathSegment);
        }
        break;
      }

      case AeroStandardPathType::HighPerchSplitS: {
        // Path 4: High Perch with Split-S Reversal
        // Complex multi-segment maneuver combining climb, perch, and reversal

        // Segment 1: Head south for 20m
        addStraightSegment(path, entryPoint, gp_vec3(-1.0f, 0.0f, 0.0f), 20.0f, totalDistance);

        // Segment 2: 180° climbing left turn (clockwise from above) with 20m radius, climbing 20m to z=-45
        gp_vec3 seg2Start = path.back().start;
        gp_scalar climbAmount = -20.0f;  // Climb from z=-25 to z=-45
        addSpiralTurn(path, seg2Start, 20.0f, M_PI, true, climbAmount, totalDistance); // true = clockwise (left turn from pilot)

        // Segment 3: Diagonal climb north from z=-45 to z=-65 over ~40m horizontal distance
        gp_vec3 seg3Start = path.back().start;
        gp_vec3 headingNorth(1.0f, 0.0f, 0.0f);
        gp_scalar horizontalDist = 40.0f;
        gp_scalar verticalClimb = -20.0f;  // From z=-45 to z=-65
        gp_vec3 climbVector = headingNorth * horizontalDist + gp_vec3(0.0f, 0.0f, verticalClimb);
        gp_scalar climbDistance = climbVector.norm();
        addStraightSegment(path, seg3Start, climbVector.normalized(), climbDistance, totalDistance);

        // Segment 4: Hard 150° right turn with 5m radius at high perch
        gp_vec3 seg4Start = path.back().start;
        addHorizontalTurn(path, seg4Start, 5.0f, 150.0f * M_PI / 180.0f, true, totalDistance); // true = right turn

        // Calculate heading after 150° right turn from north
        // After 150° clockwise turn from north (0°): heading = 150° in NED
        gp_scalar headingAngle = 150.0f * M_PI / 180.0f;
        gp_vec3 headingSW(std::cos(headingAngle), std::sin(headingAngle), 0.0f);

        // Segment 5: Fly 30m in new heading (150° from north)
        gp_vec3 seg5Start = path.back().start;
        addStraightSegment(path, seg5Start, headingSW, 30.0f, totalDistance);

        // Segment 6: Continue in same heading until x=-20
        gp_vec3 seg6Start = path.back().start;
        if (std::abs(headingSW[0]) > 0.001f) {
          gp_scalar distToX20 = (-20.0f - seg6Start[0]) / headingSW[0];
          if (distToX20 > 0) {
            addStraightSegment(path, seg6Start, headingSW, distToX20, totalDistance);
          }
        }

        // Segment 7: 180° pitch-down loop (Split-S maneuver) at x=-20
        gp_vec3 seg7Start = path.back().start;
        addPitchDownLoop(path, seg7Start, headingSW, 15.0f, totalDistance);

        // Segment 8: Head north to x=40 (heading is now reversed after Split-S)
        gp_vec3 seg8Start = path.back().start;
        gp_vec3 headingNE = -headingSW;  // Reversed heading after loop
        if (std::abs(headingNE[0]) > 0.001f) {
          gp_scalar distToX40 = (40.0f - seg8Start[0]) / headingNE[0];
          addStraightSegment(path, seg8Start, headingNE.normalized(), std::abs(distToX40), totalDistance);
        }
        break;
      }

      case AeroStandardPathType::SeededRandomB: {
        // Generate seeded random path using local mt19937 PRNG (does not affect global rand state)
        std::mt19937 rng(seed);

        std::vector<gp_vec3> controlPoints;
        int numPoints = NUM_SEGMENTS_PER_PATH;

        // Generate random control points within cylinder
        for (int i = 0; i < numPoints; ++i) {
          controlPoints.push_back(localRandomPointInCylinder(rng, radius, height, base));
        }

        // Generate smooth path through control points
        gp_scalar odometer = 0;
        gp_scalar turnmeter = 0;
        gp_vec3 lastPoint;
        gp_vec3 lastDirection;
        bool first = true;

        for (size_t i = 1; i < controlPoints.size() - 3; ++i) {
          for (gp_scalar t = 0; t <= 1; t += static_cast<gp_scalar>(0.05f)) {
            gp_vec3 interpolatedPoint = cubicInterpolate(controlPoints[i - 1], controlPoints[i],
                                                         controlPoints[i + 1], controlPoints[i + 2], t);
            if (!first) {
              gp_scalar newDistance = (interpolatedPoint - lastPoint).norm();
              gp_vec3 newDirection = (interpolatedPoint - lastPoint).normalized();
              gp_scalar dVector = lastDirection.dot(newDirection);
              gp_scalar dAngle = std::acos(std::clamp(dVector / (lastDirection.norm() * newDirection.norm()), -1.0f, 1.0f));

              gp_scalar simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
              path.push_back(Path(interpolatedPoint, gp_vec3::UnitX(), odometer, turnmeter, simTimeMsec));

              odometer += newDistance;
              turnmeter += dAngle;
              lastDirection = newDirection;
            } else {
              first = false;
              lastDirection = (interpolatedPoint - entryPoint).normalized();
            }
            lastPoint = interpolatedPoint;
          }
        }
        break;
      }

      default: {
        std::cerr << "Unknown aero standard path type: " << pathType << std::endl;
        assert(false);
      }
    }

    return path;
  }

private:
  void addStraightSegment(std::vector<Path>& path, const gp_vec3& start, const gp_vec3& direction,
                          gp_scalar distance, gp_scalar& totalDistance) {
    const gp_scalar step = 1.0f;  // 1m spacing (baseline from horizontal 8 analysis)
    gp_vec3 dir = direction.normalized();

    // Start from step if path is not empty to avoid duplicating the last point
    gp_scalar startD = path.empty() ? 0.0f : step;

    for (gp_scalar d = startD; d <= distance; d += step) {
      gp_vec3 point = start + dir * d;
      if (!path.empty()) {
        gp_scalar segmentDist = (point - path.back().start).norm();
        totalDistance += segmentDist;
      }
      gp_scalar simTimeMsec = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0f;
      path.push_back(Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsec));
    }
  }

  void addHorizontalTurn(std::vector<Path>& path, const gp_vec3& start, gp_scalar radius,
                         gp_scalar angleRadians, bool clockwise, gp_scalar& totalDistance) {
    const gp_scalar step = 0.05f;  // 0.05 rad (~3°) spacing - baseline from horizontal 8

    // Determine initial heading from last two points if possible
    gp_vec3 heading(-1.0f, 0.0f, 0.0f); // default south
    if (path.size() >= 2) {
      heading = (path.back().start - path[path.size()-2].start).normalized();
    }

    // Project heading onto XY plane (ignore Z component for horizontal turn)
    gp_vec3 headingXY(heading[0], heading[1], 0.0f);
    headingXY = headingXY.normalized();

    // Right perpendicular in xy plane: rotate heading 90° RIGHT (clockwise when viewed from above)
    // For heading=(-1,0,0) (south), right perpendicular is (0,-1,0) (west)
    // Formula: rotate (x,y) by -90° (clockwise) = (y, -x)
    gp_vec3 rightPerpendicular(-headingXY[1], headingXY[0], 0.0f);

    // Place center to the right for clockwise turn, left for counterclockwise
    gp_vec3 center = start + rightPerpendicular * (clockwise ? 1.0f : -1.0f) * radius;

    // Starting angle: we're at the start position relative to center
    gp_scalar startAngle = std::atan2(start[1] - center[1], start[0] - center[0]);

    // When looking down from above (down +z axis), clockwise motion means INCREASING angle
    gp_scalar angleSign = clockwise ? 1.0f : -1.0f;

    for (gp_scalar angle = step; angle <= angleRadians; angle += step) {
      gp_scalar totalAngle = startAngle + angleSign * angle;
      gp_vec3 point(center[0] + radius * std::cos(totalAngle), center[1] + radius * std::sin(totalAngle), start[2]);

      gp_scalar segmentDist = (point - path.back().start).norm();
      totalDistance += segmentDist;

      gp_scalar simTimeMsec = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0f;
      path.push_back(Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsec));
    }
  }

  void addSpiralTurn(std::vector<Path>& path, const gp_vec3& start, gp_scalar radius,
                     gp_scalar angleRadians, bool clockwise, gp_scalar totalClimb, gp_scalar& totalDistance) {
    const gp_scalar step = 0.05f;  // 0.05 rad (~3°) spacing - baseline from horizontal 8

    gp_vec3 heading(-1.0f, 0.0f, 0.0f);
    if (path.size() >= 2) {
      heading = (path.back().start - path[path.size()-2].start).normalized();
    }

    // Project heading onto XY plane (ignore Z component for turn calculation)
    gp_vec3 headingXY(heading[0], heading[1], 0.0f);
    headingXY = headingXY.normalized();

    // Right perpendicular in xy plane: rotate (x,y) by -90° = (y, -x)
    gp_vec3 rightPerpendicular(-headingXY[1], headingXY[0], 0.0f);

    // Place center based on turn direction
    gp_vec3 center = start + rightPerpendicular * (clockwise ? 1.0f : -1.0f) * radius;

    // Starting angle from center to start position
    gp_scalar startAngle = std::atan2(start[1] - center[1], start[0] - center[0]);

    // When looking down from above, clockwise means increasing angle
    gp_scalar angleSign = clockwise ? 1.0f : -1.0f;

    for (gp_scalar angle = step; angle <= angleRadians; angle += step) {
      gp_scalar totalAngle = startAngle + angleSign * angle;
      gp_scalar zOffset = (angle / angleRadians) * totalClimb;
      gp_vec3 point(center[0] + radius * std::cos(totalAngle), center[1] + radius * std::sin(totalAngle), start[2] + zOffset);

      gp_scalar segmentDist = (point - path.back().start).norm();
      totalDistance += segmentDist;

      gp_scalar simTimeMsec = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0f;
      path.push_back(Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsec));
    }
  }

  void addHorizontalLoop(std::vector<Path>& path, const gp_vec3& loopOrigin, gp_scalar loopRadius,
                         bool clockwise, gp_scalar& totalDistance) {
    const gp_scalar step = 0.05f;  // 0.05 rad (~3°) spacing - baseline from horizontal 8
    gp_scalar sign = clockwise ? -1.0f : 1.0f; // inverted for consistency

    for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += step) {
      gp_vec3 circleCenter = loopOrigin + gp_vec3(0.0f, sign * loopRadius, 0.0f);
      gp_vec3 point = circleCenter + gp_vec3(-loopRadius * std::sin(turn), -sign * loopRadius * std::cos(turn), 0.0f);

      if (!path.empty()) {
        gp_scalar distance = (point - path.back().start).norm();
        totalDistance += distance;
      }
      gp_scalar simTimeMsec = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0f;
      path.push_back(Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsec));
    }
  }

  void addPitchDownLoop(std::vector<Path>& path, const gp_vec3& start, const gp_vec3& heading,
                        gp_scalar loopRadius, gp_scalar& totalDistance) {
    // Pitch-down loop: 180° loop in the vertical plane perpendicular to current heading
    // This is a Split-S maneuver - half loop downward (toward +z) that reverses course
    const gp_scalar step = 0.05f;  // 0.05 rad (~3°) spacing - baseline from horizontal 8

    // Heading should be normalized (only xy components matter for horizontal heading)
    gp_vec3 headingNorm = heading.normalized();

    // Find the center of the loop: radius distance DOWN (+z direction) from start position
    // In NED, +z is down, so we move the center downward
    gp_vec3 center = start + gp_vec3(0.0f, 0.0f, loopRadius);

    // Loop from 0 to π (180°), starting from top of loop heading down
    for (gp_scalar angle = step; angle <= static_cast<gp_scalar>(M_PI); angle += step) {
      // In the vertical plane aligned with heading:
      // - At angle=0: at top of loop (start position)
      // - At angle=π/2: halfway, heading straight down
      // - At angle=π: at bottom, reversed heading

      // Position along the loop arc
      // x,y offset: moves forward along heading direction as we loop
      gp_scalar horizontalOffset = loopRadius * std::sin(angle);  // 0 at top, r at middle, 0 at bottom
      // z offset: loops downward (toward +z)
      gp_scalar verticalOffset = -loopRadius * std::cos(angle);  // -r at top (angle=0), +r at bottom (angle=π)

      // Combine: move along heading (xy) and vertically (z)
      gp_vec3 point = center + headingNorm * horizontalOffset + gp_vec3(0.0f, 0.0f, verticalOffset);

      gp_scalar segmentDist = (point - path.back().start).norm();
      totalDistance += segmentDist;

      gp_scalar simTimeMsec = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0f;
      path.push_back(Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsec));
    }
  }
};

class GenerateLongSequential : public GeneratorMethod {
public:
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base, unsigned int seed) override {
    // LongSequential path is deterministic, seed is unused
    std::vector<Path> longPath;
    gp_scalar totalDistance = 0.0f;

    // Origin point and radius
    gp_vec3 origin(0.0f, 0.0f, SIM_INITIAL_ALTITUDE);
    gp_scalar loopRadius = static_cast<gp_scalar>(20.0f);

    // 0. Lead-in: ~1s straight-and-level southbound to let MSP override settle
    const gp_scalar leadSeconds = static_cast<gp_scalar>(1.0f);
    const gp_scalar leadDistance = SIM_RABBIT_VELOCITY * leadSeconds;
    const int leadSteps = 20;
    for (int i = 0; i <= leadSteps; ++i) {
      gp_scalar frac = static_cast<gp_scalar>(i) / static_cast<gp_scalar>(leadSteps);
      gp_vec3 point = origin + gp_vec3(-leadDistance * frac, 0.0f, 0.0f); // move south (negative X)
      gp_scalar distance = (longPath.empty() ? 0.0f : (point - longPath.back().start).norm());
      totalDistance += distance;
      gp_scalar simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
      Path pathSegment = Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsecLocal);
      longPath.push_back(pathSegment);
    }

    // Use the end of the lead-in as the loop origin so we start the horizontal 8 from the current spot
    gp_vec3 loopOrigin = longPath.empty() ? origin : longPath.back().start;

    // 2. LEFT HORIZONTAL LOOP - Start at origin heading south, turn left (counter-clockwise)
    // Circle center to the west of origin so we turn left around it  
    for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
      gp_vec3 circleCenter = loopOrigin + gp_vec3(0.0f, -loopRadius, 0.0f);
      // Start at origin (turn=0): point at center + (-radius, 0, 0), heading south
      gp_vec3 point = circleCenter + gp_vec3(-loopRadius * std::sin(turn), loopRadius * std::cos(turn), 0.0f);
      gp_scalar distance = (longPath.empty() ? 0.0f : (point - longPath.back().start).norm());
      totalDistance += distance;
      gp_scalar simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
      Path pathSegment = Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsecLocal);
      longPath.push_back(pathSegment);
    }

    
    // 3. RIGHT HORIZONTAL LOOP - Start at origin heading south, turn right (clockwise)  
    // Circle center to the east of origin so we turn right around it
    for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
      gp_vec3 circleCenter = loopOrigin + gp_vec3(0.0f, loopRadius, 0.0f);
      // Start at origin (turn=π): point at center + (-radius, 0, 0), heading south
      gp_vec3 point = circleCenter + gp_vec3(-loopRadius * std::sin(turn), -loopRadius * std::cos(turn), 0.0f);
      gp_scalar distance = (point - longPath.back().start).norm();
      totalDistance += distance;
      gp_scalar simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
      Path pathSegment = Path(point, gp_vec3::UnitX(), totalDistance, 0.0f, simTimeMsecLocal);
      longPath.push_back(pathSegment);
    }

    return longPath;
  }
};


#endif
