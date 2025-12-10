// Define a structure for 3D points
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>

#include "minisim.h"
#include "gp_types.h"

#define NUM_SEGMENTS_PER_PATH 16

gp_vec3 randomPointInCylinder(gp_scalar radius, gp_scalar height, gp_scalar base = SIM_INITIAL_ALTITUDE);
gp_vec3 cubicInterpolate(const gp_vec3& p0, const gp_vec3& p1, const gp_vec3& p2, const gp_vec3& p3, gp_scalar t);
std::vector<std::vector<Path>> generateSmoothPaths(char* method, int numPaths, gp_scalar radius, gp_scalar height);

class GeneratorMethod {
public:
  virtual ~GeneratorMethod() = default;
  virtual std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base) = 0;
};

class GenerateRandom : public GeneratorMethod {
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base) override {
    std::vector<gp_vec3> controlPoints;
    std::vector<Path> path;
    int numPoints = NUM_SEGMENTS_PER_PATH;

    // Generate random control points
    for (size_t i = 0; i < numPoints; ++i) {
      controlPoints.push_back(randomPointInCylinder(radius, height));
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
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base) override {
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
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base) override {
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
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base) override {

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

class GenerateLongSequential : public GeneratorMethod {
public:
  std::vector<Path> method(int pathIndex, gp_scalar radius, gp_scalar height, gp_scalar base) override {
    std::vector<Path> longPath;
    gp_scalar totalDistance = 0.0f;

    // Origin point and radius
    gp_vec3 origin(0.0f, 0.0f, SIM_INITIAL_ALTITUDE);
    gp_scalar loopRadius = static_cast<gp_scalar>(20.0f);

    // 2. LEFT HORIZONTAL LOOP - Start at origin heading south, turn left (counter-clockwise)
    // Circle center to the west of origin so we turn left around it  
    for (gp_scalar turn = 0; turn < static_cast<gp_scalar>(M_PI * 2.0); turn += static_cast<gp_scalar>(0.05f)) {
      gp_vec3 circleCenter = origin + gp_vec3(0.0f, -loopRadius, 0.0f);
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
      gp_vec3 circleCenter = origin + gp_vec3(0.0f, loopRadius, 0.0f);
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
