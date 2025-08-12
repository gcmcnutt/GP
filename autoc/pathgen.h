// Define a structure for 3D points
#ifndef PATHGEN_H
#define PATHGEN_H

#include <vector>

#include "minisim.h"

#define NUM_SEGMENTS_PER_PATH 16

Eigen::Vector3d randomPointInCylinder(double radius, double height, double base = SIM_INITIAL_ALTITUDE);
Eigen::Vector3d cubicInterpolate(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double t);
std::vector<std::vector<Path>> generateSmoothPaths(char* method, int numPaths, double radius, double height);

class GeneratorMethod {
public:
  virtual ~GeneratorMethod() = default;
  virtual std::vector<Path> method(int pathIndex, double radius, double height, double base) = 0;
};

class GenerateRandom : public GeneratorMethod {
  std::vector<Path> method(int pathIndex, double radius, double height, double base) override {
    std::vector<Eigen::Vector3d> controlPoints;
    std::vector<Path> path;
    int numPoints = NUM_SEGMENTS_PER_PATH;

    // Generate random control points
    for (size_t i = 0; i < numPoints; ++i) {
      controlPoints.push_back(randomPointInCylinder(radius, height));
    }

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

          // add next segment with simulation timestamp
          double simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), odometer, turnmeter, simTimeMsec };
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
  std::vector<Path> method(int pathIndex, double radius, double height, double base) override {
    std::vector<Eigen::Vector3d> controlPoints;
    std::vector<Path> path;
    int numPoints = NUM_SEGMENTS_PER_PATH;

    controlPoints.push_back({ 0, 0, base });

    // Sin
    double x, y, z;
    for (size_t i = 0; i < numPoints; ++i) {
      x = -(cos(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS / 2 - SIM_PATH_BOUNDS / 2);
      y = sin(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS / 2;
      z = base - i;
      controlPoints.push_back(Eigen::Vector3d(x, y, z));
    }

    double lastX = x;
    double lastY = y;
    double lastZ = z;
    for (size_t i = 0; i < numPoints; ++i) {
      y = lastX + sin(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS / 3;
      z = lastZ - SIM_PATH_BOUNDS / 2 + cos(2 * M_PI * i / numPoints) * SIM_PATH_BOUNDS / 2;
      x = lastY - i;
      controlPoints.push_back(Eigen::Vector3d(x, y, z));
    }

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

          // add next segment with simulation timestamp
          double simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), odometer, turnmeter, simTimeMsec };
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
  std::vector<Path> method(int pathIndex, double radius, double height, double base) override {
    std::vector<Path> path;

    PathType pathType = static_cast<PathType>(pathIndex % static_cast<int>(PathType::END_MARKER));

    switch(pathType) {
      case PathType::HorizontalCircle: {
        double cSize = 0.8 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { cSize * cos(turn), cSize * sin(turn), base - SIM_INITIAL_LOCATION_DITHER };
          double simTimeMsecLocal = (turn * radius / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), turn * radius, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::VerticalCircle: {
        double cSize = 0.8 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { cSize * cos(turn), 0.0, base - SIM_INITIAL_LOCATION_DITHER - cSize * sin(turn)};
          double simTimeMsecLocal = (turn * cSize / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), turn * cSize, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::FigureEight: {
        double cSize = 0.4 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { (cSize * cos(turn)) - cSize, 0.0, base - SIM_INITIAL_LOCATION_DITHER - (cSize * sin(turn))};
          double simTimeMsecLocal = ((turn * cSize) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), (turn * cSize), 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { (-cSize * cos(turn)) + cSize, 0.0, base - SIM_INITIAL_LOCATION_DITHER - (cSize * sin(turn))};
          double simTimeMsecLocal = (((turn * cSize) + (M_PI * 2 * cSize)) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), (turn * cSize) + (M_PI * 2 * cSize), 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::HorizontalSquare: {
        double sSize = 0.7 * radius;
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { span, -sSize, base - SIM_INITIAL_LOCATION_DITHER };
          double simTimeMsecLocal = (span / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { sSize, span, base - SIM_INITIAL_LOCATION_DITHER };
          double simTimeMsecLocal = ((2 * sSize + span) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 2 * sSize + span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        } 
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -span, sSize, base - SIM_INITIAL_LOCATION_DITHER };
          double simTimeMsecLocal = ((4 * sSize + span) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 4 * sSize + span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -sSize, -span, base - SIM_INITIAL_LOCATION_DITHER };
          double simTimeMsecLocal = ((6 * sSize + span) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 6 * sSize + span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }      
        break;
      }
      case PathType::VerticalSquare: {
        double sSize = 0.7 * radius;
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { span, 0, base - SIM_INITIAL_LOCATION_DITHER - sSize };
          double simTimeMsecLocal = (span / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { sSize, 0, base - SIM_INITIAL_LOCATION_DITHER + span };
          double simTimeMsecLocal = ((2 * sSize + span) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 2 * sSize + span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        } 
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -span, 0, base - SIM_INITIAL_LOCATION_DITHER + sSize };
          double simTimeMsecLocal = ((4 * sSize + span) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 4 * sSize + span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -sSize, 0, base - SIM_INITIAL_LOCATION_DITHER -span };
          double simTimeMsecLocal = ((6 * sSize + span) / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 6 * sSize + span, 0.0, simTimeMsecLocal };
          path.push_back(pathSegment);
        }      
        break;
      }
      case PathType::FortyFiveLoop: {
        double cSize = 0.8 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { cSize * cos(turn), cSize * sin(turn), base - SIM_INITIAL_LOCATION_DITHER - cSize * sin(turn)};
          double simTimeMsecLocal = (turn * cSize / SIM_RABBIT_VELOCITY) * 1000.0;
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), turn * cSize, 0.0, simTimeMsecLocal };
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
  std::vector<Path> method(int pathIndex, double radius, double height, double base) override {

    std::vector<Path> path;
    double distance = 0.0;
    for (double i = -SIM_PATH_RADIUS_LIMIT; i <= SIM_PATH_RADIUS_LIMIT; i += 5) {
      Eigen::Vector3d interpolatedPoint = { -30, i, SIM_INITIAL_ALTITUDE };
      double simTimeMsecLocal = (distance / SIM_RABBIT_VELOCITY) * 1000.0;
      Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), distance, 0.0, simTimeMsecLocal };
      path.push_back(pathSegment);
      distance += 5;
    }
    return path;
  }
};

class GenerateLongSequential : public GeneratorMethod {
public:
  std::vector<Path> method(int pathIndex, double radius, double height, double base) override {
    std::vector<Path> longPath;
    double totalDistance = 0.0;

    // Origin point and radius
    Eigen::Vector3d origin(0, 0, SIM_INITIAL_ALTITUDE);
    double loopRadius = 20.0;

    // // 1. VERTICAL LOOP - Enter at origin heading south, pitch up through full loop
    // // Start at origin, loop goes down first then up and around
    // for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
    //   Eigen::Vector3d circleCenter = origin + Eigen::Vector3d(0, 0, -loopRadius);
    //   // Start at origin (turn=0): point at center + (0, 0, radius), heading south
    //   Eigen::Vector3d point = circleCenter + Eigen::Vector3d(-loopRadius * sin(turn), 0, loopRadius * cos(turn));
    //   double distance = (longPath.empty() ? 0 : (point - longPath.back().start).norm());
    //   totalDistance += distance;
    //   double simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0;
    //   Path pathSegment = { point, Eigen::Vector3d::UnitX(), totalDistance, 0.0, simTimeMsecLocal };
    //   longPath.push_back(pathSegment);
    // }

    // 2. LEFT HORIZONTAL LOOP - Start at origin heading south, turn left (counter-clockwise)
    // Circle center to the west of origin so we turn left around it  
    for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
      Eigen::Vector3d circleCenter = origin + Eigen::Vector3d(0, -loopRadius, 0);
      // Start at origin (turn=0): point at center + (-radius, 0, 0), heading south
      Eigen::Vector3d point = circleCenter + Eigen::Vector3d(-loopRadius * sin(turn), loopRadius * cos(turn), 0);
      double distance = (longPath.empty() ? 0 : (point - longPath.back().start).norm());
      totalDistance += distance;
      double simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0;
      Path pathSegment = { point, Eigen::Vector3d::UnitX(), totalDistance, 0.0, simTimeMsecLocal };
      longPath.push_back(pathSegment);
    }

    
    // 3. RIGHT HORIZONTAL LOOP - Start at origin heading south, turn right (clockwise)  
    // Circle center to the east of origin so we turn right around it
    for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
      Eigen::Vector3d circleCenter = origin + Eigen::Vector3d(0, loopRadius, 0);
      // Start at origin (turn=π): point at center + (-radius, 0, 0), heading south
      Eigen::Vector3d point = circleCenter + Eigen::Vector3d(-loopRadius * sin(turn), -loopRadius * cos(turn), 0);
      double distance = (point - longPath.back().start).norm();
      totalDistance += distance;
      double simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0;
      Path pathSegment = { point, Eigen::Vector3d::UnitX(), totalDistance, 0.0, simTimeMsecLocal };
      longPath.push_back(pathSegment);
    }

    // STILL COMMENTED OUT: Keep remaining loops for later testing
    /*
    // 4. RIGHT TILT LOOP - Roll right 45°, pitch up through tilted circle
    // Same as vertical loop is good enough, not rotated just stretched in y
    for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
      Eigen::Vector3d circleCenter = origin + Eigen::Vector3d(0, loopRadius, -loopRadius);
      // Start at origin (turn=π): point at center + (0, 0, radius), heading south
      Eigen::Vector3d point = circleCenter + Eigen::Vector3d(-loopRadius * sin(turn), -loopRadius * cos(turn), loopRadius * cos(turn));
      double distance = (longPath.empty() ? 0 : (point - longPath.back().start).norm());
      totalDistance += distance;
      double simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0;
      Path pathSegment = { point, Eigen::Vector3d::UnitX(), totalDistance, 0.0, simTimeMsecLocal };
      longPath.push_back(pathSegment);
    }

    // 5. LEFT TILT LOOP - Roll left 45°, pitch up through tilted circle
    // Start at origin, perform vertical loop tilted 45° to the left (toward -Y/west)
    for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
      Eigen::Vector3d circleCenter = origin + Eigen::Vector3d(0, -loopRadius, -loopRadius);
      // Start at origin (turn=π): point at center + (0, 0, radius), heading south
      Eigen::Vector3d point = circleCenter + Eigen::Vector3d(-loopRadius * sin(turn), loopRadius * cos(turn), loopRadius * cos(turn));
      double distance = (longPath.empty() ? 0 : (point - longPath.back().start).norm());
      totalDistance += distance;
      double simTimeMsecLocal = (totalDistance / SIM_RABBIT_VELOCITY) * 1000.0;
      Path pathSegment = { point, Eigen::Vector3d::UnitX(), totalDistance, 0.0, simTimeMsecLocal };
      longPath.push_back(pathSegment);
    }
    */

    return longPath;
  }
};


#endif