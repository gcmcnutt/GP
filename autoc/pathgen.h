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

          // add next segment            
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), odometer, turnmeter };
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

          // add next segment            
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), odometer, turnmeter };
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
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), turn * radius, 0.0 };
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::VerticalCircle: {
        double cSize = 0.8 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { cSize * cos(turn), 0.0, base - SIM_INITIAL_LOCATION_DITHER - cSize * sin(turn)};
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), turn * cSize, 0.0 };
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::FigureEight: {
        double cSize = 0.4 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { (cSize * cos(turn)) - cSize, 0.0, base - SIM_INITIAL_LOCATION_DITHER - (cSize * sin(turn))};
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), (turn * cSize), 0.0 };
          path.push_back(pathSegment);
        }
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { (-cSize * cos(turn)) + cSize, 0.0, base - SIM_INITIAL_LOCATION_DITHER - (cSize * sin(turn))};
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), (turn * cSize) + (M_PI * 2 * cSize), 0.0 };
          path.push_back(pathSegment);
        }
        break;
      }
      case PathType::HorizontalSquare: {
        double sSize = 0.7 * radius;
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { span, -sSize, base - SIM_INITIAL_LOCATION_DITHER };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), span, 0.0 };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { sSize, span, base - SIM_INITIAL_LOCATION_DITHER };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 2 * sSize + span, 0.0 };
          path.push_back(pathSegment);
        } 
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -span, sSize, base - SIM_INITIAL_LOCATION_DITHER };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 4 * sSize + span, 0.0 };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -sSize, -span, base - SIM_INITIAL_LOCATION_DITHER };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 6 * sSize + span, 0.0 };
          path.push_back(pathSegment);
        }      
        break;
      }
      case PathType::VerticalSquare: {
        double sSize = 0.7 * radius;
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { span, 0, base - SIM_INITIAL_LOCATION_DITHER - sSize };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), span, 0.0 };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { sSize, 0, base - SIM_INITIAL_LOCATION_DITHER + span };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 2 * sSize + span, 0.0 };
          path.push_back(pathSegment);
        } 
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -span, 0, base - SIM_INITIAL_LOCATION_DITHER + sSize };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 4 * sSize + span, 0.0 };
          path.push_back(pathSegment);
        }  
        for (double span = -sSize; span < sSize; span += 0.1) {
          Eigen::Vector3d interpolatedPoint = { -sSize, 0, base - SIM_INITIAL_LOCATION_DITHER -span };
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), 6 * sSize + span, 0.0 };
          path.push_back(pathSegment);
        }      
        break;
      }
      case PathType::FortyFiveLoop: {
        double cSize = 0.8 * radius;
        for (double turn = 0; turn < M_PI * 2; turn += 0.05) {
          Eigen::Vector3d interpolatedPoint = { cSize * cos(turn), cSize * sin(turn), base - SIM_INITIAL_LOCATION_DITHER - cSize * sin(turn)};
          Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), turn * cSize, 0.0 };
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
    for (double i = 0; i <= SIM_PATH_RADIUS_LIMIT; i += 5) {
      Eigen::Vector3d interpolatedPoint = { i, i, SIM_INITIAL_ALTITUDE };
      Path pathSegment = { interpolatedPoint, Eigen::Vector3d::UnitX(), distance, 0.0 };
      path.push_back(pathSegment);
      distance += 5;
    }
    return path;
  }
};


#endif