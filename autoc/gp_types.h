#pragma once

// Portable Eigen includes — desktop uses system Eigen, Arduino uses bundled
#ifdef ARDUINO
#include "ArduinoEigenDense.h"
#else
#include <Eigen/Dense>
#include <Eigen/Geometry>
#endif

// Shared scalar + Eigen aliases; float-only path for eval performance
using gp_scalar = float;
using gp_vec3 = Eigen::Matrix<gp_scalar, 3, 1>;
using gp_quat = Eigen::Quaternion<gp_scalar>;

// Double precision for fitness accumulation during training/eval
using gp_fitness = double;
