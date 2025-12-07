#pragma once

#ifdef GP_BUILD
#include <Eigen/Dense>
#include <Eigen/Geometry>
#else
#include <ArduinoEigenDense.h>
#endif

// Shared scalar + Eigen aliases for GP eval; float-only path
using gp_scalar = float;
using gp_vec3 = Eigen::Matrix<gp_scalar, 3, 1>;
using gp_quat = Eigen::Quaternion<gp_scalar>;
