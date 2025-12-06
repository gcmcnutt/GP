#ifndef FASTMATH_ORIENTATION_MATH_H
#define FASTMATH_ORIENTATION_MATH_H

#include <cmath>
#include <Eigen/Geometry>

namespace fastmath {

inline Eigen::Vector3d rotateVectorRaw(double w, double x, double y, double z, const Eigen::Vector3d& v) {
    double tx = 2.0 * (y * v.z() - z * v.y());
    double ty = 2.0 * (z * v.x() - x * v.z());
    double tz = 2.0 * (x * v.y() - y * v.x());

    Eigen::Vector3d result;
    result.x() = v.x() + w * tx + (y * tz - z * ty);
    result.y() = v.y() + w * ty + (z * tx - x * tz);
    result.z() = v.z() + w * tz + (x * ty - y * tx);
    return result;
}

inline Eigen::Vector3d rotateWorldToBody(const Eigen::Quaterniond& q, const Eigen::Vector3d& v) {
    return rotateVectorRaw(q.w(), -q.x(), -q.y(), -q.z(), v);
}

inline double rollFromQuaternion(const Eigen::Quaterniond& q) {
    double sinr_cosp = 2.0 * (q.w() * q.x() + q.y() * q.z());
    double cosr_cosp = 1.0 - 2.0 * (q.x() * q.x() + q.y() * q.y());
    return std::atan2(sinr_cosp, cosr_cosp);
}

inline double pitchFromQuaternion(const Eigen::Quaterniond& q) {
    double sinp = 2.0 * (q.w() * q.y() - q.z() * q.x());
    if (std::abs(sinp) >= 1.0) {
        return std::copysign(M_PI / 2.0, sinp);
    }
    return std::asin(sinp);
}

}  // namespace fastmath

#endif  // FASTMATH_ORIENTATION_MATH_H
