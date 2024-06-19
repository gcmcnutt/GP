#ifndef VMATH_H
#define VMATH_H

#include <array>

namespace vmath
{
// Define a 3x3 matrix type
using Matrix3x3 = std::array<std::array<double, 3>, 3>;

// Multiply two 3x3 matrices
Matrix3x3 multiplyMatrices(const Matrix3x3 &A, const Matrix3x3 &B);

// Multiply a 3x3 matrix by a 3D vector
std::array<double, 3> multiplyMatrixVector(const Matrix3x3 &M, const std::array<double, 3> &v);

// Create a yaw rotation matrix
Matrix3x3 getYawMatrix(double yaw);

// Create a pitch rotation matrix
Matrix3x3 getPitchMatrix(double pitch);

// Create a roll rotation matrix
Matrix3x3 getRollMatrix(double roll);

// Create an inverse yaw rotation matrix
Matrix3x3 getInverseYawMatrix(double yaw);

// Create an inverse pitch rotation matrix
Matrix3x3 getInversePitchMatrix(double pitch);

// Create an inverse roll rotation matrix
Matrix3x3 getInverseRollMatrix(double roll);

// Convert world coordinates to body coordinates
std::array<double, 3> worldToLocal(double yaw, double pitch, double roll, const std::array<double, 3> &worldCoords);

// Convert local (body frame) coordinates to world coordinates
std::array<double, 3> localToWorld(double yaw, double pitch, double roll, const std::array<double, 3>& localCoords);

// Convert world coordinates to body frame
Matrix3x3 worldToBody(double yaw, double pitch, double roll);

// Convert body frame coordinates to world coordinates
Matrix3x3 bodyToWorld(double yaw, double pitch, double roll);

// Function to convert Euler angles to a rotation matrix
Matrix3x3 eulerToRotationMatrix(double yaw, double pitch, double roll);

// Function to create a rotation matrix for pitch and roll adjustments
Matrix3x3 getAdjustmentMatrix(double deltaPitch, double deltaRoll);

// Function to combine the original rotation matrix with the adjustments
Matrix3x3 applyAdjustments(const Matrix3x3& originalRotation, double deltaPitch, double deltaRoll);

// Function to convert a rotation matrix back to Euler angles
std::array<double, 3> rotationMatrixToEuler(const Matrix3x3& R);

} // namespace vmath

#endif // VMATH_H
