#include "vmath.h"

#include <cmath>

namespace vmath {
// Multiply two 3x3 matrices
Matrix3x3 multiplyMatrices(const Matrix3x3& A, const Matrix3x3& B) {
    Matrix3x3 C = {0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Multiply a 3x3 matrix by a 3D vector
std::array<double, 3> multiplyMatrixVector(const Matrix3x3& M, const std::array<double, 3>& v) {
    std::array<double, 3> result = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        result[i] = M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2];
    }
    return result;
}

// Create a yaw rotation matrix
Matrix3x3 getYawMatrix(double yaw) {
    return {{
        {cos(yaw), -sin(yaw), 0},
        {sin(yaw), cos(yaw), 0},
        {0, 0, 1}
    }};
}

// Create a pitch rotation matrix
Matrix3x3 getPitchMatrix(double pitch) {
    return {{
        {cos(pitch), 0, sin(pitch)},
        {0, 1, 0},
        {-sin(pitch), 0, cos(pitch)}
    }};
}

// Create a roll rotation matrix
Matrix3x3 getRollMatrix(double roll) {
    return {{
        {1, 0, 0},
        {0, cos(roll), -sin(roll)},
        {0, sin(roll), cos(roll)}
    }};
}

// Create an inverse yaw rotation matrix
Matrix3x3 getInverseYawMatrix(double yaw) {
    return {{
        {cos(yaw), sin(yaw), 0},
        {-sin(yaw), cos(yaw), 0},
        {0, 0, 1}
    }};
}

// Create an inverse pitch rotation matrix
Matrix3x3 getInversePitchMatrix(double pitch) {
    return {{
        {cos(pitch), 0, -sin(pitch)},
        {0, 1, 0},
        {sin(pitch), 0, cos(pitch)}
    }};
}

// Create an inverse roll rotation matrix
Matrix3x3 getInverseRollMatrix(double roll) {
    return {{
        {1, 0, 0},
        {0, cos(roll), sin(roll)},
        {0, -sin(roll), cos(roll)}
    }};
}

// Convert world coordinates to body frame
Matrix3x3 worldToBody(double yaw, double pitch, double roll) {
    // Get the individual rotation matrices
    Matrix3x3 Rz = getYawMatrix(yaw);
    Matrix3x3 Ry = getPitchMatrix(pitch);
    Matrix3x3 Rx = getRollMatrix(roll);

    // Combine the rotation matrices into one
    Matrix3x3 R = multiplyMatrices(Rx, Ry);
    R = multiplyMatrices(R, Rz);

    return R;
}

// Convert body frame coordinates to world coordinates
Matrix3x3 bodyToWorld(double yaw, double pitch, double roll) {
    // Get the individual rotation matrices
    Matrix3x3 Rz = getYawMatrix(yaw);
    Matrix3x3 Ry = getPitchMatrix(pitch);
    Matrix3x3 Rx = getRollMatrix(roll);

    // Combine the rotation matrices into one
    Matrix3x3 R = multiplyMatrices(Rz, Ry);
    R = multiplyMatrices(R, Rx);

    return R;
}

// Convert world coordinates to body coordinates
std::array<double, 3> worldToLocal(double yaw, double pitch, double roll, const std::array<double, 3>& worldCoords) {
    // Get the individual inverse rotation matrices
    Matrix3x3 Rz_inv = getInverseYawMatrix(yaw);
    Matrix3x3 Ry_inv = getInversePitchMatrix(pitch);
    Matrix3x3 Rx_inv = getInverseRollMatrix(roll);

    // Combine the inverse rotation matrices into one
    Matrix3x3 R_inv = multiplyMatrices(Rx_inv, Ry_inv);
    R_inv = multiplyMatrices(R_inv, Rz_inv);

    // Transform the world coordinates to body coordinates
    return multiplyMatrixVector(R_inv, worldCoords);
}

// Convert local (body frame) coordinates to world coordinates
std::array<double, 3> localToWorld(double yaw, double pitch, double roll, const std::array<double, 3>& localCoords) {
    // Get the individual rotation matrices
    Matrix3x3 Rz = getYawMatrix(yaw);
    Matrix3x3 Ry = getPitchMatrix(pitch);
    Matrix3x3 Rx = getRollMatrix(roll);

    // Combine the rotation matrices into one
    Matrix3x3 R = multiplyMatrices(Rz, Ry);
    R = multiplyMatrices(R, Rx);

    // Transform the local coordinates to world coordinates
    return multiplyMatrixVector(R, localCoords);
}

// Function to convert Euler angles to a rotation matrix
Matrix3x3 eulerToRotationMatrix(double yaw, double pitch, double roll) {
    Matrix3x3 Rz = getYawMatrix(yaw);
    Matrix3x3 Ry = getPitchMatrix(pitch);
    Matrix3x3 Rx = getRollMatrix(roll);
    
    Matrix3x3 R = multiplyMatrices(Rz, Ry);
    R = multiplyMatrices(R, Rx);
    
    return R;
}

// Function to create a rotation matrix for pitch and roll adjustments
Matrix3x3 getAdjustmentMatrix(double deltaPitch, double deltaRoll) {
    Matrix3x3 Ry = getPitchMatrix(deltaPitch);
    Matrix3x3 Rx = getRollMatrix(deltaRoll);
    
    return multiplyMatrices(Rx, Ry);
}

// Function to combine the original rotation matrix with the adjustments
Matrix3x3 applyAdjustments(const Matrix3x3& originalRotation, double deltaPitch, double deltaRoll) {
    Matrix3x3 adjustmentMatrix = getAdjustmentMatrix(deltaPitch, deltaRoll);
    
    return multiplyMatrices(originalRotation, adjustmentMatrix);
}

// Function to convert a rotation matrix back to Euler angles
std::array<double, 3> rotationMatrixToEuler(const Matrix3x3& R) {
    double sy = sqrt(R[0][0] * R[0][0] + R[1][0] * R[1][0]);
    
    bool singular = sy < 1e-6; // If sy is close to zero, consider it singular

    double yaw, pitch, roll;
    if (!singular) {
        yaw = atan2(R[1][0], R[0][0]);
        pitch = atan2(-R[2][0], sy);
        roll = atan2(R[2][1], R[2][2]);
    } else {
        yaw = atan2(-R[1][2], R[1][1]);
        pitch = atan2(-R[2][0], sy);
        roll = 0;
    }
    
    return {yaw, pitch, roll};
}

} // namespace vmath

