#include "blackbox_extractor_simple.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

BlackboxExtractor::BlackboxExtractor() 
    : hasReferencePosition(false), referencePosition(0, 0, 0) {
}

bool BlackboxExtractor::loadBlackboxFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open blackbox file: " << filename << std::endl;
        return false;
    }
    
    records.clear();
    hasReferencePosition = false;
    
    std::string line;
    bool headerParsed = false;
    
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Parse header line (field names) - simplified for now
        if (!headerParsed && line.find("loopIteration") != std::string::npos) {
            headerParsed = true;
            continue;
        }
        
        // Parse data lines
        if (headerParsed) {
            BlackboxRecord record;
            if (parseBlackboxLine(line, record)) {
                records.push_back(record);
                
                // Set reference position from first GPS record
                if (!hasReferencePosition && record.latitude != 0 && record.longitude != 0) {
                    referencePosition = Eigen::Vector3d(record.latitude, record.longitude, record.altitude);
                    hasReferencePosition = true;
                }
            }
        }
    }
    
    file.close();
    
    std::cout << "Loaded " << records.size() << " blackbox records" << std::endl;
    if (hasReferencePosition) {
        std::cout << "Reference position: " << referencePosition.transpose() << std::endl;
    }
    
    return !records.empty();
}

bool BlackboxExtractor::parseBlackboxLine(const std::string& line, BlackboxRecord& record) {
    std::istringstream iss(line);
    std::string field;
    std::vector<std::string> fields;
    
    // Split by comma
    while (std::getline(iss, field, ',')) {
        fields.push_back(field);
    }
    
    // Basic parsing - adjust indices based on actual blackbox format
    if (fields.size() < 10) {
        return false;
    }
    
    try {
        // Parse common fields (simplified for this demo)
        record.timestamp_ms = std::stoul(fields[0]);
        record.latitude = std::stod(fields[4]) / 10000000.0;  // GPS coordinates are typically scaled
        record.longitude = std::stod(fields[5]) / 10000000.0;
        record.altitude = std::stod(fields[6]) / 100.0;       // Altitude in cm, convert to meters
        record.roll = std::stod(fields[7]) / 10.0;            // Attitude in decidegrees
        record.pitch = std::stod(fields[8]) / 10.0;
        record.yaw = std::stod(fields[9]) / 10.0;
        
        // Set defaults for missing fields
        record.throttle = 0.5;
        record.roll_command = 0.0;
        record.pitch_command = 0.0;
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing blackbox line: " << e.what() << std::endl;
        return false;
    }
}

std::vector<SimpleAircraftState> BlackboxExtractor::extractAircraftStates() {
    std::vector<SimpleAircraftState> states;
    
    if (records.empty() || !hasReferencePosition) {
        return states;
    }
    
    for (const auto& record : records) {
        // Convert GPS coordinates to local coordinates
        Eigen::Vector3d localPos = gpsToLocal(record.latitude, record.longitude, record.altitude);
        
        // Convert Euler angles to quaternion
        Eigen::Quaterniond orientation = eulerToQuaternion(
            record.roll * DEG_TO_RAD,
            record.pitch * DEG_TO_RAD,
            record.yaw * DEG_TO_RAD
        );
        
        // Create SimpleAircraftState with blackbox data
        SimpleAircraftState state(
            localPos,                    // position
            orientation,                 // aircraft orientation
            20.0,                        // default relative velocity
            record.pitch_command,        // pitch command
            record.roll_command,         // roll command
            record.throttle,             // throttle command
            record.timestamp_ms          // time in milliseconds
        );
        
        states.push_back(state);
    }
    
    return states;
}

Eigen::Vector3d BlackboxExtractor::gpsToLocal(double lat, double lon, double alt) const {
    if (!hasReferencePosition) {
        return Eigen::Vector3d(0, 0, 0);
    }
    
    // Convert GPS coordinates to local NED (North-East-Down) coordinates
    double refLat = referencePosition[0] * DEG_TO_RAD;
    double refLon = referencePosition[1] * DEG_TO_RAD;
    double refAlt = referencePosition[2];
    
    double dLat = lat * DEG_TO_RAD - refLat;
    double dLon = lon * DEG_TO_RAD - refLon;
    double dAlt = alt - refAlt;
    
    // Convert to local coordinates (approximate for small distances)
    double north = dLat * EARTH_RADIUS;
    double east = dLon * EARTH_RADIUS * std::cos(refLat);
    double down = -dAlt; // NED convention: down is positive
    
    // Convert to renderer coordinate system (X=East, Y=North, Z=Up)
    return Eigen::Vector3d(east, north, -down);
}

Eigen::Quaterniond BlackboxExtractor::eulerToQuaternion(double roll, double pitch, double yaw) const {
    // Convert Euler angles (roll, pitch, yaw) to quaternion
    // Using ZYX rotation sequence (yaw, pitch, roll)
    Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());
    
    Eigen::Quaterniond q = yawAngle * pitchAngle * rollAngle;
    return q.normalized();
}

Eigen::Vector3d BlackboxExtractor::getReferencePosition() const {
    return referencePosition;
}

unsigned long BlackboxExtractor::getFlightDuration() const {
    if (records.empty()) {
        return 0;
    }
    return records.back().timestamp_ms - records.front().timestamp_ms;
}

size_t BlackboxExtractor::getRecordCount() const {
    return records.size();
}