#ifndef BLACKBOX_EXTRACTOR_SIMPLE_H
#define BLACKBOX_EXTRACTOR_SIMPLE_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>

/*
 * Simple blackbox data extractor for iNav flight logs
 * Extracts position and attitude data for rendering
 */
struct BlackboxRecord {
    unsigned long timestamp_ms;
    double latitude;
    double longitude;
    double altitude;
    double roll;
    double pitch;
    double yaw;
    double throttle;
    double roll_command;
    double pitch_command;
};

struct SimpleAircraftState {
    Eigen::Vector3d position;
    Eigen::Quaterniond orientation;
    double relativeVelocity;
    double pitchCommand;
    double rollCommand;
    double throttleCommand;
    unsigned long timeMsec;
    
    SimpleAircraftState(const Eigen::Vector3d& pos, const Eigen::Quaterniond& orient, 
                       double vel, double pitch, double roll, double throttle, unsigned long time)
        : position(pos), orientation(orient), relativeVelocity(vel), 
          pitchCommand(pitch), rollCommand(roll), throttleCommand(throttle), timeMsec(time) {}
};

class BlackboxExtractor {
public:
    BlackboxExtractor();
    
    // Parse blackbox log file and extract flight data
    bool loadBlackboxFile(const std::string& filename);
    
    // Convert blackbox data to SimpleAircraftState format
    std::vector<SimpleAircraftState> extractAircraftStates();
    
    // Get the reference position (first GPS point) for centering
    Eigen::Vector3d getReferencePosition() const;
    
    // Get flight duration in milliseconds
    unsigned long getFlightDuration() const;
    
    // Get number of records
    size_t getRecordCount() const;
    
private:
    std::vector<BlackboxRecord> records;
    Eigen::Vector3d referencePosition;
    bool hasReferencePosition;
    
    // Helper functions
    bool parseBlackboxLine(const std::string& line, BlackboxRecord& record);
    Eigen::Vector3d gpsToLocal(double lat, double lon, double alt) const;
    Eigen::Quaterniond eulerToQuaternion(double roll, double pitch, double yaw) const;
    
    // GPS coordinate conversion constants
    static constexpr double DEG_TO_RAD = M_PI / 180.0;
    static constexpr double EARTH_RADIUS = 6371000.0; // meters
};

#endif