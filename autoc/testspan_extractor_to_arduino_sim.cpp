/*
 * Flight Log Extractor and C++ Generator for XIAO-GP
 *
 * Consolidates flight log extraction and C++ code generation into a single tool.
 * Reads INAV blackbox CSV files, extracts test spans (MSPRCOVERRIDE periods)
 * with configurable buffers, and generates Arduino-compatible C++ flight data.
 *
 * Features:
 * - Extracts test spans with MSPRCOVERRIDE flag
 * - Adds 1.5 second buffers before/after each test
 * - Handles adjacent test collision detection
 * - Selects specific test spans from multi-test files
 * - Generates optimized C++ with pre-computed MSP structures
 *
 * Usage:
 *   ./flight_extractor_cpp -i input.csv -t 1 > flight_data.cpp
 *   ./flight_extractor_cpp -i - -t 2 < input.csv > flight_data.cpp
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <algorithm>

struct CSVFrame {
    uint32_t time_us;
    float quaternion[4];        // w, x, y, z (normalized)
    int32_t navPos[3];          // lat, lon, alt
    std::string flightModeFlags;
    int16_t attitude[3];        // roll, pitch, yaw (decidegrees)
    int16_t gyroADC[3];         // roll, pitch, yaw rates
    int16_t accSmooth[3];       // x, y, z
    float vbat;
    uint16_t servo[3];          // aileron, elevator, throttle (AET)
    uint32_t mspOverrideFlags;  // Direct from CSV - 1 = override active
    bool hasMSPRCOverride;      // Derived from mspOverrideFlags
    size_t csv_line_number;     // Source line number in CSV file (for diagnostics)
};

struct TestSpan {
    size_t start_index;
    size_t end_index;
    uint32_t start_time_us;
    uint32_t end_time_us;
    size_t frame_count;
};

// Required columns that xiao-gp fetches from INAV
const std::vector<std::string> REQUIRED_COLUMNS = {
    "loopIteration", "time (us)",
    "quaternion[0]", "quaternion[1]", "quaternion[2]", "quaternion[3]",
    "navPos[0]", "navPos[1]", "navPos[2]",
    "rcData[0]", "rcData[1]", "rcData[2]", "rcData[3]",
    "flightModeFlags (flags)",
    "mspOverrideFlags",
    "attitude[0]", "attitude[1]", "attitude[2]",
    "gyroADC[0]", "gyroADC[1]", "gyroADC[2]",
    "accSmooth[0]", "accSmooth[1]", "accSmooth[2]",
    "vbat (V)"
};

const uint32_t BUFFER_TIME_US = 1500000; // 1.5 seconds in microseconds

// Command line options
struct Options {
    std::string input_file = "-";  // stdin by default
    int test_span = 1;             // first test span by default
    bool list_spans = false;       // list available spans
};

uint32_t parseFlightModeFlags(const std::string& flagStr) {
    std::map<std::string, int> flagMap = {
        {"ARM", 0}, {"ANGLE", 1}, {"HORIZON", 2}, {"NAVALTHOLD", 3},
        {"MAG", 4}, {"HEADFREE", 5}, {"HEADADJ", 6}, {"CAMSTAB", 7},
        {"NAVRTH", 8}, {"NAVPOSHOLD", 9}, {"PASSTHRU", 10}, {"BEEPERON", 11},
        {"LEDLOW", 12}, {"LLIGHTS", 13}, {"OSD", 14}, {"TELEMETRY", 15},
        {"GTUNE", 16}, {"SONAR", 17}, {"BLACKBOX", 18}, {"FAILSAFE", 19},
        {"NAVWP", 20}, {"AIRMODE", 21}, {"HOMERESET", 22}, {"GCSNAV", 23},
        {"HEADINGLOCK", 24}, {"SURFACE", 25}, {"FLAPERON", 26}, {"TURNASSIST", 27},
        {"NAVLAUNCH", 28}, {"AUTOTRIM", 29}, {"MSPRCOVERRIDE", 30}
    };

    uint32_t bitmask = 0;
    std::istringstream ss(flagStr);
    std::string flag;

    while (std::getline(ss, flag, '|')) {
        // Trim whitespace
        flag.erase(0, flag.find_first_not_of(" \t"));
        flag.erase(flag.find_last_not_of(" \t") + 1);

        auto it = flagMap.find(flag);
        if (it != flagMap.end()) {
            bitmask |= (1 << it->second);
        }
    }

    return bitmask;
}

bool hasMSPRCOverride(const std::string& flagStr) {
    return flagStr.find("MSPRCOVERRIDE") != std::string::npos;
}

std::vector<CSVFrame> parseCSV(std::istream& input) {
    std::vector<CSVFrame> frames;
    std::string line;
    std::map<std::string, int> columnMap;

    // Parse header
    if (!std::getline(input, line)) {
        std::cerr << "Error: Cannot read CSV header" << std::endl;
        return frames;
    }

    std::istringstream ss(line);
    std::string column;
    int index = 0;

    while (std::getline(ss, column, ',')) {
        // Trim whitespace
        column.erase(0, column.find_first_not_of(" \t"));
        column.erase(column.find_last_not_of(" \t") + 1);
        columnMap[column] = index++;
    }

    // Check for required columns
    std::vector<std::string> missing_columns;
    for (const auto& req_col : REQUIRED_COLUMNS) {
        if (columnMap.find(req_col) == columnMap.end()) {
            missing_columns.push_back(req_col);
        }
    }

    if (!missing_columns.empty()) {
        std::cerr << "Warning: Missing required columns: ";
        for (size_t i = 0; i < missing_columns.size(); i++) {
            std::cerr << missing_columns[i];
            if (i < missing_columns.size() - 1) std::cerr << ", ";
        }
        std::cerr << std::endl;
    }

    // Parse data lines
    size_t line_number = 2; // First data line is line 2 (after header)
    while (std::getline(input, line)) {
        std::istringstream ss(line);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }

        if (fields.size() < columnMap.size()) {
            line_number++;
            continue;
        }

        CSVFrame frame = {};
        frame.csv_line_number = line_number;

        // Helper to get field value
        auto getField = [&](const std::string& name, int defaultValue = 0) -> int {
            auto it = columnMap.find(name);
            if (it != columnMap.end() && it->second < fields.size()) {
                std::string val = fields[it->second];
                val.erase(0, val.find_first_not_of(" \t"));
                val.erase(val.find_last_not_of(" \t") + 1);
                try {
                    return std::stoi(val);
                } catch (...) {
                    return defaultValue;
                }
            }
            return defaultValue;
        };

        auto getFieldFloat = [&](const std::string& name, float defaultValue = 0.0f) -> float {
            auto it = columnMap.find(name);
            if (it != columnMap.end() && it->second < fields.size()) {
                std::string val = fields[it->second];
                val.erase(0, val.find_first_not_of(" \t"));
                val.erase(val.find_last_not_of(" \t") + 1);
                try {
                    return std::stof(val);
                } catch (...) {
                    return defaultValue;
                }
            }
            return defaultValue;
        };

        auto getFieldString = [&](const std::string& name) -> std::string {
            auto it = columnMap.find(name);
            if (it != columnMap.end() && it->second < fields.size()) {
                std::string val = fields[it->second];
                val.erase(0, val.find_first_not_of(" \t"));
                val.erase(val.find_last_not_of(" \t") + 1);
                return val;
            }
            return "";
        };

        // Extract data
        frame.time_us = getField("time (us)");

        // Quaternion - normalize from INAV scaled integers
        int q_raw[4] = {
            getField("quaternion[0]", 10000),
            getField("quaternion[1]", 0),
            getField("quaternion[2]", 0),
            getField("quaternion[3]", 0)
        };

        float q_norm = std::sqrt(q_raw[0]*q_raw[0] + q_raw[1]*q_raw[1] +
                                q_raw[2]*q_raw[2] + q_raw[3]*q_raw[3]);
        if (q_norm > 0) {
            frame.quaternion[0] = q_raw[0] / q_norm;  // w
            frame.quaternion[1] = q_raw[1] / q_norm;  // x
            frame.quaternion[2] = q_raw[2] / q_norm;  // y
            frame.quaternion[3] = q_raw[3] / q_norm;  // z
        } else {
            frame.quaternion[0] = 1.0f;  // Default identity quaternion
            frame.quaternion[1] = 0.0f;
            frame.quaternion[2] = 0.0f;
            frame.quaternion[3] = 0.0f;
        }

        // Navigation position
        frame.navPos[0] = getField("navPos[0]");
        frame.navPos[1] = getField("navPos[1]");
        frame.navPos[2] = getField("navPos[2]");



        // Flight mode flags
        frame.flightModeFlags = getFieldString("flightModeFlags (flags)");
        frame.mspOverrideFlags = getField("mspOverrideFlags", 0);
        frame.hasMSPRCOverride = hasMSPRCOverride(frame.flightModeFlags);

        // Attitude
        for (int i = 0; i < 3; i++) {
            frame.attitude[i] = getField("attitude[" + std::to_string(i) + "]");
            frame.gyroADC[i] = getField("gyroADC[" + std::to_string(i) + "]");
            frame.accSmooth[i] = getField("accSmooth[" + std::to_string(i) + "]");
        }

        // Servo outputs from rcData (AET format) - what INAV actually received as input
        // rcData values are already in 1000-2000 range from receiver/MSP override
        frame.servo[0] = getField("rcData[0]", 1500);  // Aileron (roll input)
        frame.servo[1] = getField("rcData[1]", 1500);  // Elevator (pitch input)
        frame.servo[2] = getField("rcData[3]", 1500);  // Throttle


        // Battery
        frame.vbat = getFieldFloat("vbat (V)");

        frames.push_back(frame);
        line_number++;
    }

    return frames;
}

std::vector<TestSpan> findTestSpans(const std::vector<CSVFrame>& frames) {
    std::vector<TestSpan> spans;

    bool in_test = false;
    size_t test_start = 0;

    for (size_t i = 0; i < frames.size(); i++) {
        if (!in_test && frames[i].hasMSPRCOverride) {
            // Start of test span
            in_test = true;
            test_start = i;
        } else if (in_test && !frames[i].hasMSPRCOverride) {
            // End of test span
            in_test = false;

            TestSpan span;
            span.start_index = test_start;
            span.end_index = i - 1;
            span.start_time_us = frames[test_start].time_us;
            span.end_time_us = frames[i - 1].time_us;
            span.frame_count = span.end_index - span.start_index + 1;

            spans.push_back(span);
        }
    }

    // Handle case where file ends during a test
    if (in_test) {
        TestSpan span;
        span.start_index = test_start;
        span.end_index = frames.size() - 1;
        span.start_time_us = frames[test_start].time_us;
        span.end_time_us = frames.back().time_us;
        span.frame_count = span.end_index - span.start_index + 1;

        spans.push_back(span);
    }

    return spans;
}

std::vector<CSVFrame> extractTestSpanWithBuffer(const std::vector<CSVFrame>& frames,
                                               const TestSpan& span,
                                               const std::vector<TestSpan>& allSpans) {

    // Calculate buffer start/end times
    uint32_t buffer_start_time = span.start_time_us - BUFFER_TIME_US;
    uint32_t buffer_end_time = span.end_time_us + BUFFER_TIME_US;

    // Find actual start/end indices with collision detection
    size_t actual_start = 0;
    size_t actual_end = frames.size() - 1;

    // Find start index (with buffer, but not into other tests)
    for (size_t i = 0; i < span.start_index; i++) {
        if (frames[i].time_us >= buffer_start_time) {
            // Check if this frame is in another test span
            bool in_other_test = false;
            for (const auto& other_span : allSpans) {
                if (&other_span != &span &&
                    i >= other_span.start_index && i <= other_span.end_index) {
                    in_other_test = true;
                    break;
                }
            }

            if (!in_other_test) {
                actual_start = i;
                break;
            }
        }
    }

    // Find end index (with buffer, but not into other tests)
    for (size_t i = frames.size() - 1; i > span.end_index; i--) {
        if (frames[i].time_us <= buffer_end_time) {
            // Check if this frame is in another test span
            bool in_other_test = false;
            for (const auto& other_span : allSpans) {
                if (&other_span != &span &&
                    i >= other_span.start_index && i <= other_span.end_index) {
                    in_other_test = true;
                    break;
                }
            }

            if (!in_other_test) {
                actual_end = i;
                break;
            }
        }
    }

    // Extract frames
    std::vector<CSVFrame> extracted;
    for (size_t i = actual_start; i <= actual_end; i++) {
        extracted.push_back(frames[i]);
    }

    std::cerr << "Extracted test span " << (&span - &allSpans[0] + 1)
              << ": " << extracted.size() << " frames" << std::endl;
    std::cerr << "  Test period: " << span.frame_count << " frames" << std::endl;
    std::cerr << "  Buffer before: " << (span.start_index - actual_start) << " frames" << std::endl;
    std::cerr << "  Buffer after: " << (actual_end - span.end_index) << " frames" << std::endl;

    return extracted;
}

void generateCppFile(const std::vector<CSVFrame>& frames, std::ostream& out) {
    // Get current time for header
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);

    out << "// Auto-generated flight data for XIAO-GP simulation\n";
    out << "// Generated on " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
    out << "// Frames: " << frames.size() << "\n\n";

    out << "#include <MSP.h>\n\n";

    out << "// Flight data frame structure\n";
    out << "struct FlightDataFrame {\n";
    out << "    uint32_t time_us;\n";
    out << "    msp_status_t status;\n";
    out << "    msp_attitude_quaternion_t attitude_quaternion;\n";
    out << "    msp_waypoint_t waypoint;\n";
    out << "    uint16_t servo[3];          // servo outputs: aileron, elevator, throttle (AET)\n";
    out << "};\n\n";

    out << "// Flight data array\n";
    out << "const FlightDataFrame flight_data[] = {\n";

    for (size_t i = 0; i < frames.size(); i++) {
        const auto& frame = frames[i];
        uint32_t flags = parseFlightModeFlags(frame.flightModeFlags);

        // Add MSPRCOVERRIDE flag when present in string flags
        if (hasMSPRCOverride(frame.flightModeFlags)) {
            flags |= (1 << 30);  // Set MSPRCOVERRIDE bit
        }

        out << "    { // Frame " << i << " (CSV line " << frame.csv_line_number << ")\n";
        out << "        .time_us = " << frame.time_us << ",\n";
        out << "        .status = {\n";
        out << "            .cycleTime = 4000,\n";
        out << "            .i2cErrorCounter = 0,\n";
        out << "            .sensor = 31,\n";  // All sensors present
        out << "            .flightModeFlags = " << flags << ",\n";
        out << "            .configProfileIndex = 0\n";
        out << "        },\n";

        out << "        .attitude_quaternion = {\n";
        out << "            .q = {" << std::fixed << std::setprecision(6)
            << frame.quaternion[0] << "f, " << frame.quaternion[1] << "f, "
            << frame.quaternion[2] << "f, " << frame.quaternion[3] << "f}\n";
        out << "        },\n";

        out << "        .waypoint = {\n";
        out << "            .waypointNumber = 255,\n";  // MSP_WP_CURRENT_POSITION
        out << "            .action = 1,\n";            // MSP_NAV_WP_ACTION_WAYPOINT
        out << "            .lat = " << frame.navPos[0] << ",\n";
        out << "            .lon = " << frame.navPos[1] << ",\n";
        out << "            .alt = " << frame.navPos[2] << ",\n";
        out << "            .p1 = 0,\n";
        out << "            .p2 = 0,\n";
        out << "            .p3 = 0,\n";
        out << "            .flag = 0\n";
        out << "        },\n";
        out << "        .servo = {" << frame.servo[0] << ", " << frame.servo[1] << ", " << frame.servo[2] << "}\n";
        out << "    }";

        if (i < frames.size() - 1) {
            out << ",";
        }
        out << "\n";
    }

    out << "};\n\n";

    out << "// Array size\n";
    out << "const size_t flight_data_count = " << frames.size() << ";\n\n";

    out << "// Access functions\n";
    out << "const FlightDataFrame* get_flight_data_frame(size_t index) {\n";
    out << "    if (index >= flight_data_count) {\n";
    out << "        return nullptr;\n";
    out << "    }\n";
    out << "    return &flight_data[index];\n";
    out << "}\n\n";

    out << "size_t get_flight_data_frame_count() {\n";
    out << "    return flight_data_count;\n";
    out << "}\n\n";

    out << "const FlightDataFrame* get_flight_data_frame_at_time_msec(uint32_t time_msec) {\n";
    out << "    uint32_t time_us = time_msec * 1000; // Convert msec to usec for comparison\n";
    out << "    for (size_t i = 0; i < flight_data_count; i++) {\n";
    out << "        if (flight_data[i].time_us >= time_us) {\n";
    out << "            return &flight_data[i];\n";
    out << "        }\n";
    out << "    }\n";
    out << "    return nullptr;\n";
    out << "}\n";

    std::cerr << "Generated " << frames.size() << " frames" << std::endl;
    std::cerr << "Memory usage: ~" << (frames.size() * 100) << " bytes" << std::endl;
}

void printUsage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [options]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  -i <file>    Input CSV file (use '-' for stdin)\n";
    std::cout << "  -t <num>     Test span number to extract (1-based)\n";
    std::cout << "  -l           List available test spans and exit\n";
    std::cout << "  -h           Show this help\n";
    std::cout << "\nExamples:\n";
    std::cout << "  " << prog_name << " -i flight9.csv -t 1 > flight_data.cpp\n";
    std::cout << "  cat flight9.csv | " << prog_name << " -i - -t 2 > flight_data.cpp\n";
    std::cout << "  " << prog_name << " -i flight9.csv -l\n";
}

Options parseCommandLine(int argc, char* argv[]) {
    Options opts;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-i" && i + 1 < argc) {
            opts.input_file = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            opts.test_span = std::stoi(argv[++i]);
        } else if (arg == "-l") {
            opts.list_spans = true;
        } else if (arg == "-h") {
            printUsage(argv[0]);
            exit(0);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            printUsage(argv[0]);
            exit(1);
        }
    }

    return opts;
}

int main(int argc, char* argv[]) {
    Options opts = parseCommandLine(argc, argv);

    // Open input stream
    std::istream* input = &std::cin;
    std::ifstream file_input;

    if (opts.input_file != "-") {
        file_input.open(opts.input_file);
        if (!file_input.is_open()) {
            std::cerr << "Error: Cannot open file " << opts.input_file << std::endl;
            return 1;
        }
        input = &file_input;
    }

    std::cerr << "Parsing CSV data..." << std::endl;
    auto frames = parseCSV(*input);

    if (frames.empty()) {
        std::cerr << "Error: No data frames parsed" << std::endl;
        return 1;
    }

    std::cerr << "Parsed " << frames.size() << " total frames" << std::endl;

    // Find test spans
    auto test_spans = findTestSpans(frames);

    if (test_spans.empty()) {
        std::cerr << "Error: No MSPRCOVERRIDE test spans found" << std::endl;
        return 1;
    }

    std::cerr << "Found " << test_spans.size() << " test spans:" << std::endl;
    for (size_t i = 0; i < test_spans.size(); i++) {
        const auto& span = test_spans[i];
        float duration_sec = (span.end_time_us - span.start_time_us) / 1000000.0f;

        std::cerr << "  Test " << (i + 1) << ": " << span.frame_count
                  << " frames, " << std::fixed << std::setprecision(1)
                  << duration_sec << "s duration" << std::endl;
    }

    if (opts.list_spans) {
        return 0; // Just list spans and exit
    }

    // Validate test span selection
    if (opts.test_span < 1 || opts.test_span > (int)test_spans.size()) {
        std::cerr << "Error: Test span " << opts.test_span
                  << " not found. Available spans: 1-" << test_spans.size() << std::endl;
        return 1;
    }

    // Extract selected test span with buffer
    const auto& selected_span = test_spans[opts.test_span - 1];
    auto extracted_frames = extractTestSpanWithBuffer(frames, selected_span, test_spans);

    if (extracted_frames.empty()) {
        std::cerr << "Error: No frames extracted for test span " << opts.test_span << std::endl;
        return 1;
    }

    // Generate C++ output
    std::cerr << "Generating C++ flight data..." << std::endl;
    generateCppFile(extracted_frames, std::cout);

    std::cerr << "Extraction complete!" << std::endl;
    return 0;
}