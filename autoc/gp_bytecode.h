#ifndef GP_BYTECODE_H
#define GP_BYTECODE_H

#include <vector>
#include <string>
#include <cstdint>
#include "aircraft_state.h"

// Bytecode instruction structure
struct GPBytecode {
    uint8_t opcode;     // Operation code (maps to Operators enum)
    uint8_t argc;       // Number of arguments
    float constant;     // For literal values (PI, 0, 1, 2)
    
    GPBytecode(uint8_t op = 0, uint8_t args = 0, float val = 0.0f) 
        : opcode(op), argc(args), constant(val) {}
};

// Bytecode file header
struct GPBytecodeHeader {
    uint32_t magic;         // File format magic number
    uint32_t version;       // Bytecode format version
    uint32_t instruction_count;
    uint32_t fitness_int;   // Fitness as fixed-point (fitness * 1000000)
    uint32_t length;        // Original GP length
    uint32_t depth;         // Original GP depth
    char s3_key[256];       // Source S3 key
    uint32_t generation;    // Generation number
    
    static const uint32_t MAGIC = 0x47504243; // "GPBC"
    static const uint32_t VERSION = 1;
};

// Stack-based bytecode interpreter for autoc
class GPBytecodeInterpreter {
private:
    std::vector<GPBytecode> program;
    GPBytecodeHeader header;
    static const int MAX_STACK_SIZE = 256;
    
    // Helper functions for complex operations
    int getIndex(const std::vector<Path>& path, double arg, AircraftState& aircraftState);
    double executeGetDPhi(const std::vector<Path>& path, double arg, AircraftState& aircraftState);
    double executeGetDTheta(const std::vector<Path>& path, double arg, AircraftState& aircraftState);
    double executeGetDTarget(const std::vector<Path>& path, double arg, AircraftState& aircraftState);
    double executeGetDHome(AircraftState& aircraftState);
    
public:
    GPBytecodeInterpreter();
    ~GPBytecodeInterpreter();
    
    // Load bytecode program from file
    bool loadProgram(const std::string& filename);
    
    // Get program information
    const GPBytecodeHeader& getHeader() const { return header; }
    double getFitness() const { return header.fitness_int / 1000000.0; }
    uint32_t getLength() const { return header.length; }
    uint32_t getDepth() const { return header.depth; }
    const std::string getS3Key() const { return std::string(header.s3_key); }
    uint32_t getGeneration() const { return header.generation; }
    
    // Main evaluation function - compatible with autoc's existing interface
    double evaluate(AircraftState& aircraftState, std::vector<Path>& path, double arg);
    
    // Debug functions
    void printProgram() const;
    bool isLoaded() const { return !program.empty(); }
};

#endif // GP_BYTECODE_H