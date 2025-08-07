#include "gp_bytecode.h"
#include "gp_evaluator_desktop.h"
#include "autoc.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Geometry>

GPBytecodeInterpreter::GPBytecodeInterpreter() {
    program.clear();
    memset(&header, 0, sizeof(header));
}

GPBytecodeInterpreter::~GPBytecodeInterpreter() {
    program.clear();
}

bool GPBytecodeInterpreter::loadProgram(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open bytecode file: " << filename << std::endl;
        return false;
    }
    
    // Read header
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (file.gcount() != sizeof(header)) {
        std::cerr << "Error: Invalid bytecode file header" << std::endl;
        return false;
    }
    
    // Verify magic number and version
    if (header.magic != GPBytecodeHeader::MAGIC) {
        std::cerr << "Error: Invalid bytecode file magic number" << std::endl;
        return false;
    }
    
    if (header.version != GPBytecodeHeader::VERSION) {
        std::cerr << "Error: Unsupported bytecode version: " << header.version << std::endl;
        return false;
    }
    
    // Read instructions
    program.clear();
    program.resize(header.instruction_count);
    
    for (uint32_t i = 0; i < header.instruction_count; i++) {
        file.read(reinterpret_cast<char*>(&program[i]), sizeof(GPBytecode));
        if (file.gcount() != sizeof(GPBytecode)) {
            std::cerr << "Error: Failed to read instruction " << i << std::endl;
            program.clear();
            return false;
        }
    }
    
    file.close();
    
    std::cout << "Loaded GP bytecode program:" << std::endl;
    std::cout << "  S3 Key: " << getS3Key() << std::endl;
    std::cout << "  Generation: " << getGeneration() << std::endl;
    std::cout << "  Length: " << getLength() << ", Depth: " << getDepth() << std::endl;
    std::cout << "  Fitness: " << getFitness() << std::endl;
    std::cout << "  Instructions: " << header.instruction_count << std::endl;
    
    return true;
}

// Helper functions are now in gp_evaluator_portable.cc

double GPBytecodeInterpreter::evaluate(AircraftState& aircraftState, std::vector<Path>& path, double arg) {
    if (program.empty()) {
        std::cerr << "Error: No bytecode program loaded" << std::endl;
        return 0.0;
    }
    
    // Delegate to desktop extension which handles format conversion
    return GPEvaluatorDesktop::evaluateGPBytecode(program, path, aircraftState, arg);
}

void GPBytecodeInterpreter::printProgram() const {
    std::cout << "GP Bytecode Program (" << program.size() << " instructions):" << std::endl;
    for (size_t i = 0; i < program.size(); i++) {
        const auto& inst = program[i];
        std::cout << "  " << i << ": op=" << (int)inst.opcode << " argc=" << (int)inst.argc;
        if (inst.argc == 0) {
            std::cout << " const=" << inst.constant;
        }
        std::cout << std::endl;
    }
}