#include "gp_bytecode.h"
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

int GPBytecodeInterpreter::getIndex(const std::vector<Path>& path, double arg, AircraftState& aircraftState) {
    if (std::isnan(arg)) {
        return aircraftState.getThisPathIndex();
    }
    int steps = std::max(-5, std::min(5, (int)arg));
    double distanceSoFar = path.at(aircraftState.getThisPathIndex()).distanceFromStart;
    double distanceGoal = distanceSoFar + steps * 22.0 * (200.0 / 1000.0);
    int currentStep = aircraftState.getThisPathIndex();
    
    if (steps > 0) {
        while (path.at(currentStep).distanceFromStart < distanceGoal && currentStep < (int)path.size() - 1) {
            currentStep++;
        }
    } else {
        while (path.at(currentStep).distanceFromStart > distanceGoal && currentStep > 0) {
            currentStep--;
        }
    }
    return currentStep;
}

double GPBytecodeInterpreter::executeGetDPhi(const std::vector<Path>& path, double arg, AircraftState& aircraftState) {
    int index = getIndex(path, arg, aircraftState);
    Eigen::Vector3d craftToTarget = path.at(index).start - aircraftState.getPosition();
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    return std::atan2(projectedVector.y(), -projectedVector.z());
}

double GPBytecodeInterpreter::executeGetDTheta(const std::vector<Path>& path, double arg, AircraftState& aircraftState) {
    int index = getIndex(path, arg, aircraftState);
    Eigen::Vector3d craftToTarget = path.at(index).start - aircraftState.getPosition();
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    double rollEstimate = std::atan2(projectedVector.y(), -projectedVector.z());
    Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
    Eigen::Quaterniond virtualOrientation = aircraftState.getOrientation() * rollRotation;
    Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;
    return std::atan2(-newLocalTargetVector.z(), newLocalTargetVector.x());
}

double GPBytecodeInterpreter::executeGetDTarget(const std::vector<Path>& path, double arg, AircraftState& aircraftState) {
    int index = getIndex(path, arg, aircraftState);
    double distance = (path.at(index).start - aircraftState.getPosition()).norm();
    return std::max(-1.0, std::min(1.0, (distance - 10.0) / aircraftState.getRelVel()));
}

double GPBytecodeInterpreter::executeGetDHome(AircraftState& aircraftState) {
    return (Eigen::Vector3d(0, 0, -10.0) - aircraftState.getPosition()).norm();
}

double GPBytecodeInterpreter::evaluate(AircraftState& aircraftState, std::vector<Path>& path, double arg) {
    if (program.empty()) {
        std::cerr << "Error: No bytecode program loaded" << std::endl;
        return 0.0;
    }
    
    float stack[MAX_STACK_SIZE];
    int stack_ptr = 0;
    
    // Execute bytecode instructions
    for (const auto& instruction : program) {
        switch (instruction.opcode) {
            case ADD: {
                if (stack_ptr < 2) return 0.0;
                float b = stack[--stack_ptr];
                float a = stack[--stack_ptr];
                stack[stack_ptr++] = a + b;
                break;
            }
            case SUB: {
                if (stack_ptr < 2) return 0.0;
                float b = stack[--stack_ptr];
                float a = stack[--stack_ptr];
                stack[stack_ptr++] = -a - b;  // Note: autoc uses -a - b
                break;
            }
            case MUL: {
                if (stack_ptr < 2) return 0.0;
                float b = stack[--stack_ptr];
                float a = stack[--stack_ptr];
                stack[stack_ptr++] = a * b;
                break;
            }
            case DIV: {
                if (stack_ptr < 2) return 0.0;
                float b = stack[--stack_ptr];
                float a = stack[--stack_ptr];
                stack[stack_ptr++] = (b == 0.0f) ? 0.0f : (a / b);
                break;
            }
            case IF: {
                if (stack_ptr < 3) return 0.0;
                float falseVal = stack[--stack_ptr];
                float trueVal = stack[--stack_ptr];
                float condition = stack[--stack_ptr];
                stack[stack_ptr++] = (condition != 0.0f) ? trueVal : falseVal;
                break;
            }
            case EQ: {
                if (stack_ptr < 2) return 0.0;
                float b = stack[--stack_ptr];
                float a = stack[--stack_ptr];
                stack[stack_ptr++] = (a == b) ? 1.0f : 0.0f;
                break;
            }
            case GT: {
                if (stack_ptr < 2) return 0.0;
                float b = stack[--stack_ptr];
                float a = stack[--stack_ptr];
                stack[stack_ptr++] = (a > b) ? 1.0f : 0.0f;
                break;
            }
            case SIN: {
                if (stack_ptr < 1) return 0.0;
                stack[stack_ptr-1] = std::sin(stack[stack_ptr-1]);
                break;
            }
            case COS: {
                if (stack_ptr < 1) return 0.0;
                stack[stack_ptr-1] = std::cos(stack[stack_ptr-1]);
                break;
            }
            case SETPITCH: {
                if (stack_ptr < 1) return 0.0;
                float value = stack[--stack_ptr];
                stack[stack_ptr++] = aircraftState.setPitchCommand(value);
                break;
            }
            case SETROLL: {
                if (stack_ptr < 1) return 0.0;
                float value = stack[--stack_ptr];
                stack[stack_ptr++] = aircraftState.setRollCommand(value);
                break;
            }
            case SETTHROTTLE: {
                if (stack_ptr < 1) return 0.0;
                float value = stack[--stack_ptr];
                stack[stack_ptr++] = aircraftState.setThrottleCommand(value);
                break;
            }
            case GETPITCH: {
                stack[stack_ptr++] = aircraftState.getPitchCommand();
                break;
            }
            case GETROLL: {
                stack[stack_ptr++] = aircraftState.getRollCommand();
                break;
            }
            case GETTHROTTLE: {
                stack[stack_ptr++] = aircraftState.getThrottleCommand();
                break;
            }
            case GETVEL: {
                stack[stack_ptr++] = aircraftState.getRelVel();
                break;
            }
            case GETDPHI: {
                if (stack_ptr < 1) return 0.0;
                float argVal = stack[--stack_ptr];
                stack[stack_ptr++] = executeGetDPhi(path, argVal, aircraftState);
                break;
            }
            case GETDTHETA: {
                if (stack_ptr < 1) return 0.0;
                float argVal = stack[--stack_ptr];
                stack[stack_ptr++] = executeGetDTheta(path, argVal, aircraftState);
                break;
            }
            case GETDTARGET: {
                if (stack_ptr < 1) return 0.0;
                float argVal = stack[--stack_ptr];
                stack[stack_ptr++] = executeGetDTarget(path, argVal, aircraftState);
                break;
            }
            case GETDHOME: {
                stack[stack_ptr++] = executeGetDHome(aircraftState);
                break;
            }
            case PI: {
                stack[stack_ptr++] = M_PI;
                break;
            }
            case ZERO: {
                stack[stack_ptr++] = 0.0f;
                break;
            }
            case ONE: {
                stack[stack_ptr++] = 1.0f;
                break;
            }
            case TWO: {
                stack[stack_ptr++] = 2.0f;
                break;
            }
            case PROGN: {
                if (stack_ptr < 2) return 0.0;
                stack[--stack_ptr];  // Discard first value
                // Second value remains on stack as result
                break;
            }
            default:
                std::cerr << "Error: Unknown bytecode instruction: " << (int)instruction.opcode << std::endl;
                return 0.0;
        }
        
        // Stack overflow protection
        if (stack_ptr >= MAX_STACK_SIZE) {
            std::cerr << "Error: Stack overflow in bytecode execution" << std::endl;
            return 0.0;
        }
    }
    
    // Return final result
    if (stack_ptr != 1) {
        std::cerr << "Error: Invalid stack state after bytecode execution (stack_ptr=" << stack_ptr << ")" << std::endl;
        return 0.0;
    }
    
    double returnValue = stack[0];
    
    // Apply range limiting (same as native GP)
    const double RANGELIMIT = 1000000.0;
    if (returnValue < -RANGELIMIT) return -RANGELIMIT;
    if (returnValue > RANGELIMIT) return RANGELIMIT;
    if (std::abs(returnValue) < 0.000001) return 0.0;
    
    return returnValue;
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