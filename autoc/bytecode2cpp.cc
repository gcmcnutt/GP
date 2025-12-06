#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cmath>
#include <algorithm>

#include "gp_bytecode.h"
#include "gp_evaluator_portable.h"
#include "gp_evaluator_desktop.h"
#include "autoc.h"

class BytecodeToSourceGenerator {
private:
    const struct GPBytecode* program;  // Use GP build GPBytecode
    int program_size;
    GPBytecodeHeader header;  // Store header information for comments
    
    int analyzeStackDepth(const struct GPBytecode* program, int size) {
        int maxStack = 0;
        int currentStack = 0;
        
        for (int i = 0; i < size; i++) {
            int opcode = program[i].opcode;
            int argc = program[i].argc;
            
            // Determine stack effect for each operation
            switch (opcode) {
                // Zero-arg operations (terminals) - push 1
                case OP_PI: case ZERO: case ONE: case TWO:
                case GETPITCH: case GETROLL: case GETTHROTTLE: case GETVEL:
                case GETDHOME: case GETALPHA: case GETBETA:
                case GETVELX: case GETVELY: case GETVELZ:
                case GETROLL_RAD: case GETPITCH_RAD:
                    currentStack += 1;
                    break;
                    
                // Unary operations - pop 1, push 1 (net 0)
                case SIN: case COS: case ABS: case SQRT:
                case SETPITCH: case SETROLL: case SETTHROTTLE:
                case GETDPHI: case GETDTHETA: case GETDTARGET:
                    // net 0
                    break;
                    
                // Binary operations - pop 2, push 1 (net -1)
                case ADD: case SUB: case MUL: case DIV:
                case EQ: case GT: case ATAN2: case MIN: case MAX:
                    currentStack -= 1;
                    break;
                    
                // Ternary operations - pop 3, push 1 (net -2)
                case IF: case CLAMP:
                    currentStack -= 2;
                    break;
                    
                // PROGN - pop 2, push 1 (net -1)
                case PROGN:
                    currentStack -= 1;
                    break;
                    
                default:
                    // Unknown operation, assume no net effect
                    break;
            }
            
            if (currentStack > maxStack) {
                maxStack = currentStack;
            }
        }
        
        // Use same max stack size as desktop evaluator for consistency
        // Add safety margin but cap at reasonable maximum
        int calculatedMax = maxStack + 8;
        return std::min(calculatedMax, 256);  // Same as desktop evaluator
    }
    
    std::string getOperatorName(int opcode) {
        switch (opcode) {
            case ADD: return "ADD";
            case SUB: return "SUB";
            case MUL: return "MUL";
            case DIV: return "DIV";
            case IF: return "IF";
            case EQ: return "EQ";
            case GT: return "GT";
            case SIN: return "SIN";
            case COS: return "COS";
            case SETPITCH: return "SETPITCH";
            case SETROLL: return "SETROLL";
            case SETTHROTTLE: return "SETTHROTTLE";
            case GETPITCH: return "GETPITCH";
            case GETROLL: return "GETROLL";
            case GETTHROTTLE: return "GETTHROTTLE";
            case GETVEL: return "GETVEL";
            case GETDPHI: return "GETDPHI";
            case GETDTHETA: return "GETDTHETA";
            case GETDTARGET: return "GETDTARGET";
            case GETDHOME: return "GETDHOME";
            case GETALPHA: return "GETALPHA";
            case GETBETA: return "GETBETA";
            case GETVELX: return "GETVELX";
            case GETVELY: return "GETVELY";
            case GETVELZ: return "GETVELZ";
            case GETROLL_RAD: return "GETROLL_RAD";
            case GETPITCH_RAD: return "GETPITCH_RAD";
            case CLAMP: return "CLAMP";
            case ATAN2: return "ATAN2";
            case ABS: return "ABS";
            case SQRT: return "SQRT";
            case MIN: return "MIN";
            case MAX: return "MAX";
            case OP_PI: return "OP_PI";
            case ZERO: return "ZERO";
            case ONE: return "ONE";
            case TWO: return "TWO";
            case PROGN: return "PROGN";
            default: return "UNKNOWN";
        }
    }
    
    void generateInstruction(std::stringstream& code, const struct GPBytecode& inst, int pc) {
        std::string opName = getOperatorName(inst.opcode);
        int opcodeInt = (int)inst.opcode;  // Cast to int for safe printing
        
        switch (inst.opcode) {
            // Binary operations
            case ADD: case SUB: case MUL: case DIV:
            case EQ: case GT: case ATAN2: case MIN: case MAX:
                code << "    // " << opName << "\n";
                code << "    {\n";
                code << "        fastmath::GPScalar args[2] = {stack[sp-2], stack[sp-1]};\n";
                code << "        sp -= 2;\n";
                code << "        stack[sp++] = evaluateGPOperator(" << opcodeInt << ", pathProvider, aircraftState, args, 2, context);\n";
                code << "    }\n";
                break;
                
            // Unary operations
            case SIN: case COS: case ABS: case SQRT:
            case SETPITCH: case SETROLL: case SETTHROTTLE:
            case GETDPHI: case GETDTHETA: case GETDTARGET:
                code << "    // " << opName << "\n";
                code << "    {\n";
                code << "        fastmath::GPScalar args[1] = {stack[sp-1]};\n";
                code << "        sp -= 1;\n";
                code << "        stack[sp++] = evaluateGPOperator(" << opcodeInt << ", pathProvider, aircraftState, args, 1, context);\n";
                code << "    }\n";
                break;
                
            // Ternary operations
            case IF: case CLAMP:
                code << "    // " << opName << "\n";
                code << "    {\n";
                code << "        fastmath::GPScalar args[3] = {stack[sp-3], stack[sp-2], stack[sp-1]};\n";
                code << "        sp -= 3;\n";
                code << "        stack[sp++] = evaluateGPOperator(" << opcodeInt << ", pathProvider, aircraftState, args, 3, context);\n";
                code << "    }\n";
                break;
                
            // PROGN - special case
            case PROGN:
                code << "    // PROGN\n";
                code << "    {\n";
                code << "        fastmath::GPScalar args[2] = {stack[sp-2], stack[sp-1]};\n";
                code << "        sp -= 2;\n";
                code << "        stack[sp++] = evaluateGPOperator(" << opcodeInt << ", pathProvider, aircraftState, args, 2, context);\n";
                code << "    }\n";
                break;
                
            // Zero-argument operations (terminals/sensors)
            default:
                code << "    stack[sp++] = evaluateGPOperator(" << opcodeInt << ", pathProvider, aircraftState, nullptr, 0, context); // " << opName << "\n";
                break;
        }
    }
    
public:
    BytecodeToSourceGenerator(const struct GPBytecode* prog, int size, const GPBytecodeHeader& hdr) 
        : program(prog), program_size(size), header(hdr) {}
    
    std::string generateEvaluatorFunction(const std::string& functionName) {
        std::stringstream code;
        
        code << "// Auto-generated GP evaluator function\n";
        code << "//\n";
        code << "// Source GP Information:\n";
        code << "//   S3 Key: " << std::string(header.s3_key) << "\n";
        code << "//   Generation: " << header.generation << "\n";
        code << "//   Original Length: " << header.length << "\n";
        code << "//   Original Depth: " << header.depth << "\n";
        code << "//   Fitness: " << (header.fitness_int / 1000000.0) << "\n";
        code << "//   Bytecode Instructions: " << program_size << "\n";
        code << "//\n";
        code << "#include \"gp_program.h\"\n";
        code << "#include \"fastmath/gp_scalar.h\"\n\n";
        
        code << "double " << functionName << "(PathProvider& pathProvider, AircraftState& aircraftState, double arg) {\n";
        code << "    fastmath::GPScalar context = fastmath::GPScalar::fromDouble(arg);\n";
        
        // Analyze bytecode to determine max stack depth
        int maxStack = analyzeStackDepth(program, program_size);
        code << "    fastmath::GPScalar stack[" << maxStack << "];\n";
        code << "    int sp = 0;\n\n";
        
        // Generate unrolled execution
        for (int i = 0; i < program_size; i++) {
            generateInstruction(code, program[i], i);
        }
        
        code << "\n    return applyRangeLimit(stack[0]).toDouble();\n";
        code << "}\n";
        
        return code.str();
    }
    
    // Header generation removed - use stable gp_program.h instead
    
    std::string generateStaticBytecodeArray(const std::string& functionName) {
        std::stringstream code;
        
        code << "// Auto-generated GP bytecode array for embedded evaluation\n";
        code << "//\n";
        code << "// Source GP Information:\n";
        code << "//   S3 Key: " << std::string(header.s3_key) << "\n";
        code << "//   Generation: " << header.generation << "\n";
        code << "//   Original Length: " << header.length << "\n";
        code << "//   Original Depth: " << header.depth << "\n";
        code << "//   Fitness: " << (header.fitness_int / 1000000.0) << "\n";
        code << "//   Bytecode Instructions: " << program_size << "\n";
        code << "//\n";
        code << "#include \"gp_program.h\"\n\n";
        
        // Generate static bytecode array
        code << "// Static bytecode program for embedded evaluation\n";
        code << "static const GPBytecode embedded_gp_bytecode[" << program_size << "] = {\n";
        for (int i = 0; i < program_size; i++) {
            const struct GPBytecode& inst = program[i];
            code << "    {" << (int)inst.opcode << ", " << (int)inst.argc << ", ";
            if (inst.constant == 0.0f) {
                code << "0.0f";
            } else {
                code << inst.constant << "f";
            }
            code << "}";
            if (i < program_size - 1) code << ",";
            code << "  // " << getOperatorName(inst.opcode);
            if (inst.argc > 0) code << " (argc=" << (int)inst.argc << ")";
            if (inst.constant != 0.0f) code << " (const=" << inst.constant << ")";
            code << "\n";
        }
        code << "};\n\n";
        
        code << "static const int embedded_gp_bytecode_size = " << program_size << ";\n\n";
        
        // Generate evaluation function that uses bytecode interpreter
        code << "double " << functionName << "(PathProvider& pathProvider, AircraftState& aircraftState, double arg) {\n";
        code << "    // Use portable bytecode evaluator for consistent behavior\n";
        code << "    return evaluateBytecodePortable(embedded_gp_bytecode, embedded_gp_bytecode_size, pathProvider, aircraftState, arg);\n";
        code << "}\n";
        
        return code.str();
    }
    
};

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " [OPTIONS]\n";
    std::cout << "Options:\n";
    std::cout << "  -i, --input FILE     Input bytecode file (required)\n";
    std::cout << "  -o, --output FILE    Output C++ source file (default: gp_program_generated.cpp)\n";
    std::cout << "  -f, --function NAME  Generated function name (default: generatedGPProgram)\n";
    std::cout << "  -b, --bytecode       Generate static bytecode array (embedded-friendly)\n";
    std::cout << "  --help               Show this help message\n";
    std::cout << "\n";
    std::cout << "Note: Function declaration is provided by stable gp_program.h header\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << progName << " -i gp_program.dat\n";
    std::cout << "  " << progName << " -i gp_program.dat -o my_gp.cpp -f myGPFunction\n";
    std::cout << "  " << progName << " -b -i gp_program.dat\n";
}

int main(int argc, char** argv) {
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"function", required_argument, 0, 'f'},
        {"bytecode", no_argument, 0, 'b'},
        {"help", no_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    std::string inputFile = "";
    std::string outputFile = "gp_program_generated.cpp";
    std::string functionName = "generatedGPProgram";
    bool generateBytecode = false;
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "i:o:f:b", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'f':
                functionName = optarg;
                break;
            case 'b':
                generateBytecode = true;
                break;
            case 0:
                if (long_options[option_index].name == std::string("help")) {
                    printUsage(argv[0]);
                    return 0;
                }
                break;
            case '?':
                printUsage(argv[0]);
                return 1;
            default:
                break;
        }
    }
    
    if (inputFile.empty()) {
        std::cerr << "Error: Input bytecode file is required\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Load bytecode program
    GPBytecodeInterpreter interpreter;
    if (!interpreter.loadProgram(inputFile)) {
        std::cerr << "Error: Failed to load bytecode from " << inputFile << std::endl;
        return 1;
    }
    
    // Access the program data (we need to expose this in the interpreter)
    // For now, we'll create a simple loader
    std::ifstream file(inputFile, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open bytecode file: " << inputFile << std::endl;
        return 1;
    }
    
    // Read header
    GPBytecodeHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    
    // Read instructions (GP format)
    std::vector<struct GPBytecode> gpProgram(header.instruction_count);
    for (uint32_t i = 0; i < header.instruction_count; i++) {
        file.read(reinterpret_cast<char*>(&gpProgram[i]), sizeof(struct GPBytecode));
    }
    file.close();
    
    // In GP build, no conversion needed - use directly as program
    std::vector<struct GPBytecode>& program = gpProgram;
    
    // Generate source code
    BytecodeToSourceGenerator generator(program.data(), program.size(), header);
    
    // Write C++ source file
    std::ofstream cppFile(outputFile);
    if (!cppFile.is_open()) {
        std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
        return 1;
    }
    if (generateBytecode) {
        cppFile << generator.generateStaticBytecodeArray(functionName);
    } else {
        cppFile << generator.generateEvaluatorFunction(functionName);
    }
    cppFile.close();
    
    
    std::cout << "Generated C++ source: " << outputFile << std::endl;
    std::cout << "Function name: " << functionName << std::endl;
    std::cout << "Instructions: " << program.size() << std::endl;
    std::cout << "Note: Function declaration available in gp_program.h" << std::endl;
    
    return 0;
}
