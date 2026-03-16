#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cmath>
#include <iomanip>

#include "nn_serialization.h"
#include "nn_evaluator_portable.h"

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " [OPTIONS]\n";
    std::cout << "Options:\n";
    std::cout << "  -i, --input FILE     Input NN weight file (required, NN01 format)\n";
    std::cout << "  -o, --output FILE    Output C++ source file (default: nn_program_generated.cpp)\n";
    std::cout << "  -f, --function NAME  Generated function name (default: generatedNNProgram)\n";
    std::cout << "  -u, --unrolled       Generate unrolled layer code (default: use nn_forward)\n";
    std::cout << "  --help               Show this help message\n";
    std::cout << "\n";
    std::cout << "Generates C++ source with embedded NN weights for desktop and embedded deployment.\n";
    std::cout << "The generated function has the same signature as generatedGPProgram().\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << progName << " -i nn_weights.dat\n";
    std::cout << "  " << progName << " -i nn_weights.dat -o my_nn.cpp -f myNNFunction\n";
    std::cout << "  " << progName << " -i nn_weights.dat -u   # unrolled layer loops\n";
}

// Generate code that uses the portable nn_forward() function with static weight array
std::string generatePortableCode(const NNGenome& genome, const std::string& functionName,
                                  const std::string& sourceFile) {
    std::stringstream code;

    code << "// Auto-generated NN evaluator function\n";
    code << "//\n";
    code << "// Source NN Information:\n";
    code << "//   Weight file: " << sourceFile << "\n";
    code << "//   Source:      " << genome.source << "\n";
    code << "//   Topology:    ";
    for (size_t i = 0; i < genome.topology.size(); i++) {
        if (i > 0) code << " -> ";
        code << genome.topology[i];
    }
    code << "\n";
    code << "//   Weights:     " << genome.weights.size() << "\n";
    code << "//   Fitness:     " << std::fixed << std::setprecision(6) << genome.fitness << "\n";
    code << "//   Sigma:       " << std::fixed << std::setprecision(6) << genome.mutation_sigma << "\n";
    code << "//\n";
    code << "#include \"nn_program.h\"\n\n";

    // Topology
    code << "static const int nn_topology[] = {";
    for (size_t i = 0; i < genome.topology.size(); i++) {
        if (i > 0) code << ", ";
        code << genome.topology[i];
    }
    code << "};\n";
    code << "static const int nn_num_layers = " << genome.topology.size() << ";\n\n";

    // Weights — formatted 8 per line for readability
    code << "static const float nn_weights[" << genome.weights.size() << "] = {\n";
    for (size_t i = 0; i < genome.weights.size(); i++) {
        if (i % 8 == 0) code << "    ";
        code << std::scientific << std::setprecision(8) << genome.weights[i] << "f";
        if (i < genome.weights.size() - 1) code << ",";
        if (i % 8 == 7 || i == genome.weights.size() - 1) code << "\n";
        else code << " ";
    }
    code << "};\n\n";

    // Topology as std::vector (needed by nn_forward)
    code << "static const std::vector<int>& getTopology() {\n";
    code << "    static const std::vector<int> t(nn_topology, nn_topology + nn_num_layers);\n";
    code << "    return t;\n";
    code << "}\n\n";

    // Main function
    code << "gp_scalar " << functionName << "(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {\n";
    code << "    float inputs[" << genome.topology.front() << "];\n";
    code << "    nn_gather_inputs(pathProvider, aircraftState, inputs);\n\n";
    code << "    float outputs[" << genome.topology.back() << "];\n";
    code << "    nn_forward(nn_weights, getTopology(), inputs, outputs);\n\n";
    code << "    // Set control commands: pitch, roll, throttle (already [-1,1] via tanh)\n";
    code << "    aircraftState.setPitchCommand(static_cast<gp_scalar>(outputs[0]));\n";
    code << "    aircraftState.setRollCommand(static_cast<gp_scalar>(outputs[1]));\n";
    code << "    aircraftState.setThrottleCommand(static_cast<gp_scalar>(outputs[2]));\n\n";
    code << "    return static_cast<gp_scalar>(outputs[0]); // return pitch for compatibility\n";
    code << "}\n\n";

    // Source identifier
    code << "const char* generatedNNProgramSource = \"" << genome.source << "\";\n";

    return code.str();
}

// Generate code with unrolled layer loops — no std::vector, no nn_forward call
std::string generateUnrolledCode(const NNGenome& genome, const std::string& functionName,
                                  const std::string& sourceFile) {
    std::stringstream code;

    code << "// Auto-generated NN evaluator function (unrolled)\n";
    code << "//\n";
    code << "// Source NN Information:\n";
    code << "//   Weight file: " << sourceFile << "\n";
    code << "//   Source:      " << genome.source << "\n";
    code << "//   Topology:    ";
    for (size_t i = 0; i < genome.topology.size(); i++) {
        if (i > 0) code << " -> ";
        code << genome.topology[i];
    }
    code << "\n";
    code << "//   Weights:     " << genome.weights.size() << "\n";
    code << "//   Fitness:     " << std::fixed << std::setprecision(6) << genome.fitness << "\n";
    code << "//   Sigma:       " << std::fixed << std::setprecision(6) << genome.mutation_sigma << "\n";
    code << "//\n";
    code << "#include \"nn_program.h\"\n\n";

    // Weights
    code << "static const float nn_weights[" << genome.weights.size() << "] = {\n";
    for (size_t i = 0; i < genome.weights.size(); i++) {
        if (i % 8 == 0) code << "    ";
        code << std::scientific << std::setprecision(8) << genome.weights[i] << "f";
        if (i < genome.weights.size() - 1) code << ",";
        if (i % 8 == 7 || i == genome.weights.size() - 1) code << "\n";
        else code << " ";
    }
    code << "};\n\n";

    // Main function with unrolled loops
    int max_layer = *std::max_element(genome.topology.begin(), genome.topology.end());
    code << "gp_scalar " << functionName << "(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {\n";
    code << "    float inputs[" << genome.topology.front() << "];\n";
    code << "    nn_gather_inputs(pathProvider, aircraftState, inputs);\n\n";

    // Layer buffers — fixed-size on stack
    code << "    float buf_a[" << max_layer << "], buf_b[" << max_layer << "];\n\n";

    // Copy inputs to buf_a
    code << "    // Copy inputs\n";
    code << "    for (int i = 0; i < " << genome.topology.front() << "; i++) buf_a[i] = inputs[i];\n\n";

    // Unrolled layer computation
    int weight_offset = 0;
    bool a_is_input = true;
    for (size_t layer = 0; layer + 1 < genome.topology.size(); layer++) {
        int in_size = genome.topology[layer];
        int out_size = genome.topology[layer + 1];
        std::string in_buf = a_is_input ? "buf_a" : "buf_b";
        std::string out_buf = a_is_input ? "buf_b" : "buf_a";

        code << "    // Layer " << layer << ": " << in_size << " -> " << out_size << "\n";
        code << "    {\n";
        code << "        const float* W = nn_weights + " << weight_offset << ";\n";
        code << "        const float* B = nn_weights + " << (weight_offset + in_size * out_size) << ";\n";
        code << "        for (int j = 0; j < " << out_size << "; j++) {\n";
        code << "            float sum = B[j];\n";
        code << "            for (int i = 0; i < " << in_size << "; i++) {\n";
        code << "                sum += W[j * " << in_size << " + i] * " << in_buf << "[i];\n";
        code << "            }\n";
        code << "            " << out_buf << "[j] = static_cast<float>(fast_tanh(static_cast<gp_scalar>(sum)));\n";
        code << "        }\n";
        code << "    }\n\n";

        weight_offset += in_size * out_size + out_size;
        a_is_input = !a_is_input;
    }

    // After the loop, a_is_input tracks which buffer has the final result.
    // a_is_input flips each layer: for odd layer count result is in buf_b, even in buf_a.
    std::string result_buf = a_is_input ? "buf_a" : "buf_b";

    code << "    // Set control commands from " << result_buf << "\n";
    code << "    aircraftState.setPitchCommand(static_cast<gp_scalar>(" << result_buf << "[0]));\n";
    code << "    aircraftState.setRollCommand(static_cast<gp_scalar>(" << result_buf << "[1]));\n";
    code << "    aircraftState.setThrottleCommand(static_cast<gp_scalar>(" << result_buf << "[2]));\n\n";
    code << "    return static_cast<gp_scalar>(" << result_buf << "[0]); // return pitch for compatibility\n";
    code << "}\n\n";

    // Source identifier
    code << "const char* generatedNNProgramSource = \"" << genome.source << "\";\n";

    return code.str();
}

int main(int argc, char** argv) {
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"function", required_argument, 0, 'f'},
        {"unrolled", no_argument, 0, 'u'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    std::string inputFile;
    std::string outputFile = "nn_program_generated.cpp";
    std::string functionName = "generatedNNProgram";
    bool unrolled = false;
    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "i:o:f:uh", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i': inputFile = optarg; break;
            case 'o': outputFile = optarg; break;
            case 'f': functionName = optarg; break;
            case 'u': unrolled = true; break;
            case 'h': printUsage(argv[0]); return 0;
            case '?': printUsage(argv[0]); return 1;
            default: break;
        }
    }

    if (inputFile.empty()) {
        std::cerr << "Error: Input file required (-i)" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Read NN weight file
    std::ifstream file(inputFile, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open input file: " << inputFile << std::endl;
        return 1;
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<uint8_t> data(size);
    if (!file.read(reinterpret_cast<char*>(data.data()), size)) {
        std::cerr << "Error reading input file" << std::endl;
        return 1;
    }
    file.close();

    // Verify NN01 format
    if (!nn_detect_format(data.data(), data.size())) {
        std::cerr << "Error: Input file is not in NN01 format" << std::endl;
        return 1;
    }

    // Deserialize
    NNGenome genome;
    if (!nn_deserialize(data.data(), data.size(), genome)) {
        std::cerr << "Error deserializing NN genome" << std::endl;
        return 1;
    }

    // Generate code
    std::string code;
    if (unrolled) {
        code = generateUnrolledCode(genome, functionName, inputFile);
    } else {
        code = generatePortableCode(genome, functionName, inputFile);
    }

    // Write output
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
        return 1;
    }
    outFile << code;
    outFile.close();

    std::cout << "Generated " << (unrolled ? "unrolled" : "portable") << " NN evaluator: " << outputFile << std::endl;
    std::cout << "  Topology:   ";
    for (size_t i = 0; i < genome.topology.size(); i++) {
        if (i > 0) std::cout << " -> ";
        std::cout << genome.topology[i];
    }
    std::cout << std::endl;
    std::cout << "  Weights:    " << genome.weights.size() << std::endl;
    std::cout << "  Function:   " << functionName << "()" << std::endl;
    std::cout << "  Fitness:    " << std::fixed << std::setprecision(6) << genome.fitness << std::endl;

    return 0;
}
