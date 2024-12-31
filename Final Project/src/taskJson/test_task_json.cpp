#include <iostream>
#include <filesystem>
#include "../json/json_config_handler.h"

int main() {
    try {
        // Load the configuration file
        JSONConfigHandler jsonHandler("config/config.json");
        const ParsedData& config = jsonHandler.getConfig();

        // Create the output directory
        std::filesystem::create_directories("Output/taskJson");

        // Output the parsed data
        std::ofstream outputFile("Output/taskJson/config_output.txt");
        if (!outputFile.is_open()) {
            throw std::runtime_error("Failed to open output file for writing.");
        }

        // Write parameters
        outputFile << "Parameters:\n";
        for (const auto& [key, value] : config.parameters) {
            outputFile << key << ": " << value << "\n";
        }

        // Write curves
        outputFile << "\nCurves:\n";
        for (const auto& curve : config.curves) {
            outputFile << "Curve:\n";
            for (const auto& [key, value] : curve) {
                outputFile << "  " << key << ": " << value << "\n";
            }
        }

        outputFile.close();
        std::cout << "Configuration processed and saved to 'Output/taskJson/config_output.txt'.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
