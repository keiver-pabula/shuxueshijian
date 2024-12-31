#include "json_config_handler.h"

JSONConfigHandler::JSONConfigHandler(const std::string& filePath) {
    std::ifstream configFile(filePath);
    if (!configFile.is_open()) {
        throw std::runtime_error("Unable to open configuration file: " + filePath);
    }

    // Read file line by line
    std::string line;
    while (std::getline(configFile, line)) {
        rawContent += line;
    }
    configFile.close();

    // Parse the JSON content manually
    parseConfig();
}

void JSONConfigHandler::parseConfig() {
    size_t pos = 0;

    // Parse the "parameters" section
    pos = rawContent.find("\"parameters\"");
    if (pos != std::string::npos) {
        size_t start = rawContent.find("{", pos);
        size_t end = rawContent.find("}", start);
        std::string parameters = rawContent.substr(start + 1, end - start - 1);

        size_t outputPos = parameters.find("\"output_directory\"");
        if (outputPos != std::string::npos) {
            size_t valueStart = parameters.find(":", outputPos) + 1;
            size_t valueEnd = parameters.find("\"", valueStart + 1);
            parsedData.parameters["output_directory"] = parameters.substr(valueStart + 1, valueEnd - valueStart - 1);
        }
    }

    // Parse the "curves" section
    pos = rawContent.find("\"curves\"");
    if (pos != std::string::npos) {
        size_t start = rawContent.find("[", pos);
        size_t end = rawContent.find("]", start);
        std::string curves = rawContent.substr(start + 1, end - start - 1);

        size_t curveStart = 0;
        while ((curveStart = curves.find("{", curveStart)) != std::string::npos) {
            size_t curveEnd = curves.find("}", curveStart);
            std::string curve = curves.substr(curveStart + 1, curveEnd - curveStart - 1);

            std::map<std::string, std::string> curveData;

            size_t namePos = curve.find("\"name\"");
            if (namePos != std::string::npos) {
                size_t valueStart = curve.find(":", namePos) + 1;
                size_t valueEnd = curve.find("\"", valueStart + 1);
                curveData["name"] = curve.substr(valueStart + 1, valueEnd - valueStart - 1);
            }

            size_t paramPos = curve.find("\"parameterization\"");
            if (paramPos != std::string::npos) {
                size_t valueStart = curve.find(":", paramPos) + 1;
                size_t valueEnd = curve.find("\"", valueStart + 1);
                curveData["parameterization"] = curve.substr(valueStart + 1, valueEnd - valueStart - 1);
            }

            size_t resPos = curve.find("\"resolution\"");
            if (resPos != std::string::npos) {
                size_t valueStart = curve.find(":", resPos) + 1;
                size_t valueEnd = curve.find(",", resPos);
                if (valueEnd == std::string::npos) {
                    valueEnd = curve.length();
                }
                curveData["resolution"] = curve.substr(valueStart, valueEnd - valueStart);
            }

            parsedData.curves.push_back(curveData);
            curveStart = curveEnd + 1;
        }
    }
}

const ParsedData& JSONConfigHandler::getConfig() const {
    return parsedData;
}
