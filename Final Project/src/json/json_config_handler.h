#ifndef JSON_CONFIG_HANDLER_H
#define JSON_CONFIG_HANDLER_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

struct ParsedData {
    std::map<std::string, std::string> parameters;
    std::vector<std::map<std::string, std::string>> curves;
};

class JSONConfigHandler {
public:
    // Constructor to load the configuration file
    explicit JSONConfigHandler(const std::string& filePath);

    // Method to retrieve parsed configuration
    const ParsedData& getConfig() const;

private:
    std::string rawContent; // Raw content of the JSON file
    ParsedData parsedData;  // Parsed JSON data

    void parseConfig(); // Method to manually parse the JSON content
};

#endif // JSON_CONFIG_HANDLER_H
