#include "Calculation.h"
#include <fstream>
#include <sstream>

void Calculation::parseKcode(const std::string& line) {
    std::istringstream iss(line);
    int value;
    while (iss >> value) {
        kcode.push_back(value);
    }
}

void Calculation::parseKsrc(const std::string& line) {
    std::istringstream iss(line);
    float value;
    while (iss >> value) {
        ksrc.push_back(value);
    }
}

void readCalculationCards(const std::string& filePath, Calculation& calculation) {
    std::ifstream file(filePath);
    std::string line;

    while (getline(file, line)) {
        if (line.find("kcode") != std::string::npos) {
            calculation.parseKcode(line.substr(line.find("kcode") + 6));
        } else if (line.find("ksource") != std::string::npos) {
            calculation.parseKsrc(line.substr(line.find("ksource") + 7));
        }
    }
}
