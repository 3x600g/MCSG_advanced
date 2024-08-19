#ifndef CALCULATION_H
#define CALCULATION_H

#include <vector>
#include <string>

class Calculation {
public:
    std::vector<int> kcode; // kcode 정보를 저장할 벡터 (예: [3000, 60, 150])
    std::vector<float> ksrc; // ksrc 정보를 저장할 벡터 (예: [0, 0, 0])

    void parseKcode(const std::string& line);
    void parseKsrc(const std::string& line);
};

void readCalculationCards(const std::string& filePath, Calculation& calculation);

#endif // CALCULATION_H
