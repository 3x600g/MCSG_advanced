#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <unordered_map>

struct XSCache {
    std::vector<std::pair<double, double>> captureData;
    std::vector<std::pair<double, double>> fissionData;
    std::vector<std::pair<double, double>> scatteringData;
};

class Material {
public:
    int m_number;
    double m_density;
    std::vector<std::pair<int, double>> composition;
    double totalMass = 0.0;
    std::vector<double> numberDensities;
    std::unordered_map<std::string, XSCache> xsCache;

    Material(int num, double dens);

    void addElement(int atomicNumber, double weightPercent);
    void calculateTotalMass();
    void calculateNumberDensities();
    void loadXSData(const std::string& base_path);
    void calculateXsAndStore(double incidentEnergy, double& cap, double& fis, double& scat);
};

std::vector<Material> readMaterialCards(const std::string& filePath);

#endif // MATERIAL_H
