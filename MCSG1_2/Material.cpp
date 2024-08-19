#include "Material.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

Material::Material(int num, double dens) : m_number(num), m_density(dens) {}

void Material::addElement(int atomicNumber, double weightPercent) {
    composition.push_back({atomicNumber, weightPercent});
}

void Material::calculateTotalMass() {
    for (const auto& element : composition) {
        int massNumber = element.first % 1000;
        double weight = element.second;
        totalMass += massNumber * (weight / 100.0);
    }
}

void Material::calculateNumberDensities() {
    numberDensities.clear();
    const double NA = 6.02e23;
    for (const auto& element : composition) {
        int massNumber = element.first % 1000;
        double weightPercent = element.second;
        double nd = (weightPercent / 100.0) * m_density * NA / totalMass * 1e-24;
        numberDensities.push_back(nd);
    }
}

void Material::loadXSData(const std::string& base_path) {
    for (const auto& comp : composition) {
        std::string isotop = std::to_string(comp.first);
        XSCache cache;
        std::vector<std::string> reactions = {"capture", "fission", "scattering"};
        for (const auto& reaction : reactions) {
            std::string filename = base_path + isotop + "_" + reaction + ".txt";
            std::ifstream file(filename);
            double energy, xs;
            if (file.is_open()) {
                while (file >> energy >> xs) {
                    if (reaction == "capture") {
                        cache.captureData.emplace_back(energy, xs);
                    } else if (reaction == "fission") {
                        cache.fissionData.emplace_back(energy, xs);
                    } else if (reaction == "scattering") {
                        cache.scatteringData.emplace_back(energy, xs);
                    }
                }
            }
        }
        xsCache[isotop] = cache;
    }
}

void Material::calculateXsAndStore(double incidentEnergy, double& cap, double& fis, double& scat) {
    cap = fis = scat = 0.0;
    for (const auto& comp : composition) {
        std::string isotop = std::to_string(comp.first);
        double weightPercent = comp.second;
        double numberDensity = (weightPercent / 100.0) * m_density * 6.02e+23 / totalMass * 1e-24;
        const auto& cache = xsCache[isotop];

        auto getCachedXSValue = [](const std::vector<std::pair<double, double>>& xsData, double incidentEnergy) -> double {
            if (xsData.empty()) return 0.0;
            auto it = std::lower_bound(xsData.begin(), xsData.end(), incidentEnergy,
                                       [](const std::pair<double, double>& data, double energy) {
                                           return data.first < energy;
                                       });
            if (it == xsData.end()) return xsData.back().second;
            if (it == xsData.begin()) return it->second;
            auto prev_it = std::prev(it);
            double x1 = prev_it->first, y1 = prev_it->second;
            double x2 = it->first, y2 = it->second;
            return y1 + (incidentEnergy - x1) * (y2 - y1) / (x2 - x1);
        };

        cap += getCachedXSValue(cache.captureData, incidentEnergy) * numberDensity;
        fis += getCachedXSValue(cache.fissionData, incidentEnergy) * numberDensity;
        scat += getCachedXSValue(cache.scatteringData, incidentEnergy) * numberDensity;
    }
}

std::vector<Material> readMaterialCards(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    std::vector<Material> materials;
    bool readingMaterialCards = false;

    while (getline(file, line)) {
        if (line.find("material card") != std::string::npos) {
            readingMaterialCards = true;
            continue;
        }
        if (line.find("end material card") != std::string::npos) break;

        if (readingMaterialCards) {
            std::istringstream iss(line);
            int num;
            double density;
            iss >> num >> density;
            Material material(num, density);

            int atomicNumber;
            double weightPercent;
            while (iss >> atomicNumber >> weightPercent) {
                material.addElement(atomicNumber, weightPercent);
            }
            material.calculateTotalMass();
            material.calculateNumberDensities();
            materials.push_back(material);
        }
    }
    return materials;
}
