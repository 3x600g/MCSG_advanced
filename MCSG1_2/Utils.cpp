#include "Utils.h"
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

double xi, theta, phi, cos_theta, sin_theta, cos_phi, sin_phi, n_energy, incident_erg, FN_num;
const double pi = 3.1415927;

double linearInterpolation(double x0, double y0, double x1, double y1, double x) {
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}
std::mt19937& get_rng() {
    static std::random_device rd; // 시드 생성기
    static std::mt19937 rng(rd()); // 난수 생성기
    return rng;
}
double generate_random() {
    static std::uniform_real_distribution<> dis(0.0, 1.0); // 0과 1 사이의 실수 난수 생성
    return dis(get_rng());
}

double readXsValue(const std::string& filename, double incidentEnergy) {
    std::ifstream file(filename);
    double energy, xs, prevEnergy, prevXs;
    bool first = true;

    if (file.is_open()) {
        while (file >> energy >> xs) {
            if (energy == incidentEnergy) {
                return xs; // 정확한 값이 파일에 있을 경우
            } else if (energy > incidentEnergy) {
                if (first) return xs; // 에너지가 범위를 벗어난 경우 첫 번째 XS 값 반환
                return linearInterpolation(prevEnergy, prevXs, energy, xs, incidentEnergy); // 선형 보간
            }
            prevEnergy = energy;
            prevXs = xs;
            first = false;
        }
    }
    return 0.0; // 파일이 없거나 에너지가 범위를 벗어난 경우
}

void loadXSData(const std::string& base_path, const std::string& isotop, XSCache& cache) {
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
}

double getCachedXSValue(const std::vector<std::pair<double, double>>& xsData, double incidentEnergy) {
    if (xsData.empty()) return 0.0;

    // 이진 탐색을 사용하여 incidentEnergy에 대한 값을 빠르게 찾습니다.
    auto it = std::lower_bound(xsData.begin(), xsData.end(), incidentEnergy,
                               [](const std::pair<double, double>& data, double energy) {
                                   return data.first < energy;
                               });

    // incidentEnergy가 범위를 벗어난 경우
    if (it == xsData.end()) {
        return xsData.back().second;
    }
    if (it == xsData.begin()) {
        return it->second;
    }

    // 선형 보간을 사용하여 정확한 값을 계산합니다.
    auto prev_it = std::prev(it);
    double x1 = prev_it->first;
    double y1 = prev_it->second;
    double x2 = it->first;
    double y2 = it->second;

    return y1 + (incidentEnergy - x1) * (y2 - y1) / (x2 - x1);
}



MaterialInfo getMaterialInfo(int cellNumber, const std::vector<CellCard>& cellCards, const std::vector<Material>& materials) {
    MaterialInfo materialInfo = {0.0, {}};

    auto cellIt = std::find_if(cellCards.begin(), cellCards.end(),
                               [cellNumber](const CellCard& cell) { return cell.cellNumber == cellNumber; });
    if (cellIt != cellCards.end()) {
        int materialNumber = cellIt->c_materialNumber;

        auto materialIt = std::find_if(materials.begin(), materials.end(),
                                       [materialNumber](const Material& material) { return material.m_number == materialNumber; });
        if (materialIt != materials.end()) {
            materialInfo.totalMass = materialIt->totalMass;
            materialInfo.numberDensities = materialIt->numberDensities;
        }
    }

    return materialInfo;
}

int findCellNumber(double x, double y, double z, const std::vector<CellCard>& cellCards, const std::vector<Surface>& surfaces) {
    for (const auto& cellCard : cellCards) {
        if (cellCard.isPointInCell(x, y, z, surfaces)) {
            return cellCard.cellNumber;
        }
    }

    for (const auto& cellCard : cellCards) {
        if (cellCard.isFinalBoundary && cellCard.isPointOutsideFinalBoundary(x, y, z, surfaces)) {
            return cellCard.cellNumber;
        }
    }
    return -1;
}

void XI_GEN() {
    xi = generate_random(); // 0과 1 사이의 난수 반환
}

void ANGLE_GEN() {
    XI_GEN();
    theta = acos(2 * xi - 1);
    theta = theta * 180 / pi;
    XI_GEN();
    phi = 2 * pi * xi;
    phi = phi * 180 / pi;

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);
}

void FISSION_NEUTRON_NUM() {
    XI_GEN();
    if (xi < 0.55) FN_num = 2;
    else FN_num = 3;
}

void N_ENERGY_GEN() {
    double vari1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(1, 1000000);
    vari1 = dis(gen);
    n_energy = vari1 / 100000;
}

void NEUTRON_SPECTRA() {
    double p_max = 0.35745, reaction1, p_n_energy;
    while (true) {
        XI_GEN();
        N_ENERGY_GEN();
        p_n_energy = 0.4865 * sinh((sqrt(2 * n_energy)) * exp(-1 * n_energy));
        reaction1 = p_n_energy / p_max;
        if (reaction1 >= xi) {
            incident_erg = n_energy * 1000000;
            break;
        }
    }
}
