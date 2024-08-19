#include <iostream>
#include <ctime>
#include "Surface.h"
#include "Material.h"
#include "CellCard.h"
#include "Calculation.h"
#include "Utils.h"

double r = 0.0;
double scat = 0.0, cap = 0.0, fis = 0.0, n2n = 0.0, scat_fraction, fis_fraction, n2n_fraction;
double tmp_erg, x, y, z, nps, mass_num;
int cycle = 0, skipped_cycle = 0, cycle_alarm = 0, area_out = 0, N_cnt = 0;
double fis_point[1000000][4] = {0, }, fis_point_tmp[1000000][4] = {0, },
        n2n_point[1000000][5] = {0, }, n2n_point_tmp[1000000][5] = {0, },
        fission_neutron_per_cycle[1000] = {0, }, n2n_neutron_per_cycle[1000] = {0, },
        keff[1000] = {0, }, k_sum = 0;

int main() {
    int fis_cnt = 0, fis_cnt_tmp = 0, cap_cnt = 0, cap_cnt_tmp = 0, n2n_cnt = 0, n2n_cnt_tmp = 0;
    clock_t start, finish;
    double duration;
    start = clock();
    std::string filePath = "D:/input/code_develop/mcsg_inputs/test0319.txt"; // 파일 경로 수정 필요
    auto surfaces = readSurfaceCards(filePath);
    auto materials = readMaterialCards(filePath);
    auto cellCards = readCellCards(filePath);
    Calculation calculation;
    readCalculationCards(filePath, calculation);
    nps = calculation.kcode[0], skipped_cycle = calculation.kcode[1], cycle = calculation.kcode[2], incident_erg = 1000000;
    x = 0, y = 0, z = 0;
    double cap = 0.0, fis = 0.0, scat = 0.0, n2n = 0.0;

    for (int i = 0; i < nps; i++) {
        tmp_erg = incident_erg;
        int cellNumber = findCellNumber(x, y, z, cellCards, surfaces);
        std::cout << cellNumber << std::endl;

        if (cellNumber != -1) {
            auto materialInfo = getMaterialInfo(cellNumber, cellCards, materials);
            if (!materialInfo.numberDensities.empty()) {
                materials[0].loadXSData("D:/input/XS/ENDF71/");
            }

            double cap = 0.0, fis = 0.0, scat = 0.0, n2n = 0.0; // 각 반복마다 초기화
            materials[0].calculateXsAndStore(tmp_erg, cap, fis, scat);
            std::cout << "Scattering XS: " << scat << std::endl;
            std::cout << "Capture XS: " << cap << std::endl;
            std::cout << "Fission XS: " << fis << std::endl;

            scat_fraction = scat / (scat + cap + fis + n2n);
            n2n_fraction = (scat + n2n) / (scat + cap + fis + n2n);
            fis_fraction = (scat + n2n + fis) / (scat + cap + fis + n2n);
            std::cout << scat_fraction << " " << n2n_fraction << " " << fis_fraction << std::endl;
        }
    }

    // XI_GEN 함수 호출
    XI_GEN();
    std::cout << "Generated XI: " << xi << std::endl;

    int cellNumber = findCellNumber(x, y, z, cellCards, surfaces);
    std::cout << cellNumber << std::endl;

    if (cellNumber != -1) {
        auto materialInfo = getMaterialInfo(cellNumber, cellCards, materials);
        if (!materialInfo.numberDensities.empty()) {
            materials[0].loadXSData("D:/input/XS/ENDF71/");

            double cap, fis, scat;
            materials[0].calculateXsAndStore(incident_erg, cap, fis, scat);
            std::cout << "Scattering XS: " << scat << std::endl;
            std::cout << "Capture XS: " << cap << std::endl;
            std::cout << "Fission XS: " << fis << std::endl;
        }
    }

    finish = clock(), duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout << "Computing Time : " << duration << "sec" << std::endl;
    return 0;
}
