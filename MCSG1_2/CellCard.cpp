#include "CellCard.h"
#include <fstream>
#include <sstream>
#include <algorithm>

CellCard::CellCard(int cellNum, int matNum) : cellNumber(cellNum), c_materialNumber(matNum) {}

void CellCard::addSurface(int surfaceNumber) {
    if (surfaceNumber == 0) {
        isFinalBoundary = true;
    } else if (isFinalBoundary) {
        finalBoundarySurface = surfaceNumber;
    } else {
        c_surfaces.push_back(surfaceNumber);
    }
}

bool CellCard::isPointInCell(double x, double y, double z, const std::vector<Surface>& surfaces) const {
    for (int surfaceNum : c_surfaces) {
        bool shouldBeInside = surfaceNum < 0;
        surfaceNum = std::abs(surfaceNum);

        const auto& surface = *std::find_if(surfaces.begin(), surfaces.end(),
                                            [surfaceNum](const Surface& s) { return s.s_number == surfaceNum; });
        if (surface.isInside(x, y, z) != shouldBeInside) {
            return false;
        }
    }
    return true;
}

bool CellCard::isPointOutsideFinalBoundary(double x, double y, double z, const std::vector<Surface>& surfaces) const {
    if (finalBoundarySurface == -1) return false;

    const auto& surface = *std::find_if(surfaces.begin(), surfaces.end(),
                                        [this](const Surface& s) { return s.s_number == finalBoundarySurface; });

    return !surface.isInside(x, y, z);
}

std::vector<CellCard> readCellCards(const std::string& filePath) {
    std::vector<CellCard> cellCards;
    std::ifstream file(filePath);
    std::string line;
    bool readingCellCards = false;
    while (getline(file, line)) {
        if (line.find("cell card") != std::string::npos) {
            readingCellCards = true;
            continue;
        }
        if (line.find("end cell card") != std::string::npos) break;

        if (readingCellCards) {
            std::istringstream iss(line);
            int cellNum, matNum, surfaceNum;
            iss >> cellNum >> matNum;
            CellCard cellCard(cellNum, matNum);

            while (iss >> surfaceNum) {
                cellCard.addSurface(surfaceNum);
            }
            cellCards.push_back(cellCard);
        }
    }
    return cellCards;
}
