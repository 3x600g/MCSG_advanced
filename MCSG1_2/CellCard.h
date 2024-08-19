#ifndef CELLCARD_H
#define CELLCARD_H

#include <vector>
#include <string>
#include "Surface.h"

class CellCard {
public:
    int cellNumber;
    int c_materialNumber;
    std::vector<int> c_surfaces;
    bool isFinalBoundary = false;
    int finalBoundarySurface = -1;

    CellCard(int cellNum, int matNum);

    void addSurface(int surfaceNumber);
    bool isPointInCell(double x, double y, double z, const std::vector<Surface>& surfaces) const;
    bool isPointOutsideFinalBoundary(double x, double y, double z, const std::vector<Surface>& surfaces) const;
};

std::vector<CellCard> readCellCards(const std::string& filePath);

#endif // CELLCARD_H
