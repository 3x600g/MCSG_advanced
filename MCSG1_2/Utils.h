#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include "Material.h"
#include "CellCard.h"



extern double xi, theta, phi, cos_theta, sin_theta, cos_phi, sin_phi, n_energy, incident_erg, FN_num;
extern const double pi;

struct MaterialInfo {
    double totalMass;
    std::vector<double> numberDensities;
};

double linearInterpolation(double x0, double y0, double x1, double y1, double x);
double readXsValue(const std::string& filename, double incidentEnergy);
void loadXSData(const std::string& base_path, const std::string& isotop, XSCache& cache);
double getCachedXSValue(const std::vector<std::pair<double, double>>& xsData, double incidentEnergy);
int findCellNumber(double x, double y, double z, const std::vector<CellCard>& cellCards, const std::vector<Surface>& surfaces);
MaterialInfo getMaterialInfo(int cellNumber, const std::vector<CellCard>& cellCards, const std::vector<Material>& materials);

void XI_GEN();
void ANGLE_GEN();
void FISSION_NEUTRON_NUM();
void N_ENERGY_GEN();
void NEUTRON_SPECTRA();

#endif // UTILS_H
