#ifndef SURFACE_H
#define SURFACE_H
#include <vector>
#include <string>
class Surface {
public:
    int s_number;
    std::string s_type;
    std::vector<double> s_parameters;

    Surface(int num, const std::string& typ, const std::vector<double>& params);
    bool isInside(double x, double y, double z) const;
};
std::vector<Surface> readSurfaceCards(const std::string& filePath);
#endif // SURFACE_H