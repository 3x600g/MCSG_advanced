#include "Surface.h"
#include <fstream>
#include <sstream>

Surface::Surface(int num, const std::string& typ, const std::vector<double>& params)
        : s_number(num), s_type(typ), s_parameters(params) {}

bool Surface::isInside(double x, double y, double z) const {
    if (s_type == "so") {
        double radius = s_parameters[0];
        return x * x + y * y + z * z < radius * radius;
    } else if (s_type == "cy") {
        double cx = s_parameters[0], cy = s_parameters[1];
        double radius = s_parameters[6], h_vec = s_parameters[5];
        double dx = x - cx, dy = y - cy;
        return dx * dx + dy * dy < radius * radius && z >= 0 && z <= s_parameters[2] + h_vec;
    } else if (s_type == "cu") {
        double minx = s_parameters[0], maxx = s_parameters[1];
        double miny = s_parameters[2], maxy = s_parameters[3];
        double minz = s_parameters[4], maxz = s_parameters[5];
        return x >= minx && x <= maxx && y >= miny && y <= maxy && z >= minz && z <= maxz;
    }
    return false;
}

std::vector<Surface> readSurfaceCards(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    std::vector<Surface> surfaces;
    bool readingSurfaceCards = false;

    while (getline(file, line)) {
        if (line.find("surface card") != std::string::npos) {
            readingSurfaceCards = true;
            continue;
        }
        if (line.find("end surface card") != std::string::npos) break;

        if (readingSurfaceCards) {
            std::istringstream iss(line);
            int num;
            std::string type;
            std::vector<double> params;
            double param;

            iss >> num >> type;
            while (iss >> param) {
                params.push_back(param);
            }
            surfaces.emplace_back(Surface(num, type, params));
        }
    }
    return surfaces;
}
