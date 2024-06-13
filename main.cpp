#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <ctime>
#include <cmath>
#include <algorithm>

double xi, x, y, z, r = 0.0, FN_num, n_energy, scat = 0.0, cap = 0.0, fis = 0.0, n2n = 0.0, scat_fraction, fis_fraction, n2n_fraction;
double tmp_erg, incident_erg, theta, phi, cos_theta, sin_theta, cos_phi, sin_phi, nps, mass_num;
const double pi = 3.1415827;
int cycle = 0, skipped_cycle = 0, cycle_alarm = 0, area_out = 0, N_cnt = 0;
double fis_point[1000000][4] = {0, }, fis_point_tmp[1000000][4] = {0, }, n2n_point[1000000][5]={0,}, n2n_point_tmp[1000000][5]={0, },
        fission_neutron_per_cycle[1000] = {0, }, n2n_neutron_per_cycle[1000]={0, }, keff[1000]={0, }, k_sum = 0;
void XI_GEN(), ANGLE_GEN(), FISSION_NEUTRON_NUM(), N_ENERGY_GEN(), NEUTRON_SPECTRA();

class Surface {
public:
    int s_number;
    std::string s_type;
    std::vector<double> s_parameters; // 다양한 파라미터를 저장하기 위한 벡터

    // 생성자는 surface 번호, 타입(원, 원기둥, 직육면체..)과, 그리고 파라미터(반지름, 높이 ..)들을 받음.
    Surface(int num, const std::string& typ, const std::vector<double>& params)
            : s_number(num), s_type(typ), s_parameters(params) {}
    bool isInside(double x, double y, double z) const {
        if (s_type == "so") { // 구
            double radius = s_parameters[0];
            return x*x + y*y + z*z < radius*radius;
        } else if (s_type == "cy") { // 원기둥
            // 원기둥 중심 (cx, cy), 반지름 radius, 높이 h
            double cx = s_parameters[0], cy = s_parameters[1];
            double radius = s_parameters[6], h_vec = s_parameters[5];
            double dx = x - cx, dy = y - cy;
            return dx*dx + dy*dy < radius*radius && z >= 0 && z <= s_parameters[2]+h_vec; // 예시 조건, 구체적인 구현은 문제에 따라 달라질 수 있음
        } else if (s_type == "cu") { // 입방체
            // 입방체의 각 축별 최소값(minx, miny, minz)과 최대값(maxx, maxy, maxz)
            double minx = s_parameters[0], maxx = s_parameters[1];
            double miny = s_parameters[2], maxy = s_parameters[3];
            double minz = s_parameters[4], maxz = s_parameters[5];
            return x >= minx && x <= maxx && y >= miny && y <= maxy && z >= minz && z <= maxz;
        }
        // 다른 표면 유형에 대한 조건 추가 가능
        return false;
    }
};
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
        if (line.find("end surface card") != std::string::npos) break; // "end surface card" = 읽기를 중단.

        if (readingSurfaceCards) {
            std::istringstream iss(line);
            int num;
            std::string type;
            std::vector<double> params;
            double param;

            iss >> num >> type;
            while (iss >> param) {
                params.push_back(param); // 파라미터를 벡터에 추가.
            }
            surfaces.emplace_back(Surface(num, type, params));
        }
    }
    return surfaces;
}
class Material {
public:
    int m_number;
    double m_density;
    std::vector<std::pair<int, double>> composition; // 원소 번호와 weight percent
    double totalMass = 0.0;
    std::vector<double> numberDensities; // 원소별 number density

    Material(int num, double dens) : m_number(num), m_density(dens) {}

    void addElement(int atomicNumber, double weightPercent) {
        composition.push_back({atomicNumber, weightPercent});
    }
    void calculateTotalMass() {
        for (auto& element : composition) {
            int massNumber = element.first % 1000; // 원자 번호에서 질량 수를 추출
            double weight = element.second;
            totalMass += massNumber * (weight / 100.0);
        }
    }
    void calculateNumberDensities() {  // 위의 total mass를 통해 원소별 수밀도 계산
        const double NA = 6.02e23; // 아보가드로 수
        for (auto& element : composition) {
            int massNumber = element.first % 1000;
            double weightPercent = element.second;
            double nd = (weightPercent / 100.0) * m_density * NA / totalMass * 1e-24;
            numberDensities.push_back(nd);
        }
    }
    void calculateXsAndStore(double incidentEnergy);
};
std::vector<Material> readMaterialCards(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    std::vector<Material> materials;
    bool readingMaterialCards = false;

    while (getline(file, line)) {
        if (line.find("material card") != std::string::npos) {
            readingMaterialCards = true;
            // continue;
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

class CellCard {
public:
    int cellNumber;
    int c_materialNumber;
    std::vector<int> c_surfaces;
    bool isFinalBoundary = false;
    int finalBoundarySurface = -1;

    CellCard(int cellNum, int matNum) : cellNumber(cellNum), c_materialNumber(matNum) {}

    void addSurface(int surfaceNumber) {
        if (surfaceNumber == 0) {
            isFinalBoundary = true;
        } else if (isFinalBoundary) {
            finalBoundarySurface = surfaceNumber;
        } else {
            c_surfaces.push_back(surfaceNumber);
        }
    }

    bool isPointInCell(double x, double y, double z, const std::vector<Surface>& surfaces) const {
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

    bool isPointOutsideFinalBoundary(double x, double y, double z, const std::vector<Surface>& surfaces) const {
        if (finalBoundarySurface == -1) return false;

        const auto& surface = *std::find_if(surfaces.begin(), surfaces.end(),
                                            [this](const Surface& s) { return s.s_number == finalBoundarySurface; });

        return !surface.isInside(x, y, z);
    }
};
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
int detectCellForPoint(double x, double y, double z, const std::vector<CellCard>& cellCards, const std::vector<Surface>& surfaces) {
    for (const auto& cell : cellCards) {
        if (cell.isPointInCell(x, y, z, surfaces)) {
            return cell.cellNumber;
        }
    }

    for (const auto& cell : cellCards) {
        if (cell.isFinalBoundary && cell.isPointOutsideFinalBoundary(x, y, z, surfaces)) {
            return cell.cellNumber;
        }
    }

    return 0;
}

class Calculation {
public:
    std::vector<int> kcode; // kcode 정보를 저장할 벡터 (예: [3000, 60, 150])
    std::vector<float> ksrc; // ksrc 정보를 저장할 벡터 (예: [0, 0, 0])

    // kcode 정보를 파싱하여 저장하는 함수
    void parseKcode(const std::string& line) {
        std::istringstream iss(line);
        int value;
        while (iss >> value) {
            kcode.push_back(value);
        }
    }

    // ksrc 정보를 파싱하여 저장하는 함수
    void parseKsrc(const std::string& line) {
        std::istringstream iss(line);
        float value;
        while (iss >> value) {
            ksrc.push_back(value);
        }
    }
};
// 파일에서 "calculation card" 정보를 읽는 함수
void readCalculationCards(const std::string& filePath, Calculation& calculation) {
    std::ifstream file(filePath);
    std::string line;

    while (getline(file, line)) {
        if (line.find("kcode") != std::string::npos) {
            calculation.parseKcode(line.substr(line.find("kcode") + 6));
        } else if (line.find("ksource") != std::string::npos) {
            calculation.parseKsrc(line.substr(line.find("ksource") + 7));
        }
    }
}
void printMaterialInfo(int cellNumber, const std::vector<CellCard>& cellCards, const std::vector<Material>& materials) {
    // cellNumber에 해당하는 CellCard 찾기
    auto cellIt = std::find_if(cellCards.begin(), cellCards.end(),
                               [cellNumber](const CellCard& cell) { return cell.cellNumber == cellNumber; });

    if (cellIt == cellCards.end()) {
        std::cout << "Cell " << cellNumber << " not found." << std::endl;
        return;
    }

    // CellCard에서 물질 번호 찾기
    int materialNumber = cellIt->c_materialNumber;

    // 해당 물질 번호에 해당하는 Material 찾기
    auto materialIt = std::find_if(materials.begin(), materials.end(),
                                   [materialNumber](const Material& material) { return material.m_number == materialNumber; });

    if (materialIt == materials.end()) {
        std::cout << "Material for cell " << cellNumber << " not found." << std::endl;
        return;
    }

    // Material 정보 출력
    const Material& material = *materialIt;
    std::cout << "Material Number: " << material.m_number << std::endl;
    std::cout << "Density: " << material.m_density << " g/cm^3" << std::endl;
    std::cout << "Total Mass: " << material.totalMass << std::endl;
    std::cout << "Number Densities: ";
    for (double nd : material.numberDensities) {
        std::cout << nd << " ";
    }
    std::cout << std::endl;
}
// Material 정보를 가져와 저장하는 구조체 정의
struct MaterialInfo {
    double totalMass;
    std::vector<double> numberDensities;
};

// Material 정보를 가져오는 함수
MaterialInfo getMaterialInfo(int cellNumber, const std::vector<CellCard>& cellCards, const std::vector<Material>& materials) {
    MaterialInfo materialInfo = {0.0, {}};

    auto cellIt = std::find_if(cellCards.begin(), cellCards.end(),
                               [cellNumber](const CellCard& cell) { return cell.cellNumber == cellNumber; });
    if (cellIt != cellCards.end()) {
        int materialNumber = cellIt->c_materialNumber;
        std::cout << "Looking for Material Number: " << materialNumber << std::endl;

        auto materialIt = std::find_if(materials.begin(), materials.end(),
                                       [materialNumber](const Material& material) { return material.m_number == materialNumber; });
        if (materialIt != materials.end()) {
            std::cout << "Found Material. Number: " << materialIt->m_number << ", Total Mass: " << materialIt->totalMass << std::endl;
            materialInfo.totalMass = materialIt->totalMass;
            materialInfo.numberDensities = materialIt->numberDensities;
        } else {
            std::cout << "Material Number " << materialNumber << " not found in materials list." << std::endl;
        }
    } else {
        std::cout << "Cell Number " << cellNumber << " not found in cell cards." << std::endl;
    }

    return materialInfo;
}

double linearInterpolation(double x0, double y0, double x1, double y1, double x) {
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
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
void Material::calculateXsAndStore(double incidentEnergy) {
    std::vector<std::string> reactions = {"capture", "fission", "scattering"};
    std::string base_path = "D:/input/XS/ENDF71/"; // XS 파일이 있는 디렉토리 경로
    for (const auto& comp : composition) {
        std::string isotop = std::to_string(comp.first); // 동위원소 번호
        double weightPercent = comp.second; // 해당 동위원소의 weight percent
        double numberDensity = (weightPercent/100) * m_density * 6.02e+23 / totalMass * 1e-24; // number density 계산

        for (const auto& reaction : reactions) {
            std::string filename = base_path + isotop + "_" + reaction + ".txt";
            double xsValue = readXsValue(filename, incidentEnergy);

            // 계산하여 전역 변수에 저장
            if (reaction == "capture") {
                cap += xsValue * numberDensity;
            } else if (reaction == "fission") {
                fis += xsValue * numberDensity;
            } else if (reaction == "scattering") {
                scat += xsValue * numberDensity;
            }
        }
    }
}
void calculateAndStoreXS(std::vector<Material>& materials, double incidentEnergy) {
    // 결과 초기화
    cap = fis = scat = 0.0;

    for (auto& material : materials) {
        material.calculateXsAndStore(incidentEnergy);
    }
}

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
    nps = calculation.kcode[0], skipped_cycle = calculation.kcode[1], cycle = calculation.kcode[2], incident_erg=1.0;
    // x = calculation.ksrc[0], y = calculation.ksrc[1], z=calculation.ksrc[2];
    x = 8.742, y = 0, z = 0;
    int cellNum = detectCellForPoint(x, y, z, cellCards, surfaces);
    std::cout << cellNum << std::endl;


    finish = clock(), duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout << "Computing Time : " << duration << "sec" << std::endl;
    return 0;
}
void XI_GEN() {
    double vari;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(1, 9999);
    vari = dis(gen);
    xi = vari / 10000;
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
        XI_GEN(), N_ENERGY_GEN();
        p_n_energy = 0.4865 * sinh((sqrt(2 * n_energy)) * exp(-1 * n_energy));
        reaction1 = p_n_energy / p_max;
        if (reaction1 >= xi) {
            incident_erg = n_energy*1000000;
            break;
        }
    }
}