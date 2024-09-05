#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<cmath>
using namespace std;

double U235cap[229302][2], U235scat[229302][2], U235fi[229302][2];
double U238cap[154426][2], U238scat[154426][2], U238fi[154426][2];
double O16cap[2985][2], O16scat[2985][2];
double H1cap[631][2], H1scat[631][2];
double incident_erg;
double cos_theta, sin_theta, cos_phi, sin_phi;
const double pi = 3.1415926535;
double theta, phi;
double xi;
double randx, randy, randz;
double r;
int simulation;
void U235_capture();
void U235_scattering();
void U235_fission();
void U238_capture();
void U238_scattering();
void U238_fission();
void O16_capture();
void O16_scattering();
void H1_capture();
void H1_scattering();
void xi_gen();
void rand_x();
void rand_y();
void rand_z();
void angle_gen();

int main() {
    U235_capture();
    U235_scattering();
    U238_capture();
    U238_scattering();
    O16_capture();
    O16_scattering();
    H1_capture();
    H1_scattering();
    U235_fission();
    U238_fission();

    double x = 0, y = 0, z = 0;
    double cnt = 0;
    int simulation;
    double tmp_erg;
    double reaction;
    double atom_num;

    cout << "Westing House 17x17 Assembly Neutron Shielding MonteCarlo Code" << endl;
    cout << endl;
    cout << "Heigt : 365.76cm, AssemblyPitch Length : 21.42cm " << endl;
    cout << endl;
    cout << "5wt% Enriched PWR Fuel" << endl;
    cout << endl;
    cout << "Developed by SG Cho"<< endl;
    cout << endl;
    cout << "Last Update : 23/04/04" << endl;
    cout << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << endl;
    cout << "Enter the Number of Simulations(nps) : ";
    cin >> simulation;
    cout << endl;
    cout << "Enter the Initial Incident Energy of Neutron(MeV) : ";
    cin >> tmp_erg;
    cout << endl;

    for (int i = 0; i < simulation; i++) {
        int idx1 = 0, idx2=0, idx3=0, idx4=0, idx5=0, idx6=0, idx7=0, idx8=0;
        double r = 0;
        double scat = 0, cap = 0;
        x = 0.0, y = 0.0, z = 0.0;
        /*rand_x();
        rand_y();
        rand_z();
        x = randx, y = randy, z = randz;
        */
        incident_erg = 1000000*tmp_erg;

        while (1) {
            atom_num = 0;
            if (incident_erg < 0.01) {
                break;
            }
            if (abs(x) >= 10.71 || abs(y) >= 10.71 || abs(z) >= 182.88) {
                cnt++;
                // cout << "----------------------------------- " << x << " " << y << " " << z;
                break;
            }

            int xc=0, yc=0;
            int xx = (int)(x*100);
            int yy = (int)(y*100);

            for (int i = -1071; i <= 1071 - 126; i += 126) {
                if ((i <= xx) && (xx <= i + 126)) {
                    xc = i + 63;
                }
                if ((i <= yy)&& (yy <= i + 126)) {
                    yc = i + 63;
                }
            }

            if (abs(x) <= 0.63) {
                if (abs(y) <= 0.63 || ( abs(abs(y) - 3.78) <= 0.63) || (abs(abs(y) - 7.56) <= 0.63)) {
                    double tmp1 = 0;
                    int idx1 = 0;
                    while (1) {
                        tmp1 = H1cap[idx1][0] - incident_erg;
                        if (tmp1 > 0) {
                            cap = 2 * H1cap[idx1][1];
                            break;
                        }
                        idx1++;
                    }
                    double tmp2 = 0;
                    int idx2 = 0;
                    while (1) {
                        tmp2 = O16cap[idx2][0] - incident_erg;
                        if (tmp2 > 0) {
                            cap = cap + O16cap[idx2][1];
                            break;
                        }
                        idx2++;
                    }
                    double tmp3 = 0;
                    int idx3 = 0;
                    while (1) {
                        tmp3 = H1scat[idx3][0] - incident_erg;
                        if (tmp3 > 0) {
                            scat = 2*H1scat[idx3][1];
                            break;
                        }
                        idx3++;
                    }
                    double tmp4 = 0;
                    int idx4 = 0;
                    while (1) {
                        tmp4 = O16scat[idx4][0] - incident_erg;
                        if (tmp4 > 0) {
                            scat = scat+ O16scat[idx4][1];
                            break;
                        }
                        idx4++;
                    }
                }
                atom_num = 18;
            }

            else if (abs(abs(x) - 3.78) <= 0.63) {
                if (abs(y) <= 0.63 || (abs(abs(y) - 3.78) <= 0.63) || (abs(abs(y) - 7.56) <= 0.63)) {
                    double tmp1 = 0;
                    int idx1 = 0;
                    while (1) {
                        tmp1 = H1cap[idx1][0] - incident_erg;
                        if (tmp1 > 0) {
                            cap = 2*H1cap[idx1][1];
                            break;
                        }
                        idx1++;
                    }
                    double tmp2 = 0;
                    int idx2 = 0;
                    while (1) {
                        tmp2 = O16cap[idx2][0] - incident_erg;
                        if (tmp2 > 0) {
                            cap = cap + O16cap[idx2][1];
                            break;
                        }
                        idx2++;
                    }
                    double tmp3 = 0;
                    int idx3 = 0;
                    while (1) {
                        tmp3 = H1scat[idx3][0] - incident_erg;
                        if (tmp3 > 0) {
                            scat = 2*H1scat[idx3][1];
                            break;
                        }
                        idx3++;
                    }
                    double tmp4 = 0;
                    int idx4 = 0;
                    while (1) {
                        tmp4 = O16scat[idx4][0] - incident_erg;
                        if (tmp4 > 0) {
                            scat = scat + O16scat[idx4][1];
                            break;
                        }
                        idx4++;
                    }
                }
                atom_num = 18;
            }

            else if (abs(abs(x) - 7.56) <= 0.63) {
                if (abs(y) <= 0.63 || (abs(abs(y) - 3.78) <= 0.63)) {
                    double tmp1 = 0;
                    int idx1 = 0;
                    while (1) {
                        tmp1 = H1cap[idx1][0] - incident_erg;
                        if (tmp1 > 0) {
                            cap = 2 * H1cap[idx1][1];
                            break;
                        }
                        idx1++;
                    }
                    double tmp2 = 0;
                    int idx2 = 0;
                    while (1) {
                        tmp2 = O16cap[idx2][0] - incident_erg;
                        if (tmp2 > 0) {
                            cap = cap + O16cap[idx2][1];
                            break;
                        }
                        idx2++;
                    }
                    double tmp3 = 0;
                    int idx3 = 0;
                    while (1) {
                        tmp3 = H1scat[idx3][0] - incident_erg;
                        if (tmp3 > 0) {
                            scat = 2*H1scat[idx3][1];
                            break;
                        }
                        idx3++;
                    }
                    double tmp4 = 0;
                    int idx4 = 0;
                    while (1) {
                        tmp4 = O16scat[idx4][0] - incident_erg;
                        if (tmp4 > 0) {
                            scat = scat + O16scat[idx4][1];
                            break;
                        }
                        idx4++;
                    }
                }
                atom_num = 18;
            }

            else if ((abs(abs(x) - 6.3) <= 0.63) && (abs(abs(y) - 6.3) <= 0.63)) {
                double tmp1 = 0;
                int idx1 = 0;
                while (1) {
                    tmp1 = H1cap[idx1][0] - incident_erg;
                    if (tmp1 > 0) {
                        cap = 2 * H1cap[idx1][1];
                        break;
                    }
                    idx1++;
                }
                double tmp2 = 0;
                int idx2 = 0;
                while (1) {
                    tmp2 = O16cap[idx2][0] - incident_erg;
                    if (tmp2 > 0) {
                        cap = cap + O16cap[idx2][1];
                        break;
                    }
                    idx2++;
                }
                double tmp3 = 0;
                int idx3 = 0;
                while (1) {
                    tmp3 = H1scat[idx3][0] - incident_erg;
                    if (tmp3 > 0) {
                        scat = 2*H1scat[idx3][1];
                        break;
                    }
                    idx3++;
                }
                double tmp4 = 0;
                int idx4 = 0;
                while (1) {
                    tmp4 = O16scat[idx4][0] - incident_erg;
                    if (tmp4 > 0) {
                        scat += O16scat[idx4][1];
                        break;
                    }
                    idx4++;
                }
                atom_num = 18;
            }

            else if (sqrt((xx - xc) * (xx - xc) + (yy - yc) * (yy - yc)) > 41) {
                double tmp1 = 0;
                int idx1 = 0;
                while (1) {
                    tmp1 = H1cap[idx1][0] - incident_erg;
                    if (tmp1 > 0) {
                        cap = 2 * H1cap[idx1][1];
                        break;
                    }
                    idx1++;
                }
                double tmp2 = 0;
                int idx2 = 0;
                while (1) {
                    tmp2 = O16cap[idx2][0] - incident_erg;
                    if (tmp2 > 0) {
                        cap += O16cap[idx2][1];
                        break;
                    }
                    idx2++;
                }
                double tmp3 = 0;
                int idx3 = 0;
                while (1) {
                    tmp3 = H1scat[idx3][0] - incident_erg;
                    if (tmp3 > 0) {
                        scat = 2*H1scat[idx3][1];
                        break;
                    }
                    idx3++;
                }
                double tmp4 = 0;
                int idx4 = 0;
                while (1) {
                    tmp4 = O16scat[idx4][0] - incident_erg;
                    if (tmp4 > 0) {
                        scat += O16scat[idx4][1];
                        break;
                    }
                    idx4++;
                }
                atom_num = 18;
            }


            else {
                double tmp1 = 0;
                int idx1 = 0;
                while (1) {
                    tmp1 = O16cap[idx1][0] - incident_erg;
                    if (tmp1 > 0) {
                        cap = 2 * O16cap[idx1][1];
                        break;
                    }
                    idx1++;
                }
                double tmp2 = 0;
                int idx2 = 0;
                while (1) {
                    tmp2 = U235cap[idx2][0] - incident_erg;
                    if (tmp2 > 0) {
                        cap = cap+ U235cap[idx2][1];
                        break;
                    }
                    idx2++;
                }
                double tmp3 = 0;
                int idx3 = 0;
                while (1) {
                    tmp3 = U238cap[idx3][0] - incident_erg;
                    if (tmp3 > 0) {
                        cap =cap+ U238cap[idx3][1];
                        break;
                    }
                    idx3++;
                }
                double tmp4 = 0;
                int idx4 = 0;
                while (1) {
                    tmp4 = O16scat[idx4][0] - incident_erg;
                    if (tmp4 > 0) {
                        scat = 2 * O16scat[idx4][1];
                        break;
                    }
                    idx4++;
                }
                double tmp5 = 0;
                int idx5 = 0;
                while (1) {
                    tmp5 = U235scat[idx5][0] - incident_erg;
                    if (tmp5 > 0) {
                        scat = scat+U235scat[idx5][1];
                        break;
                    }
                    idx5++;
                }
                double tmp6 = 0;
                int idx6 = 0;
                while (1) {
                    tmp6 = U238scat[idx6][0] - incident_erg;
                    if (tmp6 > 0) {
                        scat = scat+U238scat[idx6][1];
                        break;
                    }
                    idx6++;
                }
                double tmp7 = 0;
                int idx7 = 0;
                while (1) {
                    tmp7 = U235fi[idx7][0] - incident_erg;
                    if (tmp7 > 0) {
                        cap = cap + U235fi[idx7][1];
                        break;
                    }
                    idx7++;
                }
                double tmp8 = 0;
                int idx8 = 0;
                while (1) {
                    tmp8 = U238fi[idx8][0] - incident_erg;
                    if (tmp8 > 0) {
                        cap = cap + U238fi[idx8][1];
                        break;
                    }
                    idx8++;
                }
                atom_num = 269.85;
            }
            reaction = scat / (scat + cap);
            xi_gen();
            if (reaction >= xi) {
                xi_gen();
                angle_gen();
                r = -1 * log(xi) / (scat + cap);
                x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                incident_erg = incident_erg * pow((cos_theta + sqrt(pow(atom_num, 2) - pow(sin_theta, 2))), 2) / pow(atom_num + 1, 2);
            }
            else {
                break;
            }
        }
    }
    cout << "-------------------------result-----------------------------" << endl;
    cout << " " << endl;
    cout << "Unshielded Neutron Ratio : " << cnt / simulation * 100 << "%" << endl;
    cout << " " << endl;
    cout << simulation << "Unshielded Number of Neutron  : " << cnt << endl;
    cout << " " << endl;
    return 0;
}

void U235_capture() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/U235_capture.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> U235cap[idx][0];
        fin >> U235cap[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void U235_scattering() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/U235_elastic.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> U235scat[idx][0];
        fin >> U235scat[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void U235_fission() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/U235_fission.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> U235fi[idx][0];
        fin >> U235fi[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void U238_capture() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/U238_capture.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> U238cap[idx][0];
        fin >> U238cap[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void U238_scattering() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/U238_elastic.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> U238scat[idx][0];
        fin >> U238scat[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void U238_fission() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/U238_fission.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> U238fi[idx][0];
        fin >> U238fi[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void O16_capture() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/O16_capture.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> O16cap[idx][0];
        fin >> O16cap[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void O16_scattering() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/O16_elastic.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> O16scat[idx][0];
        fin >> O16scat[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void H1_capture() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/H1_capture.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> H1cap[idx][0];
        fin >> H1cap[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void H1_scattering() {
    int idx = 0;
    ifstream fin("D:/input/XS/WH1717_Shielding/H1_elastic.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> H1scat[idx][0];
        fin >> H1scat[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void xi_gen() {
    double vari;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(1, 9999);
    vari = dis(gen);
    xi = vari / 10000;
}

void rand_x() {
    double vari1;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(-1071, 1071);
    vari1 = dis(gen);
    randx = vari1 / 100;
}
void rand_y() {
    double vari2;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(-1071, 1071);
    vari2 = dis(gen);
    randy = vari2 / 100;
}
void rand_z() {
    double vari3;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(-18288, 18288);
    vari3 = dis(gen);
    randz = vari3 / 100;
}

void angle_gen() {
    xi_gen();
    theta = acos(2 * xi - 1);
    theta = theta * 180 / pi;
    xi_gen();
    phi = 2 * pi * xi;
    phi = phi * 180 / pi;

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);
}