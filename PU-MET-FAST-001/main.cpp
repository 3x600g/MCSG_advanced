#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<ctime>
#include<string>
using namespace std;
std::mt19937& get_rng() {
    static std::random_device rd; // 시드 생성기
    static std::mt19937 rng(rd()); // 난수 생성기
    return rng;
}
double generate_random() {
    static std::uniform_real_distribution<> dis(0.0, 1.0); // 0과 1 사이의 실수 난수 생성
    return dis(get_rng());
}
double Pu239cap[71716][2], Pu239scat[71716][2], Pu239fi[71716][2], Pu239n2n[124][2];
double Pu240cap[48276][2], Pu240scat[48276][2], Pu240fi[48276][2], Pu240n2n[94][2];
double xi, x, y, z, r = 0.0, FN_num, N_cnt = 0, n_energy, scat, cap, fis, scat_fraction, fis_fraction, n2n_fraction, n2n;
double tmp_erg, incident_erg, theta, phi, cos_theta, sin_theta, cos_phi, sin_phi, nps;
double pi = 3.1415827, atom_num = 239.045, boundary = 40.768225;
double fis_point[1000000][4] = {0, }, fis_point_tmp[1000000][4] = {0, }, n2n_point[1000000][5]={0,}, n2n_point_tmp[1000000][5]={0, },
        fission_neutron_per_cycle[1000] = {0, }, n2n_neutron_per_cycle[1000]={0, }, keff[1000]={0, }, k_sum = 0;
int cycle, skipped, cycle_alarm = 0, area_out = 0, collision = 0, std_fis = 0, std_nps = 0;
string file_path("D:/input/XS/Jezebel_fast/");
// string file_path("/Users/geun/CLionProjects/GODIVA_231031/");
void XI_GEN(), ANGLE_GEN(), FISSION_NEUTRON_NUM();
double N_ENERGY_GEN(), P_N_ENERGY();
vector<double> FISSION_ENERGY(int num_simulations);
void Pu239_capture(), Pu240_capture(), Pu239_scattering(), Pu240_scattering(), Pu239_fission(), Pu240_fission(), XS_READER(), Pu239_n2n(), Pu240_n2n();
int main() {
    clock_t start, finish;
    double duration;
    start = clock();
    int fis_cnt = 0, fis_cnt_tmp = 0, cap_cnt = 0, cap_cnt_tmp = 0, n2n_cnt = 0, n2n_cnt_tmp = 0, total_fission = 0;
    Pu239_capture(), Pu240_capture(), Pu239_scattering(), Pu240_scattering(), Pu239_fission(), Pu240_fission(), Pu239_n2n(), Pu240_n2n();

    cout << "Enter the number of Cycles : ";
    cin >> cycle;
    cout << "Enter the number of Skipped Cycle : ";
    cin >> skipped;
    cout << "Enter the Number of Neutrons per Generation : ";
    cin >> nps;
    cout << "Enter the Initial Incident Energy (MeV) : ";
    cin >> tmp_erg;
    cout << "Importing Pu-239 ENDF7.1 Cross-section ..." << endl;
    cout << "Importing Pu-240 ENDF7.1 Cross-section ..." << endl;

    total_fission = int(nps)*cycle;
    vector<double> fission_energy = FISSION_ENERGY(total_fission);
    int fission_energy_count = 0;

    for(int i = 0; i < nps; i++){
        incident_erg = tmp_erg * 1000000;
        XS_READER();
        scat_fraction = scat / (scat + cap + fis + n2n), n2n_fraction = (scat + n2n)/(scat + cap + fis + n2n), fis_fraction = (scat + n2n + fis) / (scat + cap + fis + n2n);
        XI_GEN();
        x = 0.0, y = 0.0, z = 0.0;
        if(xi <= scat_fraction){ // scattering...
            ANGLE_GEN();
            XI_GEN();
            r = -1 * log(xi) / (scat + cap + fis + n2n);
            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
            if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                area_out++;
                continue;
            }
            incident_erg *= atom_num / (atom_num + 1);
            collision++;
            while(true){
                XS_READER();
                scat_fraction = scat / (scat + cap + fis + n2n), n2n_fraction = (scat + n2n)/(scat + cap + fis + n2n), fis_fraction = (scat + fis + n2n) / (scat + cap + fis + n2n);
                XI_GEN();
                if(xi >= fis_fraction) { //capture
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (scat + cap + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                        cap_cnt++, cap_cnt_tmp++;
                        collision++;
                    }
                    else{
                        area_out++;
                    }
                    break;
                }
                else if(xi >= n2n_fraction){ // fission
                    collision++;
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (cap + scat + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                        FISSION_NEUTRON_NUM();
                        // cout << N_cnt << endl;
                        N_cnt += FN_num;
                        fis_cnt++, fis_cnt_tmp++;
                        fis_point[fis_cnt - 1][0] = fis_cnt_tmp;
                        fis_point[fis_cnt - 1][1] = x;
                        fis_point[fis_cnt - 1][2] = y;
                        fis_point[fis_cnt - 1][3] = z;
                    }
                    else{
                        area_out++;
                    }
                    break;
                }
                else if(xi >= scat_fraction){ // n2n reaction
                    collision++;
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (cap + scat + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                        n2n_cnt++, n2n_cnt_tmp++;
                        n2n_point[n2n_cnt - 1][0] = n2n_cnt_tmp;
                        n2n_point[n2n_cnt - 1][1] = x;
                        n2n_point[n2n_cnt - 1][2] = y;
                        n2n_point[n2n_cnt - 1][3] = z;
                        n2n_point[n2n_cnt - 1][4] = incident_erg;
                    }
                    else{
                        area_out++;
                    }
                    break;
                }
                else{ // scattering ... while other reaction occur
                    collision++;
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (scat + cap + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                        incident_erg *= atom_num / (atom_num + 1);
                    }
                    else {
                        area_out++;
                        break;
                    }
                }
            }
        }
        else if(xi <= n2n_fraction){ // n2n reaction
            ANGLE_GEN();
            XI_GEN();
            r = -1 * log(xi) / (cap + scat + fis + n2n);
            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
            if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                area_out++;
                continue;
            }
            n2n_cnt++, n2n_cnt_tmp++;
            n2n_point[n2n_cnt - 1][0] = n2n_cnt_tmp;
            n2n_point[n2n_cnt - 1][1] = x;
            n2n_point[n2n_cnt - 1][2] = y;
            n2n_point[n2n_cnt - 1][3] = z;
            n2n_point[n2n_cnt - 1][4] = incident_erg;
        }
        else if(xi <= fis_fraction){ // fission
            ANGLE_GEN();
            XI_GEN();
            r = -1 * log(xi) / (cap + scat + fis + n2n);
            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
            if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                area_out++;
                continue;
            }
            FISSION_NEUTRON_NUM();
            N_cnt += FN_num;
            fis_cnt++, fis_cnt_tmp++, std_fis++;
            fis_point[fis_cnt - 1][0] = fis_cnt_tmp;
            fis_point[fis_cnt - 1][1] = x;
            fis_point[fis_cnt - 1][2] = y;
            fis_point[fis_cnt - 1][3] = z;
        }
        else{  // capture
            ANGLE_GEN();
            XI_GEN();
            r = -1 * log(xi) / (cap + scat +fis + n2n);
            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
            if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                area_out++;
                continue;
            }
            cap_cnt++, cap_cnt_tmp++;
        }
    }

    fission_neutron_per_cycle[cycle_alarm] = fis_cnt_tmp;
    n2n_neutron_per_cycle[cycle_alarm] = n2n_cnt_tmp;
    fis_cnt_tmp = 0, fis_cnt = 0, cycle_alarm = 1;
    n2n_cnt_tmp = 0, n2n_cnt = 0;
    for (int i = 0; i < 1000000; i++) {
        for (int j = 0; j < 4; j++) {
            fis_point_tmp[i][j] = fis_point[i][j];
            n2n_point_tmp[i][j] = n2n_point[i][j];
        }
    }
    for (int i = 0; i < 1000000; i++) {
        for (int j = 0; j < 4; j++) {
            fis_point[i][j] = 0;
            n2n_point[i][j] = 0;
        }
    }
    cout << "fission count : " <<fission_neutron_per_cycle[0] << " capture count : " << cap_cnt << "  godiva area out : " << area_out << "      n2n : " << n2n_neutron_per_cycle[0] << endl;
    cout << fixed;
    cout.precision(5);
    cout << "1 cycle k-effective value : " << N_cnt/nps << endl;
    // cout << n2n_neutron_per_cycle[0] << endl;
    keff[0] = N_cnt/nps;

    cap_cnt = 0, area_out = 0;
    N_cnt =0;

    for(int i = 1; i < cycle; i++){
        int fis_idx = 0, n2n_idx = 0;
        for(int j = 0; j < nps; j++){
            incident_erg = fission_energy[fission_energy_count] * 1000000;
            fission_energy_count++;
            XS_READER();
            scat_fraction = scat / (scat + cap + fis + n2n), n2n_fraction = (scat + n2n)/(scat + cap + fis + n2n), fis_fraction = (scat + n2n + fis) / (scat + cap + fis + n2n);
            XI_GEN();
            x = fis_point_tmp[fis_idx][1], y = fis_point_tmp[fis_idx][2], z = fis_point_tmp[fis_idx][3];
            fis_idx++;
            if(fis_idx>=fission_neutron_per_cycle[i-1]){
                fis_idx = 0;
            }
            if(xi <= scat_fraction){ // scattering ...
                ANGLE_GEN();
                XI_GEN();
                r = -1 * log(xi) / (scat + cap + fis + n2n);
                x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                    area_out++;
                    continue;
                }
                incident_erg *= atom_num / (atom_num + 1);
                collision++;
                while(true){
                    XS_READER();
                    scat_fraction = scat / (scat + cap + fis + n2n), n2n_fraction = (scat + n2n)/(scat + cap + fis + n2n), fis_fraction = (scat + n2n + fis) / (scat + cap + fis + n2n);
                    XI_GEN();
                    if(xi >= fis_fraction) { //capture
                        ANGLE_GEN();
                        XI_GEN();
                        r = -1 * log(xi) / (scat + cap + fis + n2n);
                        x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                        if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                            cap_cnt++, cap_cnt_tmp++;
                            collision++;
                        }
                        else{
                            area_out++;
                        }
                        break;
                    }
                    else if(xi >= n2n_fraction){ // fission
                        collision++;
                        ANGLE_GEN();
                        XI_GEN();
                        r = -1 * log(xi) / (cap + scat + fis + n2n);
                        x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                        if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                            FISSION_NEUTRON_NUM();
                            N_cnt += FN_num;
                            fis_cnt++, fis_cnt_tmp++, std_fis++;
                            fis_point[fis_cnt - 1][0] = fis_cnt_tmp;
                            fis_point[fis_cnt - 1][1] = x;
                            fis_point[fis_cnt - 1][2] = y;
                            fis_point[fis_cnt - 1][3] = z;
                        }
                        else{
                            area_out++;
                        }
                        break;
                    }
                    else if(xi >= scat_fraction){  // n2n reaction
                        collision++;
                        ANGLE_GEN();
                        XI_GEN();
                        r = -1 * log(xi) / (cap + scat + fis + n2n);
                        x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                        if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                            n2n_cnt++, n2n_cnt_tmp++;
                            n2n_point[n2n_cnt - 1][0] = n2n_cnt_tmp;
                            n2n_point[n2n_cnt - 1][1] = x;
                            n2n_point[n2n_cnt - 1][2] = y;
                            n2n_point[n2n_cnt - 1][3] = z;
                            n2n_point[n2n_cnt - 1][4] = incident_erg;
                        }
                        else{
                            area_out++;
                        }
                        break;
                    }
                    else{  // scattering
                        collision++;
                        ANGLE_GEN();
                        XI_GEN();
                        r = -1 * log(xi) / (scat + cap + fis + n2n);
                        x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                        if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                            incident_erg *= atom_num / (atom_num + 1);
                        }
                        else {
                            area_out++;
                            break;
                        }
                    }
                }
            }
            else if(xi <= n2n_fraction){ // n2n reaction
                ANGLE_GEN();
                XI_GEN();
                r = -1 * log(xi) / (cap + scat + fis + n2n);
                x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                    area_out++;
                    continue;
                }
                n2n_cnt++, n2n_cnt_tmp++;
                n2n_point[n2n_cnt - 1][0] = n2n_cnt_tmp;
                n2n_point[n2n_cnt - 1][1] = x;
                n2n_point[n2n_cnt - 1][2] = y;
                n2n_point[n2n_cnt - 1][3] = z;
                n2n_point[n2n_cnt - 1][4] = incident_erg;
            }
            else if(xi <= fis_fraction){ // fission reaction
                ANGLE_GEN();
                XI_GEN();
                r = -1 * log(xi) / (cap + scat + fis + n2n);
                x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                    area_out++;
                    continue;
                }
                FISSION_NEUTRON_NUM();
                N_cnt += FN_num;
                fis_cnt++, fis_cnt_tmp++, std_fis++;
                fis_point[fis_cnt - 1][0] = fis_cnt_tmp;
                fis_point[fis_cnt - 1][1] = x;
                fis_point[fis_cnt - 1][2] = y;
                fis_point[fis_cnt - 1][3] = z;
            }
            else{  // capture
                ANGLE_GEN();
                XI_GEN();
                r = -1 * log(xi) / (cap + scat +fis + n2n);
                x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                    area_out++;
                    continue;
                }
                cap_cnt++, cap_cnt_tmp++;
            }
        }
        for(int j = 0; j<n2n_neutron_per_cycle[i-1]; j++){  // for n2n reaction physics
            for(int k = 0; k < 2; k++){
                incident_erg = (n2n_point[n2n_idx][4]) * atom_num / (atom_num + 1);
                XS_READER();
                scat_fraction = scat / (scat + cap + fis + n2n), n2n_fraction = (scat + n2n)/(scat + cap + fis + n2n), fis_fraction = (scat + n2n + fis) / (scat + cap + fis + n2n);
                XI_GEN();
                x = n2n_point_tmp[n2n_idx][1], y = n2n_point_tmp[n2n_idx][2], z = n2n_point_tmp[n2n_idx][3];
                n2n_idx++;
                if(xi <= scat_fraction){ // scattering ...
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (scat + cap + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                        area_out++;
                        continue;
                    }
                    incident_erg *= atom_num / (atom_num + 1);
                    collision++;
                    while(true){
                        XS_READER();
                        scat_fraction = scat / (scat + cap + fis + n2n), n2n_fraction = (scat + n2n)/(scat + cap + fis + n2n), fis_fraction = (scat + n2n + fis) / (scat + cap + fis + n2n);
                        XI_GEN();
                        if(xi >= fis_fraction) { //capture
                            ANGLE_GEN();
                            XI_GEN();
                            r = -1 * log(xi) / (scat + cap + fis + n2n);
                            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                            if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                                cap_cnt++, cap_cnt_tmp++;
                                collision++;
                            }
                            else{
                                area_out++;
                            }
                            break;
                        }
                        else if(xi >= n2n_fraction){ // fission
                            collision++;
                            ANGLE_GEN();
                            XI_GEN();
                            r = -1 * log(xi) / (cap + scat + fis + n2n);
                            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                            if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                                FISSION_NEUTRON_NUM();
                                N_cnt += FN_num;
                                fis_cnt++, fis_cnt_tmp++, std_fis++;
                                fis_point[fis_cnt - 1][0] = fis_cnt_tmp;
                                fis_point[fis_cnt - 1][1] = x;
                                fis_point[fis_cnt - 1][2] = y;
                                fis_point[fis_cnt - 1][3] = z;
                            }
                            else{
                                area_out++;
                            }
                            break;
                        }
                        else if(xi >= scat_fraction){  // n2n reaction
                            collision++;
                            ANGLE_GEN();
                            XI_GEN();
                            r = -1 * log(xi) / (cap + scat + fis + n2n);
                            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                            if(pow(x,2)+pow(y,2)+pow(z,2) < boundary) {
                                n2n_cnt++, n2n_cnt_tmp++;
                                n2n_point[n2n_cnt - 1][0] = n2n_cnt_tmp;
                                n2n_point[n2n_cnt - 1][1] = x;
                                n2n_point[n2n_cnt - 1][2] = y;
                                n2n_point[n2n_cnt - 1][3] = z;
                                n2n_point[n2n_cnt - 1][4] = incident_erg;
                            }
                            else{
                                area_out++;
                            }
                            break;
                        }
                        else{  // scattering
                            ANGLE_GEN();
                            XI_GEN();
                            r = -1 * log(xi) / (cap + scat +fis + n2n);
                            x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                            if(pow(x,2)+pow(y,2)+pow(z,2) <= boundary){
                                incident_erg *= atom_num / (atom_num + 1);
                            }
                            else{
                                area_out++;
                                continue;
                            }
                        }
                    }
                }
                else if(xi <= n2n_fraction){ // n2n reaction
                    collision++;
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (cap + scat + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) > boundary){
                        area_out++;
                        continue;
                    }
                    n2n_cnt++, n2n_cnt_tmp++;
                    n2n_point[n2n_cnt - 1][0] = n2n_cnt_tmp;
                    n2n_point[n2n_cnt - 1][1] = x;
                    n2n_point[n2n_cnt - 1][2] = y;
                    n2n_point[n2n_cnt - 1][3] = z;
                    n2n_point[n2n_cnt - 1][4] = incident_erg;
                }
                else if(xi <= fis_fraction){ // fission reaction
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (cap + scat + fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) <= boundary){
                        collision++;
                        FISSION_NEUTRON_NUM();
                        N_cnt += FN_num;
                        fis_cnt++, fis_cnt_tmp++, std_fis++;
                        fis_point[fis_cnt - 1][0] = fis_cnt_tmp;
                        fis_point[fis_cnt - 1][1] = x;
                        fis_point[fis_cnt - 1][2] = y;
                        fis_point[fis_cnt - 1][3] = z;
                    }
                    else{
                        area_out++;
                    }
                }
                else{ // capture
                    ANGLE_GEN();
                    XI_GEN();
                    r = -1 * log(xi) / (cap + scat +fis + n2n);
                    x += r * sin_theta * cos_phi, y += r * sin_theta * sin_phi, z += r * cos_theta;
                    if(pow(x,2)+pow(y,2)+pow(z,2) <= boundary){
                        cap_cnt++, cap_cnt_tmp++;
                    }
                    else{
                        area_out++;
                    }
                }
            }
        }

        fission_neutron_per_cycle[i] = fis_cnt_tmp;
        n2n_neutron_per_cycle[i] = n2n_cnt_tmp;

        // cout << "fission count : " << fission_neutron_per_cycle[i] << " capture count : " << cap_cnt << "  godiva area out : " << area_out << "    n2n : "<< n2n_neutron_per_cycle[i] << endl;
        cout << fixed;
        cout.precision(5);
        cout << i+1 <<" cycle k-effective value : " << N_cnt/nps << endl;
        // cout << i+1 <<" cycle k-effective value : " << fission_neutron_per_cycle[i]*2.5 / nps << endl;
        keff[i] = N_cnt/nps;
        N_cnt = 0;
        cap_cnt = 0, area_out = 0;
        fis_cnt = 0, fis_cnt_tmp = 0, collision = 0;
        n2n_cnt = 0, n2n_cnt_tmp = 0;
        cycle_alarm++;
        // cout << "===================================================================" << endl;

        for (int l = 0; l < 1000000; l++) {
            for (int m = 0; m < 4; m++) {
                fis_point_tmp[l][m] = fis_point[l][m];
                n2n_point_tmp[l][m] = n2n_point[l][m];
            }
        }
        for (int l = 0; l < 1000000; l++) {
            for (int m = 0; m < 4; m++) {
                fis_point[l][m] = 0;
                n2n_point[l][m] = 0;
            }// 다음 cycle fission 좌표 저장 위해 bank 초기화
        }
    }
    for(int i=skipped; i<cycle; i++){
        k_sum += keff[i];
    }
    cout << "average k-effective value : " << k_sum/(cycle-skipped) << endl;
    finish = clock(), duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << "Computing Time : " << duration << "sec" << endl;
    return 0;
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
    if (xi <= 0.01) FN_num = 2;
    else FN_num = 3;
    // cout << FN_num << endl;

    // FN_num = 3;
}
double N_ENERGY_GEN() {
    return 10 * generate_random();
}
double P_N_ENERGY(double n_energy) {
    return 0.355 * sinh(sqrt(2.842 * n_energy)) * exp(-1 * n_energy / 0.966);
}
vector<double> FISSION_ENERGY(int num_simulations) {
    vector<double> n_energy_values; // n_energy 값을 저장할 벡터
    const double p_max = 0.291; // 최대 확률 밀도 값

    for(int i = 0; i < num_simulations; ++i) {
        double n_energy, p_n_energy_val;

        while (true) {
            n_energy = N_ENERGY_GEN();
            p_n_energy_val = P_N_ENERGY(n_energy);
            XI_GEN();

            if (xi * p_max <= p_n_energy_val) {
                n_energy_values.push_back(n_energy); // 조건을 만족하면 벡터에 저장
                break; // 다음 시뮬레이션으로 넘어감
            }
        }
    }

    return n_energy_values;
}
void Pu239_capture() {
    int idx = 0;
    ifstream fin(file_path+"Pu239_capture.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (1) {
        fin >> Pu239cap[idx][0];
        fin >> Pu239cap[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu239_scattering() {
    int idx = 0;
    ifstream fin(file_path+"Pu239_elastic.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu239scat[idx][0];
        fin >> Pu239scat[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu239_fission() {
    int idx = 0;
    ifstream fin(file_path+"Pu239_fission.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu239fi[idx][0];
        fin >> Pu239fi[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu240_capture() {
    int idx = 0;
    ifstream fin(file_path+"Pu240_capture.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu240cap[idx][0];
        fin >> Pu240cap[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu240_scattering() {
    int idx = 0;
    ifstream fin(file_path+"Pu240_elastic.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu240scat[idx][0];
        fin >> Pu240scat[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu240_fission() {
    int idx = 0;
    ifstream fin(file_path+"Pu240_fission.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu240fi[idx][0];
        fin >> Pu240fi[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu239_n2n(){
    int idx = 0;
    ifstream fin(file_path+"Pu239_n2n.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu239n2n[idx][0];
        fin >> Pu239n2n[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}
void Pu240_n2n(){
    int idx = 0;
    ifstream fin(file_path+"Pu240_n2n.txt");
    if (!fin) {
        cerr << "Error, no such file exists" << endl;
        exit(100);
    }
    while (true) {
        fin >> Pu240n2n[idx][0];
        fin >> Pu240n2n[idx][1];
        idx++;
        if (fin.eof())
            break;
    }
    fin.close();
}

void XS_READER() {
    int left = 0, right = 71715;
    while(left <= right){
        int mid = (left + right) / 2;
        if(incident_erg >= Pu239cap[mid][0] && incident_erg <= Pu239cap[mid+1][0]){
            cap = (Pu239cap[mid][1] + Pu239cap[mid+1][1]) / 2;
            scat = (Pu239scat[mid][1] + Pu239scat[mid+1][1]) / 2;
            fis = (Pu239fi[mid][1] + Pu239fi[mid+1][1]) / 2;
            break;
        }
        if(incident_erg >= Pu239cap[mid][0]) left = mid + 1;
        else if(incident_erg < Pu239cap[mid][0]) right = mid -1;
    }
    left = 0, right = 48276;
    while(left <= right){
        int mid = (left + right) / 2;
        if(incident_erg >= Pu240cap[mid][0] && incident_erg <= Pu240cap[mid+1][0]){
            cap += (Pu240cap[mid][1] + Pu240cap[mid+1][1]) / 2;
            scat += (Pu240scat[mid][1] + Pu240scat[mid+1][1]) / 2;
            fis += (Pu240fi[mid][1] + Pu240fi[mid+1][1]) / 2;
            break;
        }
        if(incident_erg >= Pu240cap[mid][0]) left = mid + 1;
        else if(incident_erg < Pu240cap[mid][0]) right = mid -1;
    }
    n2n = 0;
    left = 0, right = 123;
    while(left <= right){
        int mid = (left + right) / 2;
        if(incident_erg >= Pu239n2n[mid][0] && incident_erg <= Pu239n2n[mid+1][0]){
            n2n += (Pu239n2n[mid][1] + Pu239n2n[mid+1][1]) / 2;
            break;
        }
        if(incident_erg >= Pu239n2n[mid][0]) left = mid + 1;
        else if(incident_erg < Pu239n2n[mid][0]) right = mid -1;
    }
    left = 0, right = 93;
    while(left <= right){
        int mid = (left + right) / 2;
        if(incident_erg >= Pu240n2n[mid][0] && incident_erg <= Pu240n2n[mid+1][0]){
            n2n += (Pu240n2n[mid][1] + Pu240n2n[mid+1][1]) / 2;
            break;
        }
        if(incident_erg >= Pu240n2n[mid][0]) left = mid + 1;
        else if(incident_erg < Pu240n2n[mid][0]) right = mid -1;
    }
}