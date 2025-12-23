// main.cpp
// Main program for ThinCurr C++ example

#include "../include/oft_capi.h"
#include "../include/error_utils.hpp"

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <cstring>
#include <fstream>

int main(int argc, char* argv[]) {
    // --- 2. Environment Initialization ---
    const char* hdf5_file_path = "tokamak_mesh_holes_16.h5";
    const char* xml_file_path = "oft_in.xml";
    const char* oftin_path = "oftcppin";

    int32_t nthreads = 28;
    int32_t slens[4];

    // Pass the user-provided config file
    log_info("run `oftpy_init`");
    oftpy_init(nthreads, oftin_path, slens, nullptr);

    // The rest of the main function remains the same...

    void* xml_ptr = nullptr;
    void* tw_obj_ptr = nullptr;
    int32_t sizes[9];
    char error_str[OFT_ERROR_SLEN] = {0};

    log_info("run `oftpy_load_xml`");
    oftpy_load_xml(xml_file_path, &xml_ptr);

    log_info("run `thincurr_setup`");
    thincurr_setup(
        hdf5_file_path, -1, nullptr, -1, nullptr, nullptr, nullptr, 0,
        &tw_obj_ptr, sizes, error_str, xml_ptr
    );
    if (check_error(error_str, "thincurr_setup")) return -1;

    int32_t nelems = sizes[7];
    int32_t n_icoils = sizes[8];
    log_info("Model setup successful. NELEMS = " + std::to_string(nelems) + ", N_ICOILS = " + std::to_string(n_icoils));

    log_info("run `thincurr_setup_io`");
    thincurr_setup_io(tw_obj_ptr, "", false, false, error_str);
    if (check_error(error_str, "thincurr_setup_io")) return -1;

    void* Mc_ptr = nullptr;
    void* Lmat_hodlr_ptr = nullptr;

    log_info("run `thincurr_Mcoil`");
    thincurr_Mcoil(tw_obj_ptr, &Mc_ptr, "", error_str);
    if (check_error(error_str, "thincurr_Mcoil")) return -1;

    log_info("run `thincurr_Lmat`");
    thincurr_Lmat(tw_obj_ptr, true, &Lmat_hodlr_ptr, "DATA_HOLDR_L.save", error_str);
    if (check_error(error_str, "thincurr_Lmat")) return -1;

    log_info("run `thincurr_Rmat`");
    thincurr_Rmat(tw_obj_ptr, false, nullptr, error_str);
    if (check_error(error_str, "thincurr_Rmat")) return -1;

    double dt = 2.0E-4;
    int32_t nsteps = 200;
    int32_t status_freq = 10;

    const int32_t n_time_points = 4;
    if (n_icoils != 2) {
        log_error("Expected 2 I-coils from model, but coil_currs is defined for 2.");
        return -1;
    }
    const double coil_currs_transposed[(1 + n_icoils) * n_time_points] = {
        0.0     ,   4.0E-3  ,   8.0E-3  ,   1.0     ,
        1.0E6   ,   1.0E6   ,   0.0     ,   0.0     ,
        0.5E6   ,   0.5E6   ,   0.5E6   ,   0.5E6
    };

    log_info("run `thincurr_td_init` for stateful stepping");
    std::vector<double> vec_ic(nelems, 0.0);
    void* td_state_ptr = nullptr;
    thincurr_td_init(
        tw_obj_ptr, false, dt, 1.0E-6, true, status_freq, 10,
        vec_ic.data(), false, &td_state_ptr, nullptr, Lmat_hodlr_ptr, error_str
    );
    if (check_error(error_str, "thincurr_td_init")) return -1;

    // Prepare coil waveforms for interpolation
    std::vector<double> times(n_time_points);
    std::array<std::vector<double>, 2> coil_wave;
    for (int i = 0; i < n_time_points; ++i) {
        times[i] = coil_currs_transposed[i];
        coil_wave[0].push_back(coil_currs_transposed[n_time_points + i]);
        coil_wave[1].push_back(coil_currs_transposed[2 * n_time_points + i]);
    }
    auto linterp = [&](int coil_idx, double t) {
        if (t <= times.front()) return coil_wave[coil_idx].front();
        if (t >= times.back()) return coil_wave[coil_idx].back();
        for (size_t k = 0; k + 1 < times.size(); ++k) {
            if (t >= times[k] && t <= times[k + 1]) {
                double w = (t - times[k]) / (times[k + 1] - times[k]);
                return (1.0 - w) * coil_wave[coil_idx][k] + w * coil_wave[coil_idx][k + 1];
            }
        }
        return coil_wave[coil_idx].back();
    };

    std::vector<double> icoil_curr(n_icoils, 0.0);
    std::vector<double> icoil_dcurr(n_icoils, 0.0);
    for (int istep = 0; istep < nsteps; ++istep) {
        double t = istep * dt;
        for (int j = 0; j < n_icoils; ++j) {
            icoil_curr[j] = linterp(j, t);
            // Match CN scheme derivative used in Fortran implementation
            double d1 = linterp(j, t + 0.25 * dt) - linterp(j, t - 0.25 * dt);
            double d2 = linterp(j, t + 1.25 * dt) - linterp(j, t + 0.75 * dt);
            icoil_dcurr[j] = d1 + d2;
        }
        thincurr_td_step(
            td_state_ptr, n_icoils, icoil_curr.data(), icoil_dcurr.data(),
            0, nullptr, false, nullptr, error_str
        );
        if (check_error(error_str, "thincurr_td_step")) return -1;
    }

    thincurr_td_finalize(td_state_ptr, vec_ic.data(), error_str);
    if (check_error(error_str, "thincurr_td_finalize")) return -1;

    log_info("Simulation Finished Successfully");

    return 0;
}