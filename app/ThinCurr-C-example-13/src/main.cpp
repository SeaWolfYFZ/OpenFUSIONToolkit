// main.cpp
// Main program for ThinCurr C++ example

#include "../include/oft_capi.h"
#include "../include/error_utils.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <array>
#include <algorithm>

int main(int argc, char* argv[]) {
    const char* hdf5_file_path = "tokamak_mesh_holes_16.h5";
    const char* xml_file_path = "oft_in.xml";
    const char* oftin_path = "oftcppin";

    int32_t nthreads = 28;
    int32_t slens[4];

    log_info("run `oftpy_init`");
    oftpy_init(nthreads, oftin_path, slens, nullptr);

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
    const double coil_currs_transposed[12] = {
        0.0     ,   4.0E-3  ,   8.0E-3  ,   1.0     ,
        1.0E6   ,   1.0E6   ,   0.0     ,   0.0     ,
        0.5E6   ,   0.5E6   ,   0.5E6   ,   0.5E6
    };

    // Piecewise-linear current waveform sampled by the original batch interface.
    const std::array<double, 4> time_points = {0.0, 4.0E-3, 8.0E-3, 1.0};
    const std::array<std::array<double, 4>, 2> coil_currents = {{
        {1.0E6, 1.0E6, 0.0, 0.0},
        {0.5E6, 0.5E6, 0.5E6, 0.5E6}
    }};

    auto sample_current = [&](int coil_idx, double t) -> double {
        if (t <= time_points.front()) {
            return coil_currents[coil_idx][0];
        }
        if (t >= time_points.back()) {
            return coil_currents[coil_idx].back();
        }
        for (std::size_t i = 0; i + 1 < time_points.size(); ++i) {
            if (t <= time_points[i + 1]) {
                const double t0 = time_points[i];
                const double t1 = time_points[i + 1];
                const double y0 = coil_currents[coil_idx][i];
                const double y1 = coil_currents[coil_idx][i + 1];
                const double alpha = (t - t0) / (t1 - t0);
                return y0 + alpha * (y1 - y0);
            }
        }
        return coil_currents[coil_idx].back();
    };

    log_info("run `thincurr_td_init`");
    std::vector<double> vec_ic(nelems, 0.0);
    void* td_state_ptr = nullptr;
    thincurr_td_init(
        tw_obj_ptr, &td_state_ptr, false, dt, nsteps, 1.0E-6, 1.0E-6, true, status_freq, 10,
        vec_ic.data(), nullptr, n_time_points, coil_currs_transposed, 0, nullptr,
        false, nullptr, Lmat_hodlr_ptr, error_str
    );
    if (check_error(error_str, "thincurr_td_init")) return -1;

    log_info("run `thincurr_td_step` loop");
    std::vector<double> icoil_curr(n_icoils, 0.0);
    std::vector<double> icoil_dcurr(n_icoils, 0.0);
    double t_out = 0.0;
    double sol_norm = 0.0;
    int32_t nits = 0;

    for (int32_t istep = 0; istep < nsteps; ++istep) {
        const double t_now = istep * dt;
        for (int32_t coil = 0; coil < n_icoils; ++coil) {
            icoil_curr[coil] = sample_current(coil, t_now + dt);
            icoil_dcurr[coil] =
                sample_current(coil, t_now + dt / 4.0) -
                sample_current(coil, t_now - dt / 4.0) +
                sample_current(coil, t_now + 5.0 * dt / 4.0) -
                sample_current(coil, t_now + 3.0 * dt / 4.0);
        }
        thincurr_td_step(
            tw_obj_ptr, td_state_ptr, icoil_curr.data(), icoil_dcurr.data(), nullptr,
            status_freq, 10, &t_out, &sol_norm, &nits, error_str
        );
        if (check_error(error_str, "thincurr_td_step")) return -1;
    }

    log_info("run `thincurr_td_finalize`");
    thincurr_td_finalize(tw_obj_ptr, td_state_ptr, vec_ic.data(), error_str);
    if (check_error(error_str, "thincurr_td_finalize")) return -1;

    log_info("Simulation Finished Successfully");

    return 0;
}
