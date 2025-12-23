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
#include <cmath>

// Linear interpolation function (replaces Fortran linterp)
double linterp(const double* x_arr, const double* y_arr, int32_t n, double x) {
    if (n == 1) return y_arr[0];
    if (x <= x_arr[0]) return y_arr[0];
    if (x >= x_arr[n-1]) return y_arr[n-1];

    // Find interval
    int32_t i = 0;
    while (i < n - 1 && x > x_arr[i + 1]) i++;

    // Interpolate
    double t = (x - x_arr[i]) / (x_arr[i + 1] - x_arr[i]);
    return y_arr[i] + t * (y_arr[i + 1] - y_arr[i]);
}

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
    // Coil waveform: [time, coil1, coil2] transposed for column-major access
    // Time points: 0.0, 4ms, 8ms, 1.0s
    // Coil 1: 1.0E6, 1.0E6, 0.0, 0.0
    // Coil 2: 0.5E6, 0.5E6, 0.5E6, 0.5E6
    const double coil_currs_transposed[(1 + n_icoils) * n_time_points] = {
        0.0     ,   4.0E-3  ,   8.0E-3  ,   1.0     ,  // time
        1.0E6   ,   1.0E6   ,   0.0     ,   0.0     ,  // coil 1
        0.5E6   ,   0.5E6   ,   0.5E6   ,   0.5E6     // coil 2
    };

    // =========================================================================
    // NEW: External time domain control with init/advance/finalize
    // =========================================================================
    log_info("=== Using NEW external time domain control ===");

    std::vector<double> vec_ic(nelems, 0.0);
    void* state_ptr = nullptr;

    // Initialize simulation
    log_info("run `thincurr_td_init`");
    thincurr_td_init(
        tw_obj_ptr, false, dt, 1.0E-6, true, status_freq, 10,
        vec_ic.data(), nullptr, Lmat_hodlr_ptr, &state_ptr, error_str
    );
    if (check_error(error_str, "thincurr_td_init")) return -1;

    log_info("thincurr_td_init completed successfully, state_ptr = " + std::to_string(reinterpret_cast<uintptr_t>(state_ptr)));

    // Time loop with external coil current control
    std::vector<double> icoil_curr(n_icoils);
    std::vector<double> icoil_dcurr(n_icoils);

    for (int32_t i = 0; i < nsteps; i++) {
        if (i == 0) log_info("Starting time step 0...");
        double t = i * dt;

        double freq = 1.0E3/16.0;
        double I1_amplitude = 1.0E6;
        double I2_amplitude = 0.5E6;

        double omega = 2.0 * M_PI * freq;

        // Interpolate coil currents at this time step
        for (int32_t j = 0; j < n_icoils; j++) {
            icoil_curr[j] = linterp(coil_currs_transposed,
                                    &coil_currs_transposed[(j + 1) * n_time_points],
                                    n_time_points, t);
        }

        // // Compute dI/dt for Crank-Nicolson (using 4-point stencil as in original code)
        // if (i == 0) {
        //     // First step: use forward difference
        //     double t_next = (i + 1) * dt;
        //     for (int32_t j = 0; j < n_icoils; j++) {
        //         double curr_next = linterp(coil_currs_transposed,
        //                                   &coil_currs_transposed[(j + 1) * n_time_points],
        //                                   n_time_points, t_next);
        //         icoil_dcurr[j] = (curr_next - icoil_curr[j]) / dt;
        //     }
        // } else {
        //     // Crank-Nicolson: average of derivatives at start and end of step
        //     for (int32_t j = 0; j < n_icoils; j++) {
        //         double d1 = linterp(coil_currs_transposed,
        //                            &coil_currs_transposed[(j + 1) * n_time_points],
        //                            n_time_points, t + dt/4.0);
        //         double d2 = linterp(coil_currs_transposed,
        //                            &coil_currs_transposed[(j + 1) * n_time_points],
        //                            n_time_points, t - dt/4.0);
        //         double d3 = linterp(coil_currs_transposed,
        //                            &coil_currs_transposed[(j + 1) * n_time_points],
        //                            n_time_points, t + dt*5.0/4.0);
        //         double d4 = linterp(coil_currs_transposed,
        //                            &coil_currs_transposed[(j + 1) * n_time_points],
        //                            n_time_points, t + dt*3.0/4.0);
        //         icoil_dcurr[j] = (d1 - d2) + (d3 - d4);
        //     }
        // }

        icoil_curr[0]   = I1_amplitude * sin(omega * t);
        icoil_dcurr[0]  = omega * I1_amplitude * cos(omega * t) * dt;
        icoil_curr[1]   = I2_amplitude;
        icoil_dcurr[1]  = 0.0;

        // Advance one time step
        bool save_plot = ((i + 1) % 10 == 0);
        if (i == 0) log_info("Calling thincurr_td_advance for step 0...");
        thincurr_td_advance(
            state_ptr, icoil_curr.data(), icoil_dcurr.data(), nullptr,
            save_plot, error_str
        );
        if (check_error(error_str, "thincurr_td_advance at step " + std::to_string(i))) return -1;
    }

    // Finalize simulation
    log_info("run `thincurr_td_finalize`");
    thincurr_td_finalize(state_ptr, vec_ic.data(), error_str);
    if (check_error(error_str, "thincurr_td_finalize")) return -1;

    log_info("=== NEW external control simulation Finished Successfully ===");

    return 0;
}