// main.cpp — ThinCurr C++ coupling example (example-22)
//
// Demonstrates the step-by-step coupling interface:
//   thincurr_td_init / thincurr_td_advance / thincurr_td_finalize
//
// The time loop, coil-current waveform storage, and linear interpolation are
// all handled in C++.  This mirrors the external-driver pattern needed when
// coupling ThinCurr to an MHD code such as M3D-C1.

#include "../include/oft_capi.h"
#include "../include/error_utils.hpp"

#include <cstdint>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Utility: piecewise-linear interpolation
// times[0..n-1] must be monotonically increasing.
// ---------------------------------------------------------------------------
static double linterp(const double* times, const double* vals, int n, double t) {
    if (t <= times[0])     return vals[0];
    if (t >= times[n - 1]) return vals[n - 1];
    for (int k = 0; k < n - 1; ++k) {
        if (times[k] <= t && t < times[k + 1]) {
            double alpha = (t - times[k]) / (times[k + 1] - times[k]);
            return vals[k] * (1.0 - alpha) + vals[k + 1] * alpha;
        }
    }
    return vals[n - 1];
}

// ---------------------------------------------------------------------------
// Coil-current waveform accessor
//
// Layout (row-major, matching the Fortran "transposed" convention):
//   row 0          : time points        [n_time]
//   row j+1 (j>=0): coil j currents    [n_time]
// ---------------------------------------------------------------------------
static double linterp_coil(const double* waveform, int n_time, int coil_j, double t) {
    const double* times = waveform;
    const double* vals  = waveform + (coil_j + 1) * n_time;
    return linterp(times, vals, n_time, t);
}

// ---------------------------------------------------------------------------
// Compute icoil_dcurr for one step using the Crank-Nicolson formula:
//
//   dcurr_j = [I(t + dt/4)  − I(t − dt/4)]
//           + [I(t + 5dt/4) − I(t + 3dt/4)]
//
// and icoil_curr_j = I(t + dt)   (coil current at the END of the step,
//                                  used only for diagnostic restart output)
// ---------------------------------------------------------------------------
static void compute_coil_drives_cn(
    const double* waveform, int n_time, int n_icoils,
    double t, double dt,
    double* icoil_dcurr, double* icoil_curr)
{
    for (int j = 0; j < n_icoils; ++j) {
        icoil_dcurr[j] =
            linterp_coil(waveform, n_time, j, t + dt / 4.0)
          - linterp_coil(waveform, n_time, j, t - dt / 4.0)
          + linterp_coil(waveform, n_time, j, t + 5.0 * dt / 4.0)
          - linterp_coil(waveform, n_time, j, t + 3.0 * dt / 4.0);
        icoil_curr[j] = linterp_coil(waveform, n_time, j, t + dt);
    }
}

// ---------------------------------------------------------------------------
// Compute icoil_dcurr for one step using Backward Euler:
//
//   dcurr_j = 2 * [I(t + 5dt/4) − I(t + 3dt/4)]
// ---------------------------------------------------------------------------
static void compute_coil_drives_be(
    const double* waveform, int n_time, int n_icoils,
    double t, double dt,
    double* icoil_dcurr, double* icoil_curr)
{
    for (int j = 0; j < n_icoils; ++j) {
        icoil_dcurr[j] = 2.0 * (
            linterp_coil(waveform, n_time, j, t + 5.0 * dt / 4.0)
          - linterp_coil(waveform, n_time, j, t + 3.0 * dt / 4.0));
        icoil_curr[j] = linterp_coil(waveform, n_time, j, t + dt);
    }
}

// ===========================================================================
// main
// ===========================================================================
int main(int argc, char* argv[]) {

    // -----------------------------------------------------------------------
    // 1. OFT initialisation
    // -----------------------------------------------------------------------
    const char* hdf5_file_path = "tokamak_mesh_holes_16.h5";
    const char* xml_file_path  = "oft_in.xml";
    const char* oftin_path     = "oftcppin";

    int32_t nthreads = 28;
    int32_t slens[4];
    char error_str[OFT_ERROR_SLEN] = {0};

    log_info("run `oftpy_init`");
    oftpy_init(nthreads, oftin_path, slens, nullptr);

    void*   xml_ptr    = nullptr;
    void*   tw_obj_ptr = nullptr;
    int32_t sizes[9];

    log_info("run `oftpy_load_xml`");
    oftpy_load_xml(xml_file_path, &xml_ptr);

    log_info("run `thincurr_setup`");
    thincurr_setup(
        hdf5_file_path, -1, nullptr, -1, nullptr, nullptr, nullptr, 0,
        &tw_obj_ptr, sizes, error_str, xml_ptr
    );
    if (check_error(error_str, "thincurr_setup")) return -1;

    int32_t nelems   = sizes[7];
    int32_t n_icoils = sizes[8];
    log_info("Model setup successful. NELEMS = " + std::to_string(nelems)
             + ", N_ICOILS = " + std::to_string(n_icoils));

    log_info("run `thincurr_setup_io`");
    thincurr_setup_io(tw_obj_ptr, "", false, false, error_str);
    if (check_error(error_str, "thincurr_setup_io")) return -1;

    void* Mc_ptr          = nullptr;
    void* Lmat_hodlr_ptr  = nullptr;

    log_info("run `thincurr_Mcoil`");
    thincurr_Mcoil(tw_obj_ptr, &Mc_ptr, "", error_str);
    if (check_error(error_str, "thincurr_Mcoil")) return -1;

    log_info("run `thincurr_Lmat`");
    thincurr_Lmat(tw_obj_ptr, true, &Lmat_hodlr_ptr, "DATA_HOLDR_L.save", error_str);
    if (check_error(error_str, "thincurr_Lmat")) return -1;

    log_info("run `thincurr_Rmat`");
    thincurr_Rmat(tw_obj_ptr, false, nullptr, error_str);
    if (check_error(error_str, "thincurr_Rmat")) return -1;

    // -----------------------------------------------------------------------
    // 2. Coil-current waveform
    //
    //   Row 0   : time points (s)
    //   Row 1   : coil-0 current (A)   — ohmic, ramps down after 8 ms
    //   Row 2   : coil-1 current (A)   — EFC,  held constant
    // -----------------------------------------------------------------------
    if (n_icoils != 2) {
        log_error("Expected 2 I-coils from model.");
        return -1;
    }
    const int32_t n_time_points     = 4;
    const int32_t n_icoils_expected = 2;  // compile-time constant matching model check
    // Storage layout: [n_time_points × (1 + n_icoils_expected)], row-major.
    // Row j corresponds to coil j−1 (row 0 = time axis).
    const double coil_waveform[(1 + n_icoils_expected) * n_time_points] = {
        0.0,    4.0e-3, 8.0e-3, 1.0,   // time  [s]
        1.0e6,  1.0e6,  0.0,    0.0,   // coil 0 [A]
        0.5e6,  0.5e6,  0.5e6,  0.5e6  // coil 1 [A]
    };

    // -----------------------------------------------------------------------
    // 3. Time-stepping parameters
    // -----------------------------------------------------------------------
    const double  dt          = 2.0e-4;  // [s]
    const int32_t nsteps      = 200;
    const bool    use_cn      = true;    // Crank-Nicolson
    const int32_t status_freq = 10;
    const int32_t plot_freq   = 10;

    // -----------------------------------------------------------------------
    // 4. Initialise ThinCurr state
    // -----------------------------------------------------------------------
    void* state_ptr = nullptr;
    std::vector<double> vec_ic(nelems, 0.0);

    log_info("run `thincurr_td_init`");
    thincurr_td_init(
        tw_obj_ptr, /*direct=*/false, dt, /*cg_atol=*/1.0e-6, /*cg_rtol=*/0.0,
        use_cn, status_freq, plot_freq,
        vec_ic.data(), /*sensor_ptr=*/nullptr,
        Lmat_hodlr_ptr, &state_ptr, error_str
    );
    if (check_error(error_str, "thincurr_td_init")) return -1;

    // -----------------------------------------------------------------------
    // 5. Time loop — coil current interpolation done here in C++
    //    (this is where an MHD code would supply boundary conditions)
    // -----------------------------------------------------------------------
    std::vector<double> icoil_dcurr(n_icoils, 0.0);
    std::vector<double> icoil_curr (n_icoils, 0.0);

    double t = 0.0;
    for (int32_t step = 0; step < nsteps; ++step) {

        // Compute drives for this step from the waveform.
        // Replace this block with MHD-provided boundary conditions when coupling.
        if (use_cn) {
            compute_coil_drives_cn(
                coil_waveform, n_time_points, n_icoils,
                t, dt, icoil_dcurr.data(), icoil_curr.data()
            );
        } else {
            compute_coil_drives_be(
                coil_waveform, n_time_points, n_icoils,
                t, dt, icoil_dcurr.data(), icoil_curr.data()
            );
        }

        thincurr_td_advance(
            state_ptr, tw_obj_ptr,
            n_icoils, icoil_dcurr.data(), icoil_curr.data(),
            /*n_volt=*/0, /*pcoil_volt=*/nullptr,
            /*sensor_ptr=*/nullptr, Lmat_hodlr_ptr, error_str
        );
        if (check_error(error_str, "thincurr_td_advance")) return -1;

        t += dt;
    }

    // -----------------------------------------------------------------------
    // 6. Finalise: retrieve solution, close history files, free state
    // -----------------------------------------------------------------------
    std::vector<double> vec_final(nelems, 0.0);

    log_info("run `thincurr_td_finalize`");
    thincurr_td_finalize(
        state_ptr, tw_obj_ptr, vec_final.data(),
        /*sensor_ptr=*/nullptr, error_str
    );
    if (check_error(error_str, "thincurr_td_finalize")) return -1;

    log_info("Simulation Finished Successfully");
    return 0;
}
