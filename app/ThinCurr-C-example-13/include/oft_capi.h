// oft_capi.h
// OpenFUSIONToolkit C API declarations for the ThinCurr demo.
//
// IMPORTANT:
// Keep this header in sync with the Fortran `BIND(C)` interfaces in
// `src/python/wrappers/thincurr_f.F90` and `src/python/wrappers/oft_base_f.F90`.

#ifndef OFT_CAPI_H
#define OFT_CAPI_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Constants matching Fortran module definitions
#define OFT_PATH_SLEN 1024
#define OFT_ERROR_SLEN 512
#define OFT_SLEN 256

// Initialization functions
void oftpy_init(int32_t nthreads, const char* ifile, int32_t* slens, void* abort_callback);
void oftpy_load_xml(const char* xml_file, void** oft_node_ptr);

// ThinCurr setup functions
void thincurr_setup(
    const char* mesh_file, int32_t np, void* r_loc, int32_t nc, void* lc_loc,
    void* reg_loc, void* pmap_loc, int32_t jumper_start_in, void** tw_ptr,
    int32_t* sizes, char* error_str, void* xml_ptr
);

void thincurr_setup_io(void* tw_ptr, const char* basepath, bool save_debug,
                       bool legacy_hdf5, char* error_str);

// Matrix computation functions
void thincurr_Mcoil(void* tw_ptr, void** Mc_ptr, const char* cache_file,
                    char* error_str);

void thincurr_Lmat(void* tw_ptr, bool use_hodlr, void** Lmat_ptr,
                   const char* cache_file, char* error_str);

void thincurr_Rmat(void* tw_ptr, bool copy_out, void* Rmat, char* error_str);

// Time domain simulation (batch interface)
void thincurr_time_domain(
    void* tw_ptr, bool direct, double dt, int32_t nsteps, double cg_atol,
    double cg_rtol, bool timestep_cn, int32_t status_freq, int32_t plot_freq,
    void* vec_ic, void* sensor_ptr, int32_t ncurr, const double* curr_ptr,
    int32_t nvolt, const double* volt_ptr, bool volts_full, void* sensor_vals_ptr,
    void* hodlr_ptr, char* error_str
);

// Time domain simulation (wrapper using init/step/finalize)
void thincurr_time_domain_2(
    void* tw_ptr, bool direct, double dt, int32_t nsteps, double cg_atol,
    double cg_rtol, bool timestep_cn, int32_t status_freq, int32_t plot_freq,
    void* vec_ic, void* sensor_ptr, int32_t ncurr, const double* curr_ptr,
    int32_t nvolt, const double* volt_ptr, bool volts_full, void* sensor_vals_ptr,
    void* hodlr_ptr, char* error_str
);

// Step-by-step time-domain simulation interface
void thincurr_td_init(
    void* tw_ptr, void** state_ptr, bool direct, double dt, int32_t nsteps,
    double cg_atol, double cg_rtol, bool timestep_cn, int32_t status_freq,
    int32_t plot_freq, void* vec_ic, void* sensor_ptr, int32_t ncurr,
    const double* curr_ptr, int32_t nvolt, const double* volt_ptr, bool volts_full,
    void* sensor_vals_ptr, void* hodlr_ptr, char* error_str
);

void thincurr_td_step(
    void* tw_ptr, void* state_ptr, const double* icoil_curr, const double* icoil_dcurr,
    const double* pcoil_volt, int32_t nstatus, int32_t nplot, double* t_out,
    double* sol_norm_out, int32_t* nits_out, char* error_str
);

void thincurr_td_finalize(
    void* tw_ptr, void* state_ptr, double* vec_out, char* error_str
);

// ThinCurr coupling function (from thincurr_coupling executable)
void thincurr_coupling(
    void* tw_ptr, bool use_hodlr, void** Mc_ptr, const char* cache_file,
    char* error_str
);

#ifdef __cplusplus
}
#endif

#endif // OFT_CAPI_H

