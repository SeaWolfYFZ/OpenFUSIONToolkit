// oft_capi.h
// OpenFUSIONToolkit C API declarations
// This file contains all extern "C" function declarations exposed by Fortran BIND(C)

#ifndef OFT_CAPI_H
#define OFT_CAPI_H

#include <cstdint>

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

// Time domain simulation
void thincurr_time_domain(
    void* tw_ptr, bool direct, double dt, int32_t nsteps, double lin_tol,
    bool timestep_cn, int32_t status_freq, int32_t plot_freq, void* vec_ic,
    void* sensor_ptr, int32_t ncurr, const double* curr_ptr, int32_t nvolt,
    const double* volt_ptr, bool volts_full, void* sensor_vals_ptr,
    void* hodlr_ptr, char* error_str
);

// Stateful time-domain stepping (external loop)
void thincurr_td_init(
    void* tw_ptr, bool direct, double dt, double lin_tol, bool timestep_cn,
    int32_t status_freq, int32_t plot_freq, void* vec_ic, bool volt_full,
    void** td_state_ptr, void* sensor_ptr, void* hodlr_ptr, char* error_str
);

void thincurr_td_step(
    void* td_state_ptr, int32_t ncurr, const double* icoil_curr,
    const double* icoil_dcurr, int32_t nvolt, const double* pcoil_volt,
    bool volt_full, void* sensor_vals_ptr, char* error_str
);

void thincurr_td_finalize(
    void* td_state_ptr, void* vec_out, char* error_str
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