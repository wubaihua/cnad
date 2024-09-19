

#include <complex.h>
#include <stdint.h>
#include "def.h"
#include "constant.h"
#include "gmath.h"
#include <stdio.h>
#include <math.h>
#include "cJSON.h"
#include "msmodel.h"
#include "msmodelio.h"
#include <stdbool.h>
#include "def_host.h"
#include "iofun_slave.h"
#include <slave.h>


// int mpi_size, mpi_rank, mpi_ierr;

// Dynamically allocated arrays
// double complex *mpi_den; // 3D complex array: [size1][size2][size3]
// double complex *mpi_cfall; // 5D complex array: [size1][size2][size3][size4][size5]
// double complex *mpi_cfeff; // 1D complex array: [size1]
// double *real_rho; // 3D double array: [size1][size2][size3]
// double *imag_rho; // 3D double array: [size1][size2][size3]
// double *real_cfall; // 5D double array: [size1][size2][size3][size4][size5]
// double *imag_cfall; // 5D double array: [size1][size2][size3][size4][size5]
// double *real_cfeff; // 1D double array: [size1]
// double *imag_cfeff; // 1D double array: [size1]
// double *mpi_population; // 2D double array: [size1][size2]
// double *mpi_pop_fb; // 3D double array: [size1][size2][size3]
// double *mpi_R_nuc_mean; // 3D double array: [size1][size2][size3]
// double *mpi_P_nuc_mean; // 3D double array: [size1][size2][size3]
// double *mpi_R2_nuc_mean; // 3D double array: [size1][size2][size3]
// double *mpi_P2_nuc_mean; // 3D double array: [size1][size2][size3]
// unsigned long long *mpi_N_nan_sum; // 1D int array: [size1]
// double *mpi_real_den;
// double *mpi_imag_den;
// double *mpi_real_cfeff;
// double *mpi_imag_cfeff;


int *count_st; // 2D int array: [size1][size2]
int *mpi_count_st; // 2D int array: [size1][size2]
int *count_pertraj; // 1D int array: [size1]

double *R_nuc_mean; // 3D double array: [size1][size2][size3]
double *P_nuc_mean; // 3D double array: [size1][size2][size3]
double *R2_nuc_mean; // 3D double array: [size1][size2][size3]
double *P2_nuc_mean; // 3D double array: [size1][size2][size3]

// char *filepath; // Path of Model input file
// char *workpath;

int Ndof1, Ndof2;
int Nstate;
int init_occ;

double *xe, *pe; // 1D double arrays: [size1] Meyer-Miller mapping variables
double *ye, *pxe, *pye; // 1D double arrays: [size1] Li-Miller mapping variables
double complex *ce; // 1D complex array: [size1] Electronic amplitude
double complex *correfun_t; // 2D complex array: [size1][size2] Electronic reduced density matrix for evolution
double complex *gamma_cv; // 2D complex array: [size1][size2] Commutator matrix or zero-point energy parameters
double complex *den_e, *den_e4nuc; // 2D complex arrays: [size1][size2]

double *xe_mb, *pe_mb; // 2D double arrays: [size1][size2]
int *state_unsmash; // 1D int array: [size1]
// int if_typemb;
double *gc_unsmash, *gp_unsmash; // 1D double arrays: [size1]

double *xe_cv, *pe_cv; // 1D double arrays: [size1]

// More variables...

// int type_evo;
// double gamma_zpe;
// double sigma2_lsc;
// double gamma1_lsc, gamma2_lsc;
double identity_M1, identity_M2;
// int scheme_cvk, scheme_focus;
// double alpha_gdtwa, beta_gdtwa, delta_gdtwa, eps_gdtwa;
// int if_alpha;
int id_state, id_hop;
int id_state_old;
double *t_decoh, *t_coh, *L_red; // 1D double arrays: [size1]
// double w_dish;

double *xe_old, *pe_old; // 1D double arrays: [size1]
double *P_nuc_old, *R_nuc_old; // 2D double arrays: [size1][size2]
double *force_old; // 2D double array: [size1][size2]
double complex *gamma_cv_old, *den_e_old; // 2D complex arrays: [size1][size2]
double *nacv_old, *V_old, *dV_old, *E_adia_old, *dv_adia_old; // Various dimensional arrays


// 动态分配的数组声明
double *nac_check_old; // 4D double array: [size1][size2][size3][size4]
double *U_ref_old; // 2D double array: [size1][size2]

double *P_nuc_old_traj; // 2D double array: [size1][size2]
double *R_nuc_old_traj; // 2D double array: [size1][size2]

double *deltaR_afssh; // 4D double array: [size1][size2][size3][size4]
double *deltaP_afssh; // 4D double array: [size1][size2][size3][size4]
double *deltaF_afssh; // 4D double array: [size1][size2][size3][size4]

// int index_t0;
// int index_t0_1, index_t0_2;
double complex correfun_0;

double complex *correfun_0_ms2; // 1D complex array: [size1]

double complex *cf0; // 2D complex array: [size1][size2]
double complex *cfall; // 5D complex array: [size1][size2][size3][size4][size5]
double complex *cfeff; // 1D complex array: [size1]
double *weight0; // 2D double array: [size1][size2]
double *weightt; // 2D double array: [size1][size2]
// int if_allcf, allcf_times;

// double complex correfun_0_pldm1, correfun_0_pldm2;
double complex *prop_pldm; // 2D complex array: [size1][size2]

double *U_d2a; // 2D double array: [size1][size2]
double *E_adia; // 1D double array: [size1]
double *nac; // 4D double array: [size1][size2][size3][size4]
double *dv_adia; // 4D double array: [size1][size2][size3][size4]
double *P_kin; // 2D double array: [size1][size2]
double *nac_check; // 4D double array: [size1][size2][size3][size4]
double *U_ref; // 2D double array: [size1][size2]
double *overlap_adia; // 2D double array: [size1][size2]
// int rep; // 0 for diabatic, 1 for adiabatic

double *U_d2a_old; // 2D double array: [size1][size2]

// int ifcv; // 0: no adjustment; -1: adjustment without evolution; 1: cv adjustment 
// int ifid;
double complex *den; // 3D complex array: [size1][size2][size3] (total electronic reduced density matrix)
double *population; // 2D double array: [size1][size2]
double *pop_fb; // 3D double array: [size1][size2][size3]

double complex *den_traj; // 3D complex array: [size1][size2][size3]
double complex *den2_traj; // 3D complex array: [size1][size2][size3]
double *population_traj; // 2D double array: [size1][size2]
double *population2_traj; // 2D double array: [size1][size2]
double *pop_fb_traj; // 3D double array: [size1][size2][size3]
double *pop_fb2_traj; // 3D double array: [size1][size2][size3]

// int if_st_fb;

double *V; // 2D double array: [size1][size2]
double *dV; // 4D double array: [size1][size2][size3][size4]
double *V_ref; // 3D double array: [size1][size2][size3]
double *dV_ref; // 5D double array: [size1][size2][size3][size4][size5]
double complex *propagator; // 2D complex array: [size1][size2]
double complex *propagator_path; // 2D complex array: [size1][size2]
double complex *propagator_switch; // 2D complex array: [size1][size2]

double *R_nuc; // 2D double array: [size1][size2]
double *P_nuc; // 2D double array: [size1][size2]
double *mass; // 2D double array: [size1][size2]

double *force; // 2D double array: [size1][size2]
double *force_nuc; // 2D double array: [size1][size2]
double *force_ref; // 3D double array: [size1][size2][size3]
double *force_nuc_ref; // 3D double array: [size1][size2][size3]
// int type_traj_sed;

// double beta, temperature;

// char method[20];

// double dt, ttot,
double t_now;
double *timegrid; // 1D double array: [size1]
// long long Nbreak, Ngrid;
// long long Ntraj;

// char unit_t[20];
// double unittrans_t;

// int outputtype;

// int calforcetype;

// int ifoutputmpi;

// int sampletype;

// int if_st_nan;
unsigned long long *N_nan_sum; // 1D int array: [size1]

// int if_traj;

// int type_phase;

// int type_ad_fssh;

// int if_ref;

// int if_1st;

// int if_inv_focus;

// int if_Pdis, s_N;
// double s_start, s_end;
double *s; // 1D double array: [size1]
double *real_expisP; // 1D double array: [size1]
double *imag_expisP; // 1D double array: [size1]
double complex *expisP; // 1D complex array: [size1]
double complex *mpi_expisP; // 1D complex array: [size1]
double *mpi_real_expisP;
double *mpi_imag_expisP;

double *R_nuc_ref; // 3D double array: [size1][size2][size3]
double *P_nuc_ref; // 3D double array: [size1][size2][size3]
double *xe_ref; // 2D double array: [size1][size2]
double *pe_ref; // 2D double array: [size1][size2]
double *ye_ref; // 2D double array: [size1][size2]
double *pxe_ref; // 2D double array: [size1][size2]
double *pye_ref; // 2D double array: [size1][size2]
double complex *gamma_cv_ref; // 3D complex array: [size1][size2][size3]
double complex *den_e_ref; // 3D complex array: [size1][size2][size3]
double complex *inverse_kernel; // 2D complex array: [size1][size2]

// int Nref;

double *pop0; // 1D double array: [size1]

// bool if_ad_nac;

// int mean_nuc;

// bool if_occ_rand;

// int if_engconsv;
double complex *engconsv_adjmat; // 2D complex array: [size1][size2]

// int if_RBC;

double measure_mash;
// int type_mash;
double complex rho0_mash[4], rhot_mash[4];

double U0_mash[4], mea_mat_mash[4];

double complex *rho0_unsmash; // 2D complex array: [size1][size2]
double complex *rhot_unsmash; // 2D complex array: [size1][size2]
double *U0_unsmash; // 2D double array: [size1][size2]
double *mea_mat_unsmash; // 2D double array: [size1][size2]

// int ifBA;
double *E_adia_old_traj; // 1D double array: [size1]
double *nac_BAeff; // 4D double array: [size1][size2][size3][size4]
double *tdc_BA; // 2D double array: [size1][size2]
double *nac_old; // 4D double array: [size1][size2][size3][size4]
double *P_nuc_BA_old; // 2D double array: [size1][size2]

int iflbd, occ_lbd;
double *rate_lbd; // 1D double array: [size1]
double rate_para_lbd;

int ifmfsh, occ_mfsh;
double thres;

int occ_rdmfocus;

int ifrw;
double rwfactor, beta_rw, temperature_rw;

// int ifmodprop;
double complex *propagator_ref; // 2D complex array: [size1][size2]

// int ifcorreden;

int memorylength;
int itime_save, i_re_save;

// int if_st_eng;
double *energy_est; // 1D double array: [size1]
double *mpi_energy_est; // 1D double array: [size1]

// int typeevo_ele;

double *A_jump; // 2D double array: [size1][size2]
double *lambda_jump; // 1D double array: [size1]
double *U_jump; // 2D double array: [size1][size2]
double *alpha_jump; // 1D double array: [size1]
int type_jump;

int if_switchcv, occ_switchcv, max_switchcv;

double Pn0_mf2, Wn0_mf2, energy_mf2;
int if_inv_evo, if_BO_mf2, id_max_mf2sh, type_eom_mf2sh, id_init_occ_mf2, if_bak_mf2cv;
double *da_mf2; // 1D double array: [size1]
double *R_nuc_old_mf2cv; // 2D double array: [size1][size2]
double *P_nuc_old_mf2cv; // 2D double array: [size1][size2]
double *xe_old_mf2cv; // 1D double array: [size1]
double *pe_old_mf2cv; // 1D double array: [size1]
double *dE_old_mf2cv; // 1D double array: [size1]
double *dE_mf2cv; // 1D double array: [size1]

// int type_prop_adia;

double complex *G_xpconfg; // 2D complex array: [size1][size2]
double *permutation_ms2; // 3D double array: [size1][size2][size3]
double complex *G_ms2; // 2D complex array: [size1][size2]
double thres_ms2;
int *index_ms2; // 1D int array: [size1]
int *index_old_ms2; // 1D int array: [size1]

double *R_nuc_state; // 2D double array: [size1][size2]
double *P_nuc_state; // 2D double array: [size1][size2]
double *V_state; // 2D double array: [size1][size2]
double *dv_state; // 4D double array: [size1][size2][size3][size4]
double *dv_adia_state; // 4D double array: [size1][size2][size3][size4]
double *U_d2a_state; // 2D double array: [size1][size2]
double *U_ref_state; // 2D double array: [size1][size2]
double *force_nuc_state; // 2D double array: [size1][size2]
int if_statetraj;

bool ifBC_BCMF;

int scheme_cvsh;

// int ifswitchforce, ifmashforce;

// int ifscaleenergy;

// double gamma_rescale;
// int ifscalegamma;

double E_conserve;

int ifmsbranch, type_traj_msbranch;
int itime_start_msbranch, itime_end_msbranch;
int i_re_start_msbranch, igrid_start_msbranch;
int iwrong_msbranch;
double time_start_msbranch, time_end_msbranch;
double *R_nuc_brapoint; // 2D double array: [size1][size2]
double *P_nuc_brapoint; // 2D double array: [size1][size2]
double *xe_brapoint; // 1D double array: [size1]
double *pe_brapoint; // 1D double array: [size1]
double complex *gamma_cv_brapoint; // 2D complex array: [size1][size2]
double complex *den_e_brapoint; // 2D complex array: [size1][size2]

double complex *correfun_t_oldtraj; // 3D complex array: [size1][size2][size3]
double *R_nuc_oldtraj; // 3D double array: [size1][size2][size3]
double *P_nuc_oldtraj; // 3D double array: [size1][size2][size3]

// int direc_padj;

// int ifcount;

// int ifreflp;

// int ifreflp_mash;

// int ifhardwall;

double *eig_cv; // 1D double array: [size1]
double *eig_cv_mat; // 2D double array: [size1][size2]
double complex *commu_vari; // 2D complex array: [size1][size2]

// int ifzpecorr;

// int iflangevin;
// double eta_langevin;

double scale_sqc2;

// int type_algorithm;

// int type_prop_4cont;

// int scaleenergy_type;

// int n_step_algo5;

// int allow_hop;

// int if_traceless_force;






void initial_vari() {
    int i;


    

    // 分配内存
    
    R_nuc = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    P_nuc = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    mass = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    force = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    if (forcetype == 1) {
        force_nuc = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    }
    V = (double *)malloc(Nstate * Nstate * sizeof(double));
    dV = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));

    R_nuc_old = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    P_nuc_old = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    force_old = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    V_old = (double *)malloc(Nstate * Nstate * sizeof(double));
    dV_old = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    
    propagator = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
    propagator_switch = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
    if (ifmodprop == 1) {
        propagator_ref = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
    }

    Ngrid = (long long)(ttot / dt) / Nbreak + 1;
    timegrid = (double *)malloc(Ngrid * sizeof(double));
    if (if_allcf == 0) {
        if (outputtype >= 0) {
            den = (double complex *)malloc(Nstate * Nstate * Ngrid * sizeof(double complex));
            memset(den, 0, Nstate * Nstate * Ngrid * sizeof(double complex));
        }
        if (outputtype != 0) {
            population = (double *)malloc(Nstate * Ngrid * sizeof(double));
            memset(population, 0, Nstate * Ngrid * sizeof(double));
        }
    } else if (if_allcf == 1) {
        cf0 = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
        memset(cf0, 0, Nstate * Nstate * sizeof(double complex ));
        // cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double ));
        // memset(cfall, 0, Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
        type_evo = 1;
    } else if (if_allcf >= 2) {
        cf0 = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
        memset(cf0, 0, Nstate * Nstate * sizeof(double complex ));
        cfeff = (double  complex *)malloc(Ngrid * sizeof(double complex ));
        memset(cfeff, 0, Ngrid * sizeof(double complex ));
        type_evo = 1;
        weight0 = (double *)malloc(Nstate * Nstate * sizeof(double));
        weightt = (double *)malloc(Nstate * Nstate * sizeof(double));
    }
    
    if (if_st_eng == 1) {
        energy_est = (double *)malloc(Ngrid * sizeof(double));
    }

    // if (iflbd > 0) {
    //     rate_lbd = (double *)malloc(Nstate * sizeof(double));
    //     memset(rate_lbd, 0, Nstate * sizeof(double));
    //     if (type_evo != 3) type_evo = 1;
    //     occ_lbd = init_occ;
    // }

    if (ifmfsh > 0) {
        if (type_evo != 3) type_evo = 1;
    }

    correfun_t = (double complex *)malloc(Nstate * Nstate * sizeof(double complex));
    memset(correfun_t, 0, Nstate * Nstate * sizeof(double complex));
    gamma_cv = (double complex *)malloc(Nstate * Nstate * sizeof(double complex));
    memset(gamma_cv, 0, Nstate * Nstate * sizeof(double complex));
    gamma_cv_old = (double complex *)malloc(Nstate * Nstate * sizeof(double complex));
    
    // if (ifcv == 3) {
    //     commu_vari = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     eig_cv = (double *)malloc(Nstate * sizeof(double));
    //     eig_cv_mat = (double *)malloc(Nstate * Nstate * sizeof(double));
    // }

    if (strcmp(method, "GDTWA") == 0 || strcmp(method, "eGDTWA") == 0) {
        if (type_evo != 3) type_evo = 1;
        if (if_inv_focus == 1) {
            inverse_kernel = (double  complex *)malloc(Nstate * Nstate * sizeof(double  complex ));
        }
   
    // } else if (strcmp(method, "test") == 0) {
    //     // type_evo = 1;
    // } else if (strcmp(method, "DISH") == 0 || strcmp(method, "dish") == 0) {
    //     t_decoh = (double *)malloc(Nstate * sizeof(double));
    //     t_coh = (double *)malloc(Nstate * sizeof(double));
    //     L_red = (double *)malloc(Nstate * sizeof(double));
    // } else if (strcmp(method, "A-FSSH") == 0 || strcmp(method, "a-fssh") == 0 ||
    //            strcmp(method, "AFSSH") == 0 || strcmp(method, "afssh") == 0) {
    //     deltaR_afssh = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     deltaP_afssh = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     deltaF_afssh = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     if (type_evo != 3) type_evo = 1;
    // } else if (strcmp(method, "GFSH") == 0 || strcmp(method, "gfsh") == 0 ||
    //            strcmp(method, "SC-FSSH") == 0 || strcmp(method, "sc-fssh") == 0 ||
    //            strcmp(method, "CC-FSSH") == 0 || strcmp(method, "cc-fssh") == 0) {
    //     xe_old = (double *)malloc(Nstate * sizeof(double));
    //     pe_old = (double *)malloc(Nstate * sizeof(double));
    // } else if (strcmp(method, "rdmfocus") == 0) {
    //     if (type_evo != 3) type_evo = 1;
    // } else if (strcmp(method, "focus_jump") == 0) {
    //     A_jump = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     lambda_jump = (double *)malloc(Nstate * sizeof(double));
    //     U_jump = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     alpha_jump = (double *)malloc(Nstate * sizeof(double));
    //     propagator_path = (double *)malloc(Nstate * Nstate * sizeof(double));
    // } else if (strcmp(method, "mf2-sh") == 0 || strcmp(method, "MF2-SH") == 0 ||
    //            strcmp(method, "mf2cv") == 0 || strcmp(method, "MF2CV") == 0) {
    //     da_mf2 = (double *)malloc(Nstate * sizeof(double));
    //     R_nuc_old_mf2cv = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     P_nuc_old_mf2cv = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     xe_old_mf2cv = (double *)malloc(Nstate * sizeof(double));
    //     pe_old_mf2cv = (double *)malloc(Nstate * sizeof(double));
    //     dE_old_mf2cv = (double *)malloc(Nstate * sizeof(double));
    //     dE_mf2cv = (double *)malloc(Nstate * sizeof(double));
    // } else if (strcmp(method, "ms2") == 0 || strcmp(method, "MS2") == 0 ||
    //            strcmp(method, "ms3") == 0 || strcmp(method, "MS3") == 0) {
    //     type_evo = 4;
    //     correfun_0_ms2 = (double *)malloc(Nstate * sizeof(double));
    //     permutation_ms2 = (double *)malloc(Nstate * Nstate * Nstate * sizeof(double));
    //     for (i = 0; i < Nstate; i++) {
    //         gen_permutation(Nstate, &permutation_ms2[i * Nstate * Nstate], i);
    //     }
    //             index_ms2 = (double *)malloc(Nstate * sizeof(double));
    //     index_old_ms2 = (double *)malloc(Nstate * sizeof(double));
    // } else if (strcmp(method, "mstraj") == 0 || strcmp(method, "mstraj2") == 0) {
    //     R_nuc_state = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     P_nuc_state = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     dv_state = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     V_state = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     dv_adia_state = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     U_d2a_state = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     U_ref_state = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     force_nuc_state = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    // } else if (strcmp(method, "bcmf") == 0 || strcmp(method, "BCMF") == 0) {
    //     E_adia_old_traj = (double *)malloc(Nstate * sizeof(double));
    //     xe_old = (double *)malloc(Nstate * sizeof(double));
    //     pe_old = (double *)malloc(Nstate * sizeof(double));
    //     R_nuc_old_traj = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     P_nuc_old_traj = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    // } else if (strcmp(method, "sed2") == 0 || strcmp(method, "SED2") == 0 ||
    //            strcmp(method, "sed3") == 0 || strcmp(method, "SED3") == 0) {
    //     if (outputtype >= 0) {
    //         den_traj = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
    //         memset(den_traj, 0, Nstate * Nstate * Ngrid * sizeof(double));
    //         den2_traj = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
    //         memset(den2_traj, 0, Nstate * Nstate * Ngrid * sizeof(double));
    //     }
    //     if (outputtype != 0) {
    //         population_traj = (double *)malloc(Nstate * Ngrid * sizeof(double));
    //         memset(population_traj, 0, Nstate * Ngrid * sizeof(double));
    //         population2_traj = (double *)malloc(Nstate * Ngrid * sizeof(double));
    //         memset(population2_traj, 0, Nstate * Ngrid * sizeof(double));
    //     }
    //     if (if_st_fb == 1) {
    //         pop_fb_traj = (double *)malloc(Nstate * Ngrid * 2 * sizeof(double));
    //         memset(pop_fb_traj, 0, Nstate * Ngrid * 2 * sizeof(double));
    //         pop_fb2_traj = (double *)malloc(Nstate * Ngrid * 2 * sizeof(double));
    //         memset(pop_fb2_traj, 0, Nstate * Ngrid * 2 * sizeof(double));
    //     }
    // } else if (strcmp(method, "unsmash") == 0 || strcmp(method, "UNSMASH") == 0 ||
    //            strcmp(method, "unSMASH") == 0 || strcmp(method, "unsmash-mf") == 0 ||
    //            strcmp(method, "UNSMASH-MF") == 0 || strcmp(method, "unSMASH-MF") == 0) {
    //     if_typemb = 1;
    //     type_evo = 5;
    //     xe_mb = (double *)malloc(2 * (Nstate - 1) * sizeof(double));
    //     pe_mb = (double *)malloc(2 * (Nstate - 1) * sizeof(double));
    //     state_unsmash = (double *)malloc((Nstate - 1) * sizeof(double));
    //     rho0_unsmash = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     rhot_unsmash = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     U0_unsmash = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     gp_unsmash = (double *)malloc(Nstate * sizeof(double));
    //     gc_unsmash = (double *)malloc(Nstate * sizeof(double));
    //     mea_mat_unsmash = (double *)malloc(Nstate * Nstate * sizeof(double));
    }


    // if (ifmsbranch > 0) {
    //     R_nuc_state = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     P_nuc_state = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     dv_state = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     V_state = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     dv_adia_state = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     U_d2a_state = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     U_ref_state = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     force_nuc_state = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     R_nuc_brapoint = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     P_nuc_brapoint = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
    //     xe_brapoint = (double *)malloc(Nstate * sizeof(double));
    //     pe_brapoint = (double *)malloc(Nstate * sizeof(double));
    //     gamma_cv_brapoint = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     correfun_t_oldtraj = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
    // }
    // if (ifcount == 1) {
    //     count_st = (int *)malloc(5 * Ngrid * sizeof(int));
    //     count_pertraj = (int *)malloc(5 * sizeof(int));
    //     memset(count_st, 0, 5 * Ngrid * sizeof(int));
    //     memset(count_pertraj, 0, 5 * sizeof(int));
    // }
    init_occ = init_occ_4read;
    if (init_occ == 0) {
        if_occ_rand = 1;
    }
    // if (if_ref == 1) {
    //     R_nuc_ref = (double *)malloc(Nref * Ndof1 * Ndof2 * sizeof(double));
    //     P_nuc_ref = (double *)malloc(Nref * Ndof1 * Ndof2 * sizeof(double));
    //     dV_ref = (double *)malloc(Nref * Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
    //     V_ref = (double *)malloc(Nref * Nstate * Nstate * sizeof(double));
    //     force_ref = (double *)malloc(Nref * Ndof1 * Ndof2 * sizeof(double));
    //     force_nuc_ref = (double *)malloc(Nref * Ndof1 * Ndof2 * sizeof(double));
    //     xe_ref = (double *)malloc(Nref * Nstate * sizeof(double));
    //     pe_ref = (double *)malloc(Nref * Nstate * sizeof(double));
    //     den_e_ref = (double *)malloc(Nref * Nstate * Nstate * sizeof(double));
    //     gamma_cv_ref = (double *)malloc(Nref * Nstate * Nstate * sizeof(double));
    //     prop_pldm = (double *)malloc(Nstate * Nstate * sizeof(double));
    //     memset(prop_pldm, 0, Nstate * Nstate * sizeof(double));
    //     for (i = 0; i < Nstate; i++) {
    //         prop_pldm[i * Nstate + i] = 1;
    //     }
    //     if (type_evo != 3) type_evo = 1;
    // }
    switch (type_evo) {
        case 0:
        case 5:
            xe = (double *)malloc(Nstate * sizeof(double));
            pe = (double *)malloc(Nstate * sizeof(double));
            xe_old = (double *)malloc(Nstate * sizeof(double));
            pe_old = (double *)malloc(Nstate * sizeof(double));
            break;
        case 1:
        case 2:
            xe = (double *)malloc(Nstate * sizeof(double));
            pe = (double *)malloc(Nstate * sizeof(double));
            den_e = (double  complex *)malloc(Nstate * Nstate * sizeof(double  complex ));
            xe_old = (double *)malloc(Nstate * sizeof(double));
            pe_old = (double *)malloc(Nstate * sizeof(double));
            den_e_old = (double  complex *)malloc(Nstate * Nstate * sizeof(double  complex ));
            break;
        // case 3:
        //     xe = (double *)malloc(Nstate * sizeof(double));
        //     pe = (double *)malloc(Nstate * sizeof(double));
        //     den_e4nuc = (double *)malloc(Nstate * Nstate * sizeof(double));
        //     den_e = (double *)malloc(Nstate * Nstate * sizeof(double));
        //     break;
        // case 4:
        //     xe = (double *)malloc(Nstate * sizeof(double));
        //     pe = (double *)malloc(Nstate * sizeof(double));
        //     G_xpconfg = (double *)malloc(Nstate * Nstate * sizeof(double));
        //     break;
    }
    
    xe_cv = (double *)malloc(Nstate * sizeof(double));
    pe_cv = (double *)malloc(Nstate * sizeof(double));
    if (sampletype == 2) {
        rep = 1;
    }
    if (rep == 1 || sampletype == 3) {
        U_d2a = (double *)malloc(Nstate * Nstate * sizeof(double));
        U_ref = (double *)malloc(Nstate * Nstate * sizeof(double));
        U_d2a_old = (double *)malloc(Nstate * Nstate * sizeof(double));
        U_ref_old = (double *)malloc(Nstate * Nstate * sizeof(double));
        memset(U_ref, 0, Nstate * Nstate * sizeof(double));
        E_adia = (double *)malloc(Nstate * sizeof(double));
        dv_adia = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        nac = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        E_adia_old = (double *)malloc(Nstate * sizeof(double));
        dv_adia_old = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        nac_old = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        P_kin = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
        nac_check = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        nac_check_old = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        if (type_prop_adia > 0) {
            overlap_adia = (double *)malloc(Nstate * Nstate * sizeof(double));
        }
        memset(nac_check, 0, Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        // if (ifBA == 1) {
        //     E_adia_old_traj = (double *)malloc(Nstate * sizeof(double));
        //     nac_BAeff = (double *)malloc(Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));
        //     tdc_BA = (double *)malloc(Nstate * Nstate * sizeof(double));
        //     P_nuc_BA_old = (double *)malloc(Ndof1 * Ndof2 * sizeof(double));
        //     memset(tdc_BA, 0, Nstate * Nstate * sizeof(double));
        // }
    }
    
    if (if_st_fb == 1) {
        pop_fb = (double *)malloc(Nstate * Ngrid * 2 * sizeof(double));
        memset(pop_fb, 0, Nstate * Ngrid * 2 * sizeof(double));
    }
    if (if_Pdis == 1) {
        s = (double *)malloc(s_N * sizeof(double));
        expisP = (double  complex *)malloc(s_N * sizeof(double complex ));
        memset(expisP, 0, s_N * sizeof(double complex ));
        for (i = 0; i < s_N; i++) {
            s[i] = s_start + i * (s_end - s_start) / s_N;
        }
    }
   
    if (temperature != 0.0 && beta == 0.0) {
        beta = 1 / (kb * temperature);
        if (temperature < 1.0) beta = 10000000;
    }
    
    if (ifrw > 0) {
        if (beta_rw < 0) {
            beta_rw = 1.0 / (kb * temperature_rw);
            if (temperature_rw < 1.0) beta_rw = 10000000;
        }
    }
    if (strcmp(unit_t, "au") == 0) {
        unittrans_t = 1.0;
    } else if (strcmp(unit_t, "fs") == 0) {
        unittrans_t = 1.0 / au_2_fs;
    }
    
    ttot *= unittrans_t;
    dt *= unittrans_t;

    if (mean_nuc == 1) {
        R_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        P_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        R2_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        P2_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memset(R_nuc_mean, 0, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memset(P_nuc_mean, 0, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memset(R2_nuc_mean, 0, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memset(P2_nuc_mean, 0, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        if (ifmsbranch > 0) {
            R_nuc_oldtraj = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
            P_nuc_oldtraj = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        }
    }

    // printf("111166\n");
    // if (if_engconsv == 1) {
    //     if (type_evo != 3) type_evo = 1;
    //     engconsv_adjmat = (double *)malloc(Nstate * Nstate * sizeof(double));
    // }
    // printf("111166\n");
    //  printf("Ngrid=%d\n",Ngrid);
    // unsigned long long Ntemp=Ngrid;
    N_nan_sum = (unsigned long long *)malloc(Ngrid * sizeof(unsigned long long));
    // printf("22266\n");
    mpi_N_nan_sum = (unsigned long long *)malloc(Ngrid * sizeof(unsigned long long));
    // printf("333366\n");
    memset(N_nan_sum, 0, Ngrid * sizeof(unsigned long long));
    // printf("444466\n");
    if (ifswitchforce > 0) {
        if (rep == 0) {
            U_d2a = (double *)malloc(Nstate * Nstate * sizeof(double));
            E_adia = (double *)malloc(Nstate * sizeof(double));
            U_d2a_old = (double *)malloc(Nstate * Nstate * sizeof(double));
            E_adia_old = (double *)malloc(Nstate * sizeof(double));
        }
    }

    // printf("666666\n");

}



void sample_ele() {
    double theta[Nstate], action[Nstate];
    double thetaref[Nstate], actionref[Nstate];
    double *R0 = NULL, *Es = NULL, *C = NULL;
    double complex CC[Nstate * Nstate];
    double E[Nstate], E_diag[Nstate * Nstate];
    int i, j, k, l, nid, iref;
    double p1, x1, x2, ps1, ps2;
    double alpha_mash, beta_mash;
    double complex c_main[Nstate], sumc_main[Nstate];
    double matA[Nstate * Nstate], matB[Nstate * Nstate];
    double xe_save[Nstate], pe_save[Nstate];
    double complex gamma_cv_save[Nstate * Nstate], den_e_save[Nstate * Nstate];
    double tempdm1[Nstate * Nstate], tempdm2[Nstate * Nstate], tempdv1[Nstate], tempdv2[Nstate];
    double complex tempcm1[Nstate * Nstate], tempcm2[Nstate * Nstate];

    // Generate random numbers for theta
    for (i = 0; i < Nstate; i++) {
        theta[i] = ((double) rand() / RAND_MAX) * 2 * M_PI;
    }

    if (if_occ_rand) {
        x2 = (double) rand() / RAND_MAX;
        init_occ = 1 + floor(x2 * Nstate);
    }

    if (strcmp(method, "MFT") == 0 || strcmp(method, "mft") == 0 ||
        strcmp(method, "BCMF") == 0 || strcmp(method, "bcmf") == 0) {
        
        if (if_occ_rand) {
            x2 = (double) rand() / RAND_MAX;
            init_occ = 1 + floor(x2 * Nstate);
        }

        for (i = 0; i < Nstate; i++) {
            action[i] = 0;
        }
        action[init_occ-1] = 1;
    
        for (i = 0; i < Nstate; i++) {
            xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            
            pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
            
        }
        
        memset(gamma_cv,0,Nstate*Nstate*sizeof(double complex));
        correfun_0 = 0.5 * (xe[init_occ-1] * xe[init_occ-1] + pe[init_occ-1] * pe[init_occ-1]);

        if (type_evo >= 1) {
            for (i = 0; i < Nstate * Nstate; i++) {
                den_e[i] = 0.0 + 0.0 * I;
            }
            den_e[(init_occ-1) * Nstate + (init_occ-1)] = 1.0 + 0.0 * I;
        }

        if (if_allcf != 0) {
            for (i = 0; i < Nstate * Nstate; i++) {
                cf0[i] = den_e[i] - gamma_cv[i];
            }
        }

        // if (if_ref == 1) {
        //     for (iref = 0; iref < Nref; iref++) {
        //         for (i = 0; i < Nstate * Nstate; i++) {
        //             gamma_cv_ref[iref * Nstate * Nstate + i] = 0.0 + 0.0 * I;
        //         }

        //         if (type_evo >= 1) {
        //             for (i = 0; i < Nstate * Nstate; i++) {
        //                 den_e_ref[iref * Nstate * Nstate + i] = 0.0 + 0.0 * I;
        //             }
        //             den_e_ref[iref * Nstate * Nstate + init_occ * Nstate + init_occ] = 1.0 + 0.0 * I;
        //         }
        //     }
        // }

    } else if (strcmp(method, "eCMM") == 0 || strcmp(method, "ecmm") == 0 ||
        strcmp(method, "CMM") == 0 || strcmp(method, "cmm") == 0) {
        
        // 调用 random_prob 函数
        random_prob(Nstate, action);

        // 计算 xe 和 pe
        for (int i = 0; i < Nstate; i++) {
            xe[i] = sqrt(2 * (1 + Nstate * gamma_zpe) * action[i]) * cos(theta[i]);
            pe[i] = sqrt(2 * (1 + Nstate * gamma_zpe) * action[i]) * sin(theta[i]);
        }

        // 初始化 gamma_cv
        for (int i = 0; i < Nstate; i++) {
            for (int j = 0; j < Nstate; j++) {
                gamma_cv[i * Nstate + j] = (i == j) ? gamma_zpe : 0.0;
            }
        }

        // 计算 correfun_0
        if (index_t0 == 0) {
            correfun_0 = Nstate * (0.5 * (xe[init_occ-1] * xe[init_occ-1] + pe[init_occ-1] * pe[init_occ-1]) - gamma_zpe);
        } else if (index_t0 == 1) {
            correfun_0 = Nstate * 0.5 * (xe[index_t0_1-1] - I * pe[index_t0_1-1]) * (xe[index_t0_2-1] + I * pe[index_t0_2-1]);
        }

        // 计算 den_e
        if (type_evo >= 1) {
            for (int i = 0; i < Nstate; i++) {
                for (int j = 0; j < Nstate; j++) {
                    den_e[i * Nstate + j] = 0.5 * (xe[i] + I * pe[i]) * (xe[j] - I * pe[j]);
                }
            }
        }

        // 计算 cf0
        if (if_allcf != 0) {
            for (int i = 0; i < Nstate; i++) {
                for (int j = 0; j < Nstate; j++) {
                    cf0[i * Nstate + j] = Nstate * (den_e[i * Nstate + j] - gamma_cv[i * Nstate + j]);
                }
            }
        }
    } 

     if (ifcv == -1 || ifcv == 1) {
        for (int i = 0; i < Nstate; i++) {
            if (i == init_occ - 1) {
                gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) - 1;
                if (ifscalegamma == 1) {
                    gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) * (1 + Nstate * gamma_rescale) / (1 + Nstate * gamma_zpe) - 1;
                }
            } else {
                gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]);
                if (ifscalegamma == 1) {
                    gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) * (1 + Nstate * gamma_rescale) / (1 + Nstate * gamma_zpe);
                }
            }
        }
        
    }

    if (ifcv == -2 || ifcv == 2) {
        for (int i = 0; i < Nstate; i++) {
            if (i == init_occ - 1) {
                gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) - 1;
                if (ifscalegamma == 1) {
                    gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) * (1 + Nstate * gamma_rescale) / (1 + Nstate * gamma_zpe) - 1;
                }
            } else {
                gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]);
                if (ifscalegamma == 1) {
                    gamma_cv[i * Nstate + i] = 0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) * (1 + Nstate * gamma_rescale) / (1 + Nstate * gamma_zpe);
                }
            }
        }
        for (int i = 0; i < Nstate; i++) {
            for (int j = 0; j < Nstate; j++) {
                if (i == j || i == init_occ || j == init_occ) continue;
                gamma_cv[i * Nstate + j] = 0.5 * (xe[i] + I * pe[i]) * (xe[j] - I * pe[j]);
                if (ifscalegamma == 1) {
                    gamma_cv[i * Nstate + j] = 0.5 * (xe[i] + I * pe[i]) * (xe[j] - I * pe[j]) * (1 + Nstate * gamma_rescale) / (1 + Nstate * gamma_zpe);
                }
            }
        }
    }



    if_ad_nac = 0;
    if (sampletype == 2){
        
        dV_msmodel(R_nuc, dV);
        V_msmodel(R_nuc, V, 0.0);
        cal_NACV();

        // xe = matmul(transpose(U_d2a), xe)
        transpose(U_d2a, tempdm1, Nstate);
        memcpy(tempdv1, xe, Nstate * sizeof(double));
        dd_matmul(tempdm1, tempdv1, xe, Nstate, Nstate, 1);
        memcpy(tempdv1, pe, Nstate * sizeof(double));
        dd_matmul(tempdm1, tempdv1, pe, Nstate, Nstate, 1);
        // pe = matmul(transpose(U_d2a), pe)
        dc_matmul(tempdm1,gamma_cv,tempcm1, Nstate, Nstate, Nstate);
        cd_matmul(tempcm1, U_d2a, gamma_cv, Nstate, Nstate, Nstate);
        // gamma_cv = matmul(matmul(transpose(U_d2a), gamma_cv), U_d2a)
        
        // 这里需要实现矩阵乘法和转置操作

        if (type_evo == 1 || type_evo == 3) {
            // den_e = matmul(matmul(transpose(U_d2a), den_e), U_d2a)
            dc_matmul(tempdm1,den_e,tempcm1, Nstate, Nstate, Nstate);
            cd_matmul(tempcm1, U_d2a, den_e, Nstate, Nstate, Nstate);
        }
        // if (den_e4nuc != NULL) {
        //     // den_e4nuc = matmul(matmul(transpose(U_d2a), den_e4nuc), U_d2a)
        // }
        // if (type_evo == 4) {
        //     // G_xpconfg = matmul(transpose(U_d2a), G_xpconfg)
        // }

        if (strcmp(method, "fssh") == 0 || strcmp(method, "FSSH") == 0 || strcmp(method, "SC-FSSH") == 0 || strcmp(method, "sc-fssh") == 0 ||
            strcmp(method, "CC-FSSH") == 0 || strcmp(method, "cc-fssh") == 0 || strcmp(method, "fsshswitch") == 0 || strcmp(method, "pcsh") == 0 ||
            strcmp(method, "PCSH") == 0 || strcmp(method, "PCSH-NAF") == 0 || strcmp(method, "pcsh-naf") == 0 || strcmp(method, "BCSH") == 0 ||
            strcmp(method, "bcsh") == 0 || strcmp(method, "BCSH-NAF") == 0 || strcmp(method, "bcsh-naf") == 0) {

            double x2 = (double)rand() / RAND_MAX;
            double ps1, ps2;
            for (int i = 0; i < Nstate; i++) {
                id_state = i;
                
                if (i == 0) {
                    ps1 = 0;
                    ps2 = U_d2a[init_occ * Nstate + 0] * U_d2a[init_occ * Nstate + 0];
                } else {
                    ps1 = 0;
                    for (int k = 0; k < i - 1; k++) {
                        ps1 += U_d2a[init_occ * Nstate + k] * U_d2a[init_occ * Nstate + k];
                    }
                    ps2 = ps1 + U_d2a[init_occ * Nstate + i] * U_d2a[init_occ * Nstate + i];
                }
                if (x2 >= ps1 && x2 < ps2) {
                    break;
                }
            }
        } else if (strcmp(method, "mash") == 0 || strcmp(method, "MASH") == 0) {
            
            if (0.5 * (xe[0] * xe[0] + pe[0] * pe[0]) >= 0.5) {
                id_state = 0;
            } else {
                id_state = 1;
            }
            measure_mash = fabs(xe[0] * xe[0] + pe[0] * pe[0] - xe[1] * xe[1] - pe[1] * pe[1]);
           
            memcpy(U0_mash, U_d2a, sizeof(U_d2a));
            memset(rho0_mash, 0, sizeof(rho0_mash));
            if (0.5 * (xe[0] * xe[0] + pe[0] * pe[0]) >= 0.5) {
                rho0_mash[0] = 1;
            } else {
                rho0_mash[3] = 1;
            }
            rho0_mash[1] = 0.5 * (xe[0] + I * pe[0]) * (xe[1] - I * pe[1]);
            rho0_mash[2] = 0.5 * (xe[0] - I * pe[0]) * (xe[1] + I * pe[1]);
        } else if (strcmp(method, "unsmash") == 0 || strcmp(method, "UNSMASH") == 0 || strcmp(method, "unSMASH") == 0) {
            memcpy(U0_unsmash, U_d2a, sizeof(U0_unsmash));
        } else if (strcmp(method, "mash-mf") == 0 || strcmp(method, "MASH-MF") == 0) {
            
            if (0.5 * (xe[0] * xe[0] + pe[0] * pe[0]) >= 0.5) {
                id_state = 0;
            } else {
                id_state = 1;
            }
        } else if (strcmp(method, "ms-mash") == 0 || strcmp(method, "MS-MASH") == 0 || strcmp(method, "ms-mash-focus") == 0 || strcmp(method, "MS-MASH-focus") == 0 ||
                   strcmp(method, "mash-mf3") == 0 || strcmp(method, "MASH-MF3") == 0) {
           
            // double max_val = xe[0] * xe[0] + pe[0] * pe[0];
            // for (int i = 1; i < Nstate; i++) {
            //     double val = xe[i] * xe[i] + pe[i] * pe[i];
            //     if (val > max_val) {
            //         max_val = val;
            //         id_state = i;
            //     }
            // }
            for (int i = 0; i < Nstate; i++){
                tempdv1[i] = xe[i] * xe[i] + pe[i] * pe[i];
            }
            id_state = maxloc(tempdv1, Nstate);

        // } else if (strcmp(method, "mf2-sh") == 0 || strcmp(method, "MF2-SH") == 0) {
        //     int id_state = 0;
        //     double max_val = xe[0] * xe[0] + pe[0] * pe[0];
        //     for (int i = 1; i < Nstate; i++) {
        //         double val = xe[i] * xe[i] + pe[i] * pe[i];
        //         if (val > max_val) {
        //             max_val = val;
        //             id_state = i;
        //         }
        //     }
        // } else if (strcmp(method, "cvsh") == 0 || strcmp(method, "CVSH") == 0) {
        //     double action[Nstate];
        //     for (int i = 0; i < Nstate; i++) {
        //         if (type_evo == 1) {
        //             action[i] = den_e[i * Nstate + i];
        //         } else {
        //             action[i] = (xe[i] * xe[i] + pe[i] * pe[i]) / 2 - gamma_cv[i * Nstate + i];
        //         }
        //     }
        //     int id_state = 0;
        //     double max_val = action[0];
        //     for (int i = 1; i < Nstate; i++) {
        //         if (action[i] > max_val) {
        //             max_val = action[i];
        //             id_state = i;
        //         }
        //     }
        }
    }
}



void cal_correfun() {
    int i, j, k, l;
    double x1, x2;
    double ct_sqc;
    double x0[Nstate], p0[Nstate], pop_test[Nstate], rad_test;
    double theta_cycle[2], theta_sw, varphi_sw, psi_gp;
    double basisv1[3], basisv2[3], pvec_sw[3];
    double weight_switch[2], A[2 * 2], b[2], E_vec[2], E_mat[2 * 2], U[2 * 2];
    double complex c_st[Nstate], sum_cst[Nstate];
    double action[Nstate], theta[Nstate];
    int i_st;
    double rdm_st;
    double alpha_mash, beta_mash;
    double tempv[Nstate];
    double complex tempcm[Nstate*Nstate];
    double tempdm[Nstate*Nstate],tempdm2[Nstate*Nstate];
    double complex tempcm2[Nstate*Nstate];
    int slavecore_id = athread_get_id(-1);
    
    
    // tempmat1[Nstate*Nstate],tempmat2[Nstate*Nstate],tempmat3[Nstate*Nstate]



    if (sampletype == 2 && typeevo_ele == 0) {
        // xe = matmul(U_d2a, xe)
        // matmul(U_d2a, xe, Nstate, Nstate, 1);
        
        memcpy(tempv,xe,Nstate*sizeof(double));
        dd_matmul(U_d2a,tempv,xe,Nstate,Nstate,1);
        // pe = matmul(U_d2a, pe)
        // matmul(U_d2a, pe, Nstate, Nstate, 1);
        memcpy(tempv,pe,Nstate*sizeof(double));
        dd_matmul(U_d2a,tempv,pe,Nstate,Nstate,1);
        // gamma_cv = matmul(matmul(U_d2a, gamma_cv), transpose(U_d2a))
        // matmul(U_d2a, gamma_cv, Nstate, Nstate, Nstate);
        
        memcpy(tempcm,gamma_cv,Nstate*Nstate*sizeof(double complex));
        
        transpose(U_d2a,tempdm,Nstate);
        cd_matmul(gamma_cv,tempdm,tempcm2,Nstate,Nstate,Nstate);
        dc_matmul(U_d2a,tempcm2,gamma_cv,Nstate,Nstate,Nstate);

        
        if (type_evo == 1 || type_evo == 3) {
            memcpy(tempcm,den_e,Nstate*Nstate*sizeof(double complex));
            cd_matmul(den_e,tempdm,tempcm2,Nstate,Nstate,Nstate);
            dc_matmul(U_d2a,tempcm2,den_e,Nstate,Nstate,Nstate);
        }
        // if (type_evo == 4) {
        //     // G_xpconfg = matmul(U_d2a, G_xpconfg)
        //     matmul(U_d2a, G_xpconfg, Nstate, Nstate, 1);
        // }
    }

    
    if (strcmp(method, "MFT") == 0 || strcmp(method, "mft") == 0 ||
        strcmp(method, "BCMF") == 0 || strcmp(method, "bcmf") == 0) {
        
        if (type_evo == 1 || type_evo == 3) {
            for (i = 0; i < Nstate * Nstate; i++) {
                correfun_t[i] = den_e[i];
            }
        } else {
            for (i = 0; i < Nstate; i++) {
                for (j = 0; j < Nstate; j++) {
                    correfun_t[i * Nstate + j] = 0.5 * (xe[i] + I * pe[i]) * (xe[j] - I * pe[j]);
                }
            }
        }

    } else if (strcmp(method, "eCMM") == 0 || strcmp(method, "ecmm") == 0 ||
        strcmp(method, "CMM") == 0 || strcmp(method, "cmm") == 0) {
        
        if (type_evo == 1 || type_evo == 3) {
            for (int i = 0; i < Nstate * Nstate; i++) {
                correfun_t[i] = den_e[i] * (1.0 + Nstate) / pow(1 + Nstate * gamma_zpe, 2);
            }
            for (int i = 0; i < Nstate; i++) {
                correfun_t[i * Nstate + i] -= (1.0 - gamma_zpe) / (1 + Nstate * gamma_zpe);
            }
        } else {
            for (int i = 0; i < Nstate; i++) {
                for (int j = 0; j < Nstate; j++) {
                    correfun_t[i * Nstate + j] = (1.0 + Nstate) / (2 * pow(1 + Nstate * gamma_zpe, 2)) * (xe[i] + I * pe[i]) * (xe[j] - I * pe[j]) - (1.0 - gamma_zpe) / (1 + Nstate * gamma_zpe) * (i == j ? 1 : 0);
                }
            }
        }

    }





    if (sampletype == 2 && typeevo_ele == 0) {
        // xe = matmul(U_d2a, xe)
        // matmul(U_d2a, xe, Nstate, Nstate, 1);
        transpose(U_d2a,tempdm,Nstate);
        memcpy(tempv,xe,Nstate*sizeof(double));
        dd_matmul(tempdm,tempv,xe,Nstate,Nstate,1);
        // pe = matmul(U_d2a, pe)
        // matmul(U_d2a, pe, Nstate, Nstate, 1);
        memcpy(tempv,pe,Nstate*sizeof(double));
        dd_matmul(tempdm,tempv,pe,Nstate,Nstate,1);
        // gamma_cv = matmul(matmul(U_d2a, gamma_cv), transpose(U_d2a))
        // matmul(U_d2a, gamma_cv, Nstate, Nstate, Nstate);
       
        memcpy(tempcm,gamma_cv,Nstate*Nstate*sizeof(double complex));
        cd_matmul(gamma_cv,U_d2a,tempcm2,Nstate,Nstate,Nstate);
        dc_matmul(tempdm,tempcm2,gamma_cv,Nstate,Nstate,Nstate);

        
        if (type_evo == 1 || type_evo == 3) {
            memcpy(tempcm,den_e,Nstate*Nstate*sizeof(double complex));
            cd_matmul(den_e,U_d2a,tempcm2,Nstate,Nstate,Nstate);
            dc_matmul(tempdm,tempcm2,den_e,Nstate,Nstate,Nstate);
        }
        // if (type_evo == 4) {
        //     // G_xpconfg = matmul(U_d2a, G_xpconfg)
        //     matmul(U_d2a, G_xpconfg, Nstate, Nstate, 1);
        // }
    }


}


void cal_propagator(int Nstate, double *H, double dt, double complex *U) {
    int i, j;
    double E[Nstate], C[Nstate * Nstate];
    double sineig[Nstate * Nstate], coseig[Nstate * Nstate];
    double real_pro[Nstate * Nstate], img_pro[Nstate * Nstate];
    
    
    // 调用dia_symmat函数
    dia_symmat(Nstate, H, E, C);
    
    // 初始化sineig和coseig
    // for (i = 0; i < Nstate * Nstate; i++) {
    //     sineig[i] = 0.0;
    //     coseig[i] = 0.0;
    // }
    memset(sineig,0,Nstate*Nstate*sizeof(double));
    memset(coseig,0,Nstate*Nstate*sizeof(double));
    
    // 计算sineig和coseig
    for (i = 0; i < Nstate; i++) {
        sineig[i * Nstate + i] = -1.0*sin(E[i] * dt);
        coseig[i * Nstate + i] = cos(E[i] * dt);
    }
    
    // 计算real_pro和img_pro
    double tempdm1[Nstate*Nstate],tempdm2[Nstate*Nstate];

    transpose(C, tempdm1, Nstate);
    dd_matmul(coseig,tempdm1,tempdm2,Nstate,Nstate,Nstate);
    dd_matmul(C,tempdm2,real_pro,Nstate,Nstate,Nstate);
    dd_matmul(sineig,tempdm1,tempdm2,Nstate,Nstate,Nstate);
    dd_matmul(C,tempdm2,img_pro,Nstate,Nstate,Nstate);
    
    
    
    // 计算U
    for (i = 0; i < Nstate * Nstate; i++) {
        U[i] = real_pro[i] + I * img_pro[i];
    }
    
    // 检查E_adia是否已分配
    if (E_adia != NULL) {
        memcpy(U_d2a,C,Nstate*Nstate*sizeof(double));
        memcpy(E_adia,E,Nstate*sizeof(double));
        // for (i = 0; i < Nstate * Nstate; i++) {
        //     U_d2a[i] = C[i];
        //     if(i<Nstate) {E_adia[i] = E[i]};
        // }
        // for (i = 0; i < Nstate; i++) {
        //     E_adia[i] = E[i];
        // }
    }
}
        
    
void evo_traj_ele(double deltat) {
    double x0[Nstate], p0[Nstate];
    double x0_mb[2 * (Nstate - 1)], p0_mb[2 * (Nstate - 1)];
    double complex propagator_unsmash[2 * 2 * (Nstate - 1)];
    int i;
    double complex tempv1[Nstate],tempv2[Nstate];
    double complex tempcm1[Nstate*Nstate],tempcm2[Nstate*Nstate];
    // if (U_d2a_old != NULL) memcpy(U_d2a_old, U_d2a, sizeof(U_d2a));
    // printf("11111\n");
    switch (type_evo) {
        case 0:
            
            if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
            if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
            memcpy(x0, xe, Nstate * sizeof(double));
            memcpy(p0, pe, Nstate * sizeof(double));
            // matmul_real_imag(Nstate, propagator, x0, p0, xe, pe);
            // printf("22222\n");
            
            cd_matmul(propagator,x0,tempv1,Nstate,Nstate,1);
            cd_matmul(propagator,p0,tempv2,Nstate,Nstate,1);
            
            for(i=0;i<Nstate;i++){
                xe[i]=creal(tempv1[i])-cimag(tempv2[i]);
                pe[i]=creal(tempv2[i])+cimag(tempv1[i]);
            }
           
            break;
        
        case 1:
            if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
            if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
            // matmul_complex(Nstate, propagator, den_e, den_e);
            
            diagger(propagator,tempcm1,Nstate);
            cc_matmul(den_e,tempcm1,tempcm2,Nstate,Nstate,Nstate);
            cc_matmul(propagator,tempcm2,den_e,Nstate,Nstate,Nstate);
            // if (inverse_kernel != NULL) {
            //     matmul_complex(Nstate, propagator, inverse_kernel, inverse_kernel);
            // }
            // break;
        
            
        // case 2:
        //     if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
        //     if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
        //     memcpy(x0, xe, Nstate * sizeof(double));
        //     memcpy(p0, pe, Nstate * sizeof(double));
        //     matmul_real_imag(Nstate, propagator, x0, p0, xe, pe);
        //     matmul_complex(Nstate, propagator, den_e, den_e);
        //     break;
        // case 3:
        //     if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
        //     if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
        //     matmul_complex(Nstate, propagator, den_e, den_e);
        //     matmul_complex(Nstate, propagator, den_e4nuc, den_e4nuc);
        //     break;
        // case 4:
        //     if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
        //     if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
        //     memcpy(x0, xe, Nstate * sizeof(double));
        //     memcpy(p0, pe, Nstate * sizeof(double));
        //     matmul_real_imag(Nstate, propagator, x0, p0, xe, pe);
        //     matmul_complex(Nstate, propagator, G_xpconfg, G_xpconfg);
        //     break;
        // case 5:
        //     if (rep == 0) cal_propagator_dia_unsmash(Nstate, V, id_state, state_unsmash, deltat, propagator_unsmash);
        //     if (rep == 1 && typeevo_ele == 0) cal_propagator_adia_unsmash(Nstate, id_state, state_unsmash, deltat, propagator_unsmash);
        //     memcpy(x0_mb, xe_mb, 2 * (Nstate - 1) * sizeof(double));
        //     memcpy(p0_mb, pe_mb, 2 * (Nstate - 1) * sizeof(double));
        //     for (i = 0; i < Nstate - 1; i++) {
        //         matmul_real_imag(2, &propagator_unsmash[4 * i], &x0_mb[2 * i], &p0_mb[2 * i], &xe_mb[2 * i], &pe_mb[2 * i]);
        //     }
        //     break;
    }

    if (ifcv == 1 || ifcv == 2) {
        // matmul_complex(Nstate, propagator, gamma_cv, gamma_cv);
        diagger(propagator,tempcm1,Nstate);
        cc_matmul(gamma_cv,tempcm1,tempcm2,Nstate,Nstate,Nstate);
        cc_matmul(propagator,tempcm2,gamma_cv,Nstate,Nstate,Nstate);
    }
}



// void evo_traj_ele_prop(double deltat){

// TODO !!!!!!!!!!!!!!!!!!!!!!!!!!

// }

void evo_traj_nucP(double deltat) {
    // 假设 P_nuc 和 force 是一维数组
   
    for (int i = 0; i < Ndof1*Ndof2; i++) {
        P_nuc[i] += force[i] * deltat;
    }

    // if (ifscaleenergy == 3) energy_conserve_naf_3(deltat);
}

void evo_traj_nucR(double deltat) {
    double x1, x2;

    if (iflangevin == 1) {
        // for (int i = 0; i < Ndof1*Ndof2; i++) {
        //     R_nuc[i] += P_nuc[i] * 0.5 * deltat / mass[i];
        // }

        // for (int i = 0; i < Ndof1; i++) {
        //     for (int j = 0; j < Ndof2; j++) {
        //         box_muller(&x2, &x1, 1.0, 0.0);
        //         int index = i * Ndof2 + j;
        //         P_nuc[index] = exp(-eta_langevin * deltat) * P_nuc[index] + sqrt((1 - exp(-2 * eta_langevin * deltat)) * mass[index] / beta) * x2;
        //     }
        // }

        // for (int i = 0; i < Ndof; i++) {
        //     R_nuc[i] += P_nuc[i] * 0.5 * deltat / mass[i];
        // }
    } else {
        // if (ifswitchforce == 3 && sum_abs_dv_adia(id_state, id_state) > 1e-30) {
        //     x2 = 0.0;
        //     for (int i = 0; i < Nstate; i++) {
        //         for (int j = 0; j < Nstate; j++) {
        //             if (i == j) continue;
        //             x2 += (0.5 * (xe[i] * xe[j] + pe[i] * pe[j]) - gamma_cv[i * Nstate + j]) * (E_adia[j] - E_adia[i]) * sum_nac_P_nuc_mass(i, j);
        //         }
        //     }
        //     for (int i = 0; i < Ndof; i++) {
        //         R_nuc[i] += P_nuc[i] * deltat / mass[i] + x2 / dv_adia[id_state * Ndof + id_state] * deltat;
        //     }
        // } else {
            for (int i = 0; i < Ndof1*Ndof2; i++) {
                R_nuc[i] += P_nuc[i] * deltat / mass[i];
            }
        // }
    }
}


void evo_traj_calProp(int igrid_cal) {
    int i, j, icfall;
    double x2;

    cal_correfun();

    // if (ifmsbranch > 0) {
    //     if (if_statetraj == -1) {
    //         for (i = 0; i < Nstate; i++) {
    //             for (j = 0; j < Nstate; j++) {
    //                 correfun_t[i * Nstate + j] -= correfun_t_oldtraj[i * Nstate + j + igrid_cal * Nstate * Nstate];
    //                 correfun_t_oldtraj[i * Nstate + j + igrid_cal * Nstate * Nstate] = correfun_t[i * Nstate + j] + correfun_t_oldtraj[i * Nstate + j + igrid_cal * Nstate * Nstate];
    //             }
    //         }
    //     } else {
    //         memcpy(&correfun_t_oldtraj[igrid_cal * Nstate * Nstate], correfun_t, Nstate * Nstate * sizeof(double));
    //     }
    // }

    if (any_isnan(correfun_t, Nstate * Nstate)) {
        N_nan_sum[igrid_cal]++;
        if (if_st_nan == 1) {
            memset(correfun_t, 0, Nstate * Nstate * sizeof(double complex));
        }
    } else {
        N_nan_sum[igrid_cal] = N_nan_sum[igrid_cal];
    }

    // printf("ii=%d\n",igrid_cal);

    if (den != NULL) {
        if (strcmp(trim(adjustl(method)), "gauss") == 0 || strcmp(trim(adjustl(method)), "genLSC") == 0 || strcmp(trim(adjustl(method)), "genlsc") == 0) {
            if (ifid == 0) {
                for (i = 0; i < Nstate; i++) {
                    for (j = 0; j < Nstate; j++) {
                        den[i * Nstate * Ngrid + j * Ngrid  + igrid_cal] += correfun_0 * correfun_t[i * Nstate + j];
                    }
                }
            } else if (ifid == 1) {
                for (i = 0; i < Nstate; i++) {
                    for (j = 0; j < Nstate; j++) {
                        den[i * Nstate * Ngrid + j * Ngrid  + igrid_cal] += correfun_0 * correfun_t[i * Nstate + j];
                    }
                }
                for (j = 0; j < Nstate; j++) {
                    den[j * Nstate * Ngrid + j * Ngrid + igrid_cal] += 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                }
            }
        } else {
            for (i = 0; i < Nstate; i++) {
                for (j = 0; j < Nstate; j++) {
                    den[i * Nstate * Ngrid + j * Ngrid + igrid_cal] += correfun_0 * correfun_t[i * Nstate + j];
                }
            }
        }
    }

    if (population!=NULL) {
        for (i = 0; i < Nstate; i++) {
            if (strcmp(trim(adjustl(method)), "gauss") == 0 || strcmp(trim(adjustl(method)), "genLSC") == 0 || strcmp(trim(adjustl(method)), "genlsc") == 0) {
                if (ifid == 0) {
                    population[i * Ngrid  + igrid_cal] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                    if (if_st_fb == 1) {
                        if (P_nuc[0] > 0) {
                            pop_fb[i * Ngrid *2  + igrid_cal*2] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                        } else {
                            pop_fb[i * Ngrid *2  + igrid_cal*2 + 1] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                        }
                    }
                } else if (ifid == 1) {
                    population[i * Ngrid  + igrid_cal] += creal(correfun_0 * correfun_t[i * Nstate + i]) + 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                    if (if_st_fb == 1) {
                        if (P_nuc[0] > 0) {
                            pop_fb[i * Ngrid *2  + igrid_cal*2] += creal(correfun_0 * correfun_t[i * Nstate + i]) + 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                        } else {
                            pop_fb[i * Ngrid *2  + igrid_cal*2 + 1] += creal(correfun_0 * correfun_t[i * Nstate + i]) + 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                        }
                    }
                }
            } else {
                population[i * Ngrid  + igrid_cal] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                if (if_st_fb == 1) {
                    if (P_nuc[0] > 0) {
                        pop_fb[i * Ngrid *2  + igrid_cal*2] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                    } else {
                        pop_fb[i * Ngrid *2  + igrid_cal*2 + 1] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                    }
                }
            }
        }
    }

    // if (allocated(cfall)) {
    //     for (i = 0; i < Nstate; i++) {
    //         for (j = 0; j < Nstate; j++) {
    //             for (int k = 0; k < Ndof1; k++) {
    //                 for (int l = 0; l < Ndof2; l++) {
    //                     cfall[i * Nstate + j + k * Nstate * Nstate + l * Nstate * Nstate * Ndof1 + igrid_cal * Nstate * Nstate * Ndof1 * Ndof2] += cf0[i * Nstate + j] * correfun_t[k * Nstate + l];
    //                 }
    //             }
    //         }
    //     }
    // }

    // if (cfeff!=NULL) {
    //     if (if_allcf == 2) {
    //         cfeff[igrid_cal] += sum(weight0, Nstate) * sum(weightt, Nstate) * sum(correfun_t, Nstate * Nstate);
    //     } else if (if_allcf == 3) {
    //         for (icfall = 0; icfall < allcf_times; icfall++) {
    //             cfweight_msmodel(weight0, weightt, beta, icfall);
    //             cfeff[igrid_cal] += sum(weight0, Nstate) * sum(weightt, Nstate) * sum(correfun_t, Nstate * Nstate);
    //         }
    //     } else if (if_allcf == 4) {
    //         for (icfall = 0; icfall < allcf_times; icfall++) {
    //             cfweight_msmodel(weight0, weightt, beta, icfall);
    //             cfeff[igrid_cal] += sum(weight0, Nstate) * sum(weightt, Nstate) * sum(correfun_t, Nstate * Nstate);
    //         }
    //     }
    // }

    if (R_nuc_mean!=NULL) {
        if (strcmp(trim(adjustl(method)), "sqc") == 0 
        || strcmp(trim(adjustl(method)), "SQC") == 0 
        || strcmp(trim(adjustl(method)), "mf3") == 0 
        || strcmp(trim(adjustl(method)), "MF3") == 0 
        || strcmp(trim(adjustl(method)), "sqc2") == 0 
        || strcmp(trim(adjustl(method)), "SQC2") == 0 
        || strcmp(trim(adjustl(method)), "sqc3") == 0 
        || strcmp(trim(adjustl(method)), "SQC3") == 0) {
            x2 = 0;
            for (i = 0; i < Nstate; i++) {
                x2 += correfun_t[i * Nstate + i];
            }
            for (i = 0; i < Ndof1*Ndof2; i++) {
                R_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * creal(correfun_0) * x2;
                P_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * creal(correfun_0) * x2;
                R2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * R_nuc[i] * creal(correfun_0) * x2;
                P2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * P_nuc[i] * creal(correfun_0) * x2;
            }
        } else if (strcmp(trim(adjustl(method)), "mash") == 0 
        || strcmp(trim(adjustl(method)), "MASH") == 0) {
            for (i = 0; i < Ndof1*Ndof2; i++) {
                R_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * creal(correfun_0) * 2 * measure_mash * creal(rho0_mash[(init_occ-1) * Nstate + (init_occ-1)]);
                P_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * creal(correfun_0) * 2 * measure_mash * creal(rho0_mash[(init_occ-1) * Nstate + (init_occ-1)]);
                R2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * R_nuc[i] * creal(correfun_0) * 2 * measure_mash * creal(rho0_mash[(init_occ-1) * Nstate + (init_occ-1)]);
                P2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * P_nuc[i] * creal(correfun_0) * 2 * measure_mash * creal(rho0_mash[(init_occ-1) * Nstate + (init_occ-1)]);
            }
        // } else if (strcmp(trim(adjustl(method)), "unsmash") == 0 || strcmp(trim(adjustl(method)), "UNSMASH") == 0 || strcmp(trim(adjustl(method)), "unSMASH") == 0) {
        //     for (i = 0; i < Ndof; i++) {
        //         R_nuc_mean[i + igrid_cal * Ndof] += R_nuc[i] * real(correfun_0) * Nstate * measure_mash * rho0_unsmash[init_occ * Nstate + init_occ];
        //         P_nuc_mean[i + igrid_cal * Ndof] += P_nuc[i] * real(correfun_0) * Nstate * measure_mash * rho0_unsmash[init_occ * Nstate + init_occ];
        //         R2_nuc_mean[i + igrid_cal * Ndof] += R_nuc[i] * R_nuc[i] * real(correfun_0) * Nstate * measure_mash * rho0_unsmash[init_occ * Nstate + init_occ];
        //         P2_nuc_mean[i + igrid_cal * Ndof] += P_nuc[i] * P_nuc[i] * real(correfun_0) * Nstate * measure_mash * rho0_unsmash[init_occ * Nstate + init_occ];
        //     }
        } else if (strcmp(trim(adjustl(method)), "mash-mf") == 0 || strcmp(trim(adjustl(method)), "MASH-MF") == 0) {
            for (i = 0; i < Ndof1*Ndof2; i++) {
                R_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * creal(correfun_0) * 2 * measure_mash;
                P_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * creal(correfun_0) * 2 * measure_mash;
                R2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * R_nuc[i] * creal(correfun_0) * 2 * measure_mash;
                P2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * P_nuc[i] * creal(correfun_0) * 2 * measure_mash;
            }
        } else {
            for (i = 0; i < Ndof1*Ndof2; i++) {
                R_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * creal(correfun_0);
                P_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * creal(correfun_0);
                R2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += R_nuc[i] * R_nuc[i] * creal(correfun_0);
                P2_nuc_mean[i + igrid_cal * Ndof1*Ndof2] += P_nuc[i] * P_nuc[i] * creal(correfun_0);
            }
        }
    }

    // if (ifmsbranch > 0) {
    //     if (if_statetraj == -1) {
    //         for (i = 0; i < Ndof; i++) {
    //             R_nuc_mean[i + igrid_cal * Ndof] -= R_nuc_oldtraj[i + igrid_cal * Ndof];
    //             P_nuc_mean[i + igrid_cal * Ndof] -= P_nuc_oldtraj[i + igrid_cal * Ndof];
    //         }
    //     }
    //     memcpy(&R_nuc_oldtraj[igrid_cal * Ndof], R_nuc, Ndof * sizeof(double));
    //     memcpy(&P_nuc_oldtraj[igrid_cal * Ndof], P_nuc, Ndof * sizeof(double));
    // }

    
    if (energy_est != NULL) {
        double x2 = 0.0;
        double Ekin=0.0;
        // V_msmodel(R_nuc, V, t_now); // Uncomment if needed
        for(int i=0;i<Ndof1*Ndof2;i++){
            Ekin+=0.5*P_nuc[i]*P_nuc[i]/mass[i];
        }
        if (rep == 1 && sampletype == 0) {
            // cal_NACV(); // Uncomment if needed
            for (int i = 0; i < Nstate; i++) {
                x2 += (E_adia[i] + Ekin) * correfun_0 * correfun_t[i * Nstate + i];
            }
        } else {
            for (int i = 0; i < Nstate; i++) {
                for (int j = 0; j < Nstate; j++) {
                    if (i == j) {
                        x2 += Ekin * correfun_0 * correfun_t[i * Nstate + i];
                    }
                    x2 += V[i * Nstate + j] * correfun_0 * correfun_t[j * Nstate + i];
                }
            }
        }

        if (ifswitchforce >= 1) {
            // cal_NACV(); // Uncomment if needed
            energy_est[igrid_cal] += Ekin + E_adia[id_state];
        } else {
            if (strcmp(trim(adjustl(method)), "FSSH") == 0 || strcmp(trim(adjustl(method)), "fssh") == 0 ||
                strcmp(trim(adjustl(method)), "mash") == 0 || strcmp(trim(adjustl(method)), "MASH") == 0 ||
                strcmp(trim(adjustl(method)), "ms-mash") == 0 || strcmp(trim(adjustl(method)), "MS-MASH") == 0 ||
                strcmp(trim(adjustl(method)), "unsmash") == 0 || strcmp(trim(adjustl(method)), "UNSMASH") == 0 ||
                strcmp(trim(adjustl(method)), "unSMASH") == 0) {
                energy_est[igrid_cal] += Ekin + E_adia[id_state];
            } else {
                energy_est[igrid_cal] += x2;
            }
        }
    }

    // if (count_st != NULL) {
    //     for (int i = 0; i < Nstate; i++) {
    //         count_st[i * igrid_cal] += count_pertraj;
    //     }
    //     count_pertraj = 0;
    

}


void energy_conserve_naf_1(double E0, double *dE_naf) {
    *dE_naf = E0 - E_adia[id_state];
    if (*dE_naf >= 0) {
        double P_nuc_sum = 0.0;
        for (int i = 0; i < Ndof1*Ndof2; i++) {
            P_nuc_sum += P_nuc[i] * P_nuc[i] / mass[i];
        }
        double factor = sqrt(*dE_naf / (0.5 * P_nuc_sum));
        for (int i = 0; i < Ndof1*Ndof2; i++) {
            P_nuc[i] *= factor;
        }
    }
}


void energy_conserve_naf_3(double deltat) {
    double dE_naf, x2 = 0.0;
    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) {
            if (i == j) continue;
            double sum_nac = 0.0;
            for (int k = 0; k < Ndof1*Ndof2; k++) {
                sum_nac += nac[i*Nstate*Ndof1*Ndof2+j*Ndof1*Ndof2+k] * P_nuc[k] / mass[k];
            }
            x2 += (0.5 * (xe[i] * xe[j] + pe[i] * pe[j]) - creal(gamma_cv[i*Nstate+j])) * (E_adia[j] - E_adia[i]) * sum_nac;
        }
    }
    double P_nuc_sum = 0.0;
    for (int i = 0; i < Ndof1*Ndof2; i++) {
        P_nuc_sum += P_nuc[i] * P_nuc[i] / mass[i];
    }
    double factor = x2 / P_nuc_sum * deltat;
    for (int i = 0; i < Ndof1*Ndof2; i++) {
        P_nuc[i] += P_nuc[i] * factor;
    }
}

void evo_traj_algorithm1(double deltat) {
    int slavecore_id=athread_get_id(-1);
    double tempv[Nstate];
    double tempdm[Nstate*Nstate];
    
    evo_traj_nucP(deltat / 2);
   
    if (scaleenergy_type == 3) energy_conserve_naf_3(deltat / 2);
    
    // if (strcmp(trim(adjustl(msmodelname)), "LZ") == 0) {
    //     memcpy(P_nuc, Pinit_LZ, sizeof(double));
    // }
    evo_traj_nucR(deltat);

    // if(slavecore_id==0){
    //     if(rep==1){
    //         memcpy(tempv,xe,Nstate*sizeof(double));
    //         dd_matmul(U_d2a,tempv,xe,Nstate,Nstate,1);
    //         memcpy(tempv,pe,Nstate*sizeof(double));
    //         dd_matmul(U_d2a,tempv,pe,Nstate,Nstate,1);
    //     }
    //     printf("s11 %18.8E  %18.8E  %18.8E  %18.8E\n", R_nuc[0], P_nuc[0], xe[0], pe[0]);
    //     if(rep==1){
    //         transpose(U_d2a,tempdm,Nstate);
    //         memcpy(tempv,xe,Nstate*sizeof(double));
    //         dd_matmul(tempdm,tempv,xe,Nstate,Nstate,1);
    //         memcpy(tempv,pe,Nstate*sizeof(double));
    //         dd_matmul(tempdm,tempv,pe,Nstate,Nstate,1);
    //     }
    // }

    
   
    dV_msmodel(R_nuc, dV);
    V_msmodel(R_nuc, V, t_now);
    if (rep == 1) cal_NACV();
    evo_traj_ele(deltat);
    // exit(-1);
    // if(slavecore_id==0){
    //     if(rep==1){
    //         memcpy(tempv,xe,Nstate*sizeof(double));
    //         dd_matmul(U_d2a,tempv,xe,Nstate,Nstate,1);
    //         memcpy(tempv,pe,Nstate*sizeof(double));
    //         dd_matmul(U_d2a,tempv,pe,Nstate,Nstate,1);
    //     }
    //     printf("s22 %18.8E  %18.8E  %18.8E  %18.8E\n", R_nuc[0], P_nuc[0], xe[0], pe[0]);
    //     if(rep==1){
    //         transpose(U_d2a,tempdm,Nstate);
    //         memcpy(tempv,xe,Nstate*sizeof(double));
    //         dd_matmul(tempdm,tempv,xe,Nstate,Nstate,1);
    //         memcpy(tempv,pe,Nstate*sizeof(double));
    //         dd_matmul(tempdm,tempv,pe,Nstate,Nstate,1);
    //     }
    // }
    
    cal_force();
    
    evo_traj_nucP(deltat / 2);

    if (scaleenergy_type == 3) energy_conserve_naf_3(deltat / 2);
}



void evo_traj_savetraj() {
    // 假设所有变量已经在def.h中声明
    memcpy(R_nuc_old, R_nuc, sizeof(R_nuc));
    memcpy(P_nuc_old, P_nuc, sizeof(P_nuc));
    memcpy(xe_old, xe, sizeof(xe));
    memcpy(pe_old, pe, sizeof(pe));
    if (den_e != NULL) memcpy(den_e_old, den_e, sizeof(den_e));
    memcpy(gamma_cv_old, gamma_cv, sizeof(gamma_cv));
    memcpy(V_old, V, sizeof(V));
    memcpy(dV_old, dV, sizeof(dV));
    if (E_adia != NULL) memcpy(E_adia_old, E_adia, sizeof(E_adia));
    if (dv_adia != NULL) memcpy(dv_adia_old, dv_adia, sizeof(dv_adia));
    if (nac != NULL) memcpy(nac_old, nac, sizeof(nac));
    if (U_d2a != NULL) memcpy(U_d2a_old, U_d2a, sizeof(U_d2a));
    if (U_ref != NULL) memcpy(U_ref_old, U_ref, sizeof(U_ref));
    if (nac_check != NULL) memcpy(nac_check_old, nac_check, sizeof(nac_check));
    id_state_old = id_state;

    // if (strcmp(trim(adjustl(msmodelname)), "mole") == 0) {
    //     if (file_exists("imomap.dat")) {
    //         system("cp imomap.dat imomap_save.dat");
    //     }
    // }
}

void evo_traj_back() {
    // 假设所有变量已经在def.h中声明
    memcpy(R_nuc, R_nuc_old, sizeof(R_nuc));
    memcpy(P_nuc, P_nuc_old, sizeof(P_nuc));
    memcpy(xe, xe_old, sizeof(xe));
    memcpy(pe, pe_old, sizeof(pe));
    if (den_e_old != NULL) memcpy(den_e, den_e_old, sizeof(den_e));
    memcpy(gamma_cv, gamma_cv_old, sizeof(gamma_cv));
    memcpy(V, V_old, sizeof(V));
    memcpy(dV, dV_old, sizeof(dV));
    if (E_adia_old != NULL) {
        memcpy(E_adia, E_adia_old, sizeof(E_adia));
        memcpy(dv_adia, dv_adia_old, sizeof(dv_adia));
        memcpy(nac, nac_old, sizeof(nac));
        memcpy(U_d2a, U_d2a_old, sizeof(U_d2a));
        memcpy(U_ref, U_ref_old, sizeof(U_ref));
        memcpy(nac_check, nac_check_old, sizeof(nac_check));
    }
    id_state = id_state_old;

    // if (strcmp(trim(adjustl(msmodelname)), "mole") == 0) {
    //     if (file_exists("imomap_save.dat")) {
    //         system("cp imomap_save.dat imomap.dat");
    //     }
    // }
}


void evo_traj_new(int itraj) {
    int i_re, igrid;
    int nstep, itime, icfall, iref;
    int i, j, k;//, i_lbd, hop_lbd;
    // double x0[Nstate], p0[Nstate], bound, deltaE, deltaE_all[Nstate];
    double tt1, tt2, x1, x2, sumpop, diagden[Nstate];
    // double E_diag[Nstate * Nstate];
    // double P_s[Ndof1 * Ndof2], P_s_all[Nstate * Ndof1 * Ndof2], P_s_main[Ndof1 * Ndof2], R_s[Ndof1 * Ndof2], f_s[Ndof1 * Ndof2], pex, pt, proj[Nstate];
    // double Etot, Ekin, Epot, dE;
    // double vmp[Nstate * Nstate];
    // double Rinit, Pinit;
    // double complex rho_save[Nstate * Nstate], Aforsave[Nstate * Nstate];
    // double P_nuc_BA_old[Ndof1 * Ndof2];
    // double P_nuc_old[Ndof1 * Ndof2 * memorylength], R_nuc_old[Ndof1 * Ndof2 * memorylength], E_adia_old[Nstate * memorylength], nac_old[Nstate * Nstate * Ndof1 * Ndof2 * memorylength], force_old[Ndof1 * Ndof2 * memorylength], deltaE_all_old[Nstate * memorylength];
    // int loc_bak, id_memory, jstate;
    bool if_bak;
    // long long *approachtimes;
    // double action_switchcv[Nstate], off_gamma_cv, Vij;
    // double complex c_main[Nstate];
    // int i_max, npack, index_state[Nstate];
    // double *depack;
    // int *index_pack;
    double dt_evo;
    // int nstep_small, istep_small;
    // double t_now_small;
    // bool alive;
    int slavecore_id;
    slavecore_id=athread_get_id(-1);
    
    if_bak = false;
    itime_save = 0;

    count_pertraj = 0;

    t_now = 0;
    nstep = (int)(ttot / dt);

    V_msmodel(R_nuc, V, 0.0);
    dV_msmodel(R_nuc, dV);
    if (rep == 1) cal_NACV();
    //debug
    // int slavecore_id;
    // slavecore_id=athread_get_id(-1);
    // if(slavecore_id == 0){
        // for (i=0; i<Nstate; i++){
            // printf("E_adia(%d)=%18.8E\n",i,E_adia[i]);
            // for (j=0;j<Nstate;j++){
                // printf("U(%d,%d)=%18.8E\n",i,j,U_d2a[i*Nstate+j]);
    //             for (k=0;k<1;k++){
    //                 for (int l=0;l<1;l++){
    //                     printf("dv_adia(%d,%d,%d,%d)=%18.8E\n",i,j,k,l,dv_adia[i*Nstate*Ndof1*Ndof2+j*Ndof1*Ndof2+k*Ndof2+l]);
    //                     printf("nac(%d,%d,%d,%d)=%18.8E\n",i,j,k,l,nac[i*Nstate*Ndof1*Ndof2+j*Ndof1*Ndof2+k*Ndof2+l]);
    //                 }
    //             }
            
    //         }
    //     }
    // }
    
   

    i_re = Nbreak;
    igrid = 0;

    if (ifscaleenergy > 0) {
        E_conserve = 0.0;
        for(i=0;i<Ndof1*Ndof2;i++){
            E_conserve += P_nuc[i] * P_nuc[i] / mass[i];
        }
        E_conserve *= 0.5;
        E_conserve += E_adia[id_state];
    }

    if (ifscaleenergy > 0) {
        switch (ifscaleenergy) {
            case 1:
            case 4:
            case 5:
                scaleenergy_type = 1;
                break;
            case 3:
                scaleenergy_type = 3;
                break;
        }
    }

    

    cal_force();

    // debug
    // int slavecore_id;
    // slavecore_id=athread_get_id(-1);
    // if(slavecore_id == 0){
    //     for (i=0; i<Ndof1; i++){
            
    //         for (j=0;j<Ndof2;j++){
    //             printf("force(%d,%d)=%18.8E\n",i,j,force[i*Ndof1+j]);
    //         }
    //     }
    // }
    //  exit(1);


    itime = 0;


    //debug
    // int slavecore_id;
    // slavecore_id=athread_get_id(-1);
    // if(slavecore_id == 0) printf("itime=%d,force=%18.8E\n",itime,force[299]);

   
    
    while (itime <= nstep) {
        if (i_re >= Nbreak) {
            evo_traj_calProp(igrid);

            // if(slavecore_id == 10) printf("%18.8E %18.8E %18.8E\n", t_now, R_nuc[0], P_nuc[0]);
            
            timegrid[igrid] = t_now;
            igrid++;
            i_re = 0;

            //debug
            //  if(slavecore_id == 10) printf("%18.8E %18.8E %18.8E %18.8E %18.8E \n", R_nuc[0],P_nuc[0],xe[0],pe[0],force[0]);
                
            
            // if(slavecore_id == 10) printf("%18.8E %18.8E\n", force[0],force[299]);
            

        }

        evo_traj_savetraj();

        dt_evo = dt;


        

        switch (type_algorithm) {
            case 1:
            
                evo_traj_algorithm1(dt_evo);
                break;
            // case 2:
            //     evo_traj_algorithm2(dt_evo);
            //     break;
            // case 3:
            //     evo_traj_algorithm3(dt_evo);
            //     break;
            // case 4:
            //     evo_traj_algorithm4(dt_evo);
            //     break;
            // case 5:
            //     evo_traj_algorithm5(dt_evo, n_step_algo5);
            //     break;
            // case 6:
            //     evo_traj_algorithm6(dt_evo);
            //     break;
            // case 7:
            //     evo_traj_algorithm7(dt_evo);
            //     break;
            // case 8:
            //     evo_traj_algorithm8(dt_evo);
            //     break;
            // case 9:
            //     evo_traj_algorithm9(dt_evo);
            //     break;
            // case 10:
            //     evo_traj_algorithm10(dt_evo);
            //     break;
        }

        //debug
        // if(slavecore_id == 0) printf("itime=%d,force=%18.8E\n",itime,force[299]);

        

        // if (ifscaleenergy > 0) {
        //     if (scaleenergy_type == 1) energy_conserve_naf_1(E_conserve, deltaE);
        //     if (ifscaleenergy == 4 && scaleenergy_type == 1 && deltaE < 0) {
        //         nstep_small = 2;
        //         dt_evo /= 2;
        //         t_now_small = 0;
        //         evo_traj_back();
        //         for (istep_small = 1; istep_small <= nstep_small; istep_small++) {
        //             switch (type_algorithm) {
        //                 case 1:
        //                     evo_traj_algorithm1(dt_evo);
        //                     break;
        //                 case 2:
        //                     evo_traj_algorithm2(dt_evo);
        //                     break;
        //             }
        //             if (scaleenergy_type == 1) energy_conserve_naf_1(E_conserve, deltaE);
        //             if (scaleenergy_type == 1 && deltaE < 0) {
        //                 if (dt_evo > dt / 1024) {
        //                     nstep_small *= 2;
        //                     dt_evo /= 2;
        //                     evo_traj_back();
        //                 }
        //             }
        //         }
        //         if (scaleenergy_type == 3) scaleenergy_type = 1;
        //     }
        // }

        t_now = (itime + 1) * dt;
        i_re++;
        itime++;
        // if (ifzpecorr > 0) zpecorr_msmodel(P_nuc, R_nuc, ifzpecorr);

        // if (strcmp(trim(adjustl(msmodelname)), "morse3") == 0 && ifhardwall == 1) {
        //     if (P_nuc[0] < 0 && R_nuc[0] < 0) P_nuc[0] = -P_nuc[0];
        //     if (if_ref == 1) {
        //         for (iref = 0; iref < Nref; iref++) {
        //             if (P_nuc_ref[iref * Ndof1 * Ndof2] < 0 && R_nuc_ref[iref * Ndof1 * Ndof2] < 0) P_nuc_ref[iref * Ndof1 * Ndof2] = -P_nuc_ref[iref * Ndof1 * Ndof2];
        //         }
        //     }
        // }
    }

    if (if_Pdis == 1) {
        if (strcmp(method, "mash") == 0 || strcmp(method, "MASH") == 0) {
    
            for(i=0;i<s_N;i++){
                // expisP[i] += correfun_0 * cexp(I * P_nuc[0] * s[i]) * 2 * rho0_mash[(init_occ-1) * Nstate + (init_occ-1)] * measure_mash;
                expisP[i] += correfun_0 *(cos(P_nuc[0] * s[i]) + I * sin(P_nuc[0] * s[i])) * 2 * rho0_mash[(init_occ-1) * Nstate + (init_occ-1)] * measure_mash;
            }
                // break;
        } else if (strcmp(method, "unsmash") == 0 ||
                  strcmp(method, "UNSMASH") == 0 ||
                  strcmp(method, "unSMASH") == 0 ){
            for(i=0;i<s_N;i++){
                // expisP[i] += correfun_0 * cexp(I * P_nuc[0] * s[i]) * Nstate * rho0_unsmash[(init_occ-1) * Nstate + (init_occ-1)] * measure_mash;
                expisP[i] += correfun_0 * (cos(P_nuc[0] * s[i]) + I * sin(P_nuc[0] * s[i])) * Nstate * rho0_unsmash[(init_occ-1) * Nstate + (init_occ-1)] * measure_mash;
            }
                // break;
        } else if (strcmp(method, "mash-mf") == 0 ||
                    strcmp(method, "MASH-MF") == 0 ){
            for(i=0;i<s_N;i++){
                // expisP[i] += correfun_0 * cexp(I * P_nuc[0] * s[i]) * 2 * measure_mash;
                expisP[i] += correfun_0 * (cos(P_nuc[0] * s[i]) + I* sin(P_nuc[0] * s[i])) * 2 * measure_mash;
            }
                // break;
        } else if (strcmp(method, "sqc") == 0 ||
                   strcmp(method, "SQC") == 0 ||                 
                   strcmp(method, "mf3") == 0 ||
                   strcmp(method, "MF3") == 0 ||
                   strcmp(method, "sqc2") == 0 ||
                   strcmp(method, "SQC2") == 0 ||
                   strcmp(method, "sqc3") == 0 ||
                   strcmp(method, "SQC3") == 0 ){
            x2 = 0;
            for (i = 0; i < Nstate; i++) {
                x2 += creal(correfun_t[i * Nstate + i]);
            }
            for(i=0;i<s_N;i++){
                // expisP[i] += correfun_0 * cexp(I * P_nuc[0] * s[i]) * x2;
                expisP[i] += correfun_0 * ( cos( P_nuc[0] * s[i]) + I * sin( P_nuc[0] * s[i]) ) * x2;
            }
                // break;

        } else {
                for(i=0;i<s_N;i++){
                    // expisP[i] += correfun_0 * cexp(I * P_nuc[0] * s[i]);
                    expisP[i] += correfun_0 * ( cos( P_nuc[0] * s[i]) + I * sin( P_nuc[0] * s[i]) );
                }
                // break;
        }
    }
}


void cal_force() {
    int iref, i, j;
    double frdm, x2;

    // force = 0.0;
    
    memset(force,0,Ndof1*Ndof2*sizeof(double));

    if (strcmp(trim(adjustl(method)), "FSSH") == 0 || strcmp(trim(adjustl(method)), "fssh") == 0) {
        // cal_force_fssh();
        // force = -dv_adia[id_state][id_state];
    } else {
        // if (if_ref == 1) {
        //     cal_force_mf_ref();
        // } else if (if_1st > 0) {
        //     cal_force_mf_1st();
        // } else if (ifmsbranch > 0) {
        //     cal_force_msbranch();
        // } else if (ifswitchforce == 1) {
        //     cal_force_switch();
        // } else if (ifswitchforce == 2) {
        //     cal_force_switch2();
        // } else if (ifswitchforce == 3) {
        //     cal_force_switch();
        // } else if (ifswitchforce == 4) {
        //     cal_force_switch4();
        // } else if (ifswitchforce == 5) {
        //     cal_force_switch5();
        // } else if (ifswitchforce == 6) {
        //     cal_force_switch6();
        // } else if (ifswitchforce == 7) {
        //     cal_force_switch7();
        // } else if (ifswitchforce == 8) {
        //     cal_force_switch8();
        // } else if (ifmashforce > 0) {
        //     cal_force_mashforce();
        // } else {
            cal_force_mf();
        // }
    }
    
    if (forcetype == 1) {
        nucforce_msmodel(R_nuc, force_nuc);
        for(i=0 ; i<Ndof1*Ndof2; i++){
            force[i] -= force_nuc[i];
        }
        // if (if_ref == 1) {
        //     for (iref = 1; iref <= Nref; iref++) {
        //         nucforce_msmodel(R_nuc_ref[iref], force_nuc_ref[iref]);
        //         force_ref[iref] -= force_nuc_ref[iref];
        //     }
        // }
    }
    
}



void cal_force_mf() {
    int i, j, k;
    double force_trace[Ndof1 * Ndof2];

    if (if_traceless_force == 1) {
        // for (i = 0; i < Ndof1 * Ndof2; i++) {
        //     force_trace[i] = 0;
        // }
        memset(force_trace,0,Ndof1*Ndof2*sizeof(double));
        for (i = 0; i < Nstate; i++) {
            for (j = 0; j < Ndof1 * Ndof2; j++) {
                force_trace[j] += dV[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + j] / Nstate;
            }
        }
        for (j = 0; j < Ndof1 * Ndof2; j++) {
            for (i = 0; i < Nstate; i++) {
                dV[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + j] -= force_trace[j];
            }
            force[j] = -force_trace[j];
        }
        // for (j = 0; j < Ndof1 * Ndof2; j++) {
        //     force[j] = -force_trace[j];
        // }
    }

     
   

    switch (type_evo) {
        case 0:
        case 2:
            if (rep == 0) {
                if (calforcetype == 1) {
                    for (i = 0; i < Nstate; i++) {
                        for (j = 0; j < Ndof1 * Ndof2; j++) {
                            force[j] -= dV[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + j] * ((xe[i] * xe[i] + pe[i] * pe[i]) * 0.5 - creal(gamma_cv[i * Nstate + i]));
                        }
                    }
                } else {
                    for (i = 0; i < Nstate; i++) {
                        for (j = 0; j < Nstate; j++) {
                            for (k = 0; k < Ndof1 * Ndof2; k++) {
                                force[k] -= dV[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k] * ((xe[i] * xe[j] + pe[i] * pe[j]) * 0.5 - creal(gamma_cv[i * Nstate + j]));
                            }
                        }
                    }
                }
            } else if (rep == 1) {
                for (i = 0; i < Nstate; i++) {
                    for (j = 0; j < Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < Ndof1 * Ndof2; k++) {
                                force[k] -= (0.5 * (xe[i] * xe[i] + pe[i] * pe[i]) - creal(gamma_cv[i * Nstate + i])) * dv_adia[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + k];
                                // printf("FORCE(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,dv_adia[i*Nstate*Ndof1*Ndof2+j*Ndof1*Ndof2+k*Ndof2+l]);
                            }
                        } else {
                            for (k = 0; k < Ndof1 * Ndof2; k++) {
                                force[k] -= (0.5 * (xe[i] * xe[j] + pe[i] * pe[j]) - creal(gamma_cv[i * Nstate + j])) * (E_adia[j] - E_adia[i]) * nac[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k];
                                // force[k] -= (0.5 * (xe[i] * xe[j] + pe[i] * pe[j]) - creal(gamma_cv[i * Nstate + j])) * dv_adia[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k];
                            }
                        }
                    }
                }
            }
            break;
        case 1:
            if (rep == 0) {
                if (calforcetype == 1) {
                    for (i = 0; i < Nstate; i++) {
                        for (j = 0; j < Ndof1 * Ndof2; j++) {
                            force[j] -= dV[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + j] * creal(den_e[i * Nstate + i] - gamma_cv[i * Nstate + i]);
                        }
                    }
                } else {
                    for (i = 0; i < Nstate; i++) {
                        for (j = 0; j < Nstate; j++) {
                            for (k = 0; k < Ndof1 * Ndof2; k++) {
                                force[k] -= dV[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k] * creal(den_e[i * Nstate + j] - gamma_cv[i * Nstate + j]);
                            }
                        }
                    }
                }
            } else if (rep == 1) {
                for (i = 0; i < Nstate; i++) {
                    for (j = 0; j < Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < Ndof1 * Ndof2; k++) {
                                force[k] -= creal(den_e[i * Nstate + i] - gamma_cv[i * Nstate + i]) * dv_adia[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + k];
                            }
                        } else {
                            for (k = 0; k < Ndof1 * Ndof2; k++) {
                                force[k] -= creal(den_e[i * Nstate + j] - gamma_cv[i * Nstate + j]) * (E_adia[j] - E_adia[i]) * nac[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k];
                            }
                        }
                    }
                }
            }
            break;
        // case 3:
        //     if (rep == 0) {
        //         if (calforcetype == 1) {
        //             for (i = 0; i < Nstate; i++) {
        //                 for (j = 0; j < Ndof1 * Ndof2; j++) {
        //                     force[j] -= dV[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + j] * creal(den_e4nuc[i * Nstate + i]);
        //                 }
        //             }
        //         } else {
        //             for (i = 0; i < Nstate; i++) {
        //                 for (j = 0; j < Nstate; j++) {
        //                     for (k = 0; k < Ndof1 * Ndof2; k++) {
        //                         force[k] -= dV[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k] * creal(den_e4nuc[i * Nstate + j]);
        //                     }
        //                 }
        //             }
        //         }
        //     } else if (rep == 1) {
        //         for (i = 0; i < Nstate; i++) {
        //             for (j = 0; j < Nstate; j++) {
        //                 if (i == j) {
        //                     for (k = 0; k < Ndof1 * Ndof2; k++) {
        //                         force[k] -= creal(den_e4nuc[i * Nstate + i]) * dv_adia[i * Nstate * Ndof1 * Ndof2 + i * Ndof1 * Ndof2 + k];
        //                     }
        //                 } else {
        //                     for (k = 0; k < Ndof1 * Ndof2; k++) {
        //                         force[k] -= creal(den_e4nuc[i * Nstate + j]) * (E_adia[j] - E_adia[i]) * nac[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k];
        //                     }
        //                 }
        //             }
        //         }
        //     }
        //     break;
    }
}

//  void cal_force_switch(){}

// void cal_force_fssh(){}


// void cal_force_eld(){}


void cal_NACV(){
    int i, j, imax;
    double dotp, norm_nac1, norm_nac2, cosphase;
    double overlap[Nstate * Nstate], overlap2[Nstate * Nstate], vmax;
    double gij[Ndof1 * Ndof2], alpha_BA;
    int id_max[Nstate], idloc, id1, id2;
    double complex eiet[Nstate * Nstate];

    double tempdm1[Nstate*Nstate],tempdm2[Nstate*Nstate],tempdm3[Nstate*Nstate],tempdm4[Nstate*Nstate], tempdv1[Nstate];

    dia_symmat(Nstate, V, E_adia, U_d2a);

    


    if (if_ad_nac) {
        transpose(U_d2a,tempdm1,Nstate);
        dd_matmul(tempdm1,U_ref,overlap,Nstate,Nstate,Nstate);
        // matmul(transpose(U_d2a), U_ref, overlap);
        memset(overlap2, 0, sizeof(overlap2));
       
        for (i = 0; i < Nstate * Nstate; i++){
            tempdm1[i] = fabs(overlap[i]);
        }

        for (i = 0; i < Nstate; i++) {
            idloc=maxloc(tempdm1, Nstate * Nstate);
            id1 = idloc / Nstate; 
            id2 = idloc % Nstate; 
            overlap2[idloc] = (overlap[idloc] >= 0.0) ? 1.0 : -1.0;
            for (j = 0; j < Nstate; j++) {
                tempdm1[id1 * Nstate + j] = 0;
                tempdm1[j * Nstate + id2] = 0;
            }
        }

        dd_matmul(U_d2a, overlap2, tempdm1, Nstate, Nstate, Nstate);
        memcpy(U_d2a,tempdm1,Nstate * Nstate * sizeof(double));

        for (i = 0; i < Nstate * Nstate; i++){
            tempdm2[i] = fabs(overlap2[i]);
        }
        dd_matmul(E_adia, tempdm2, tempdv1, 1, Nstate, Nstate);
        memcpy(E_adia, tempdv1, Nstate * sizeof(double));
    }

    if (type_prop_adia > 0) {
        transpose(U_d2a,tempdm1,Nstate);
        dd_matmul(tempdm1, U_ref, overlap_adia, Nstate, Nstate, Nstate);
    }

    memcpy(U_d2a_old, U_ref, Nstate * Nstate * sizeof(double));
    memcpy(U_ref, U_d2a, Nstate * Nstate * sizeof(double));
    if_ad_nac = 1;


    transpose(U_d2a,tempdm1,Nstate);
    for (i = 0; i < Ndof1; i++) {
        for (j = 0; j < Ndof2; j++) {
            for (int k = 0; k < Nstate * Nstate; k++){
                tempdm2[k] = dV[k * Ndof1 * Ndof2 + i * Ndof2 + j];
            }
            dd_matmul(tempdm1,tempdm2,tempdm3,Nstate,Nstate,Nstate);
            dd_matmul(tempdm3,U_d2a,tempdm4,Nstate,Nstate,Nstate);
            for (int k = 0; k < Nstate * Nstate; k++){
                dv_adia[k * Ndof1 * Ndof2 + i * Ndof2 + j] = tempdm4[k];
            }
        }
    }

    memset(nac, 0, Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double));

    for (i = 0; i < Nstate; i++) {
        for (j = 0; j < Nstate; j++) {
            if (i == j) continue;
            for (int k = 0; k < Ndof1 * Ndof2; k++){
                nac[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k] = dv_adia[i * Nstate * Ndof1 * Ndof2 + j * Ndof1 * Ndof2 + k] / (E_adia[j] - E_adia[i]);
            } 
        }
    }



    // int slavecore_id;
    // slavecore_id = athread_get_id(-1);
    // if(slavecore_id == 0){
    //     for(i = 0; i< Nstate;i++){
    //         printf("%d E=%18.8E\n",i,E_adia[i]);
    //     }
    //     for(i = 0; i< Nstate*Nstate;i++){
    //         printf("%d U=%18.8E\n",i,U_d2a[i]);
    //     }
    //     for(i = 0; i< Ndof1*Ndof2;i++){
    //         printf("%d dvd=%18.8E\n",i,dV[0*Nstate*Ndof1*Ndof2+0*Ndof1*Ndof2+i]);
    //     }
    //     printf("%18.8E\n",nac[0*Nstate*Ndof1*Ndof2+1*Ndof1*Ndof2+0]);
    //     for(i = 0; i< Ndof1*Ndof2;i++){
    //         printf("%d dv1=%18.8E\n",i,dv_adia[0*Nstate*Ndof1*Ndof2+0*Ndof1*Ndof2+i]);
    //     }
    //     for(i = 0; i< Ndof1*Ndof2;i++){
    //         printf("%d dv2=%18.8E\n",i,dv_adia[1*Nstate*Ndof1*Ndof2+1*Ndof1*Ndof2+i]);
    //     }
    //     for(i = 0; i< Ndof1*Ndof2;i++){
    //         printf("%d dv12=%18.8E\n",i,dv_adia[0*Nstate*Ndof1*Ndof2+1*Ndof1*Ndof2+i]);
    //     }
    //     for(i = 0; i< Ndof1*Ndof2;i++){
    //         printf("%d dv21=%18.8E\n",i,dv_adia[1*Nstate*Ndof1*Ndof2+0*Ndof1*Ndof2+i]);
    //     }
    //     exit(-1);
        // for(i = 0; i< Nstate*Nstate;i++){
        //     printf("%d nac=%18.8E\n",i,nac[i*Ndof1*Ndof2+0]);
        // }
    // }

    // if (ifBA == 1) {
    //     memset(nac_BAeff, 0, Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double complex));
    //     for (i = 0; i < Nstate; i++) {
    //         for (j = 0; j < Nstate; j++) {
    //             if (i == j) continue;
    //             if (fabs(tdc_BA[i * Nstate + j]) < 1E-10) {
    //                 memset(&nac[i * Nstate + j], 0, Ndof1 * Ndof2 * sizeof(double complex));
    //             } else {
    //                 gij = dv_adia[i * Nstate + i] - dv_adia[j * Nstate + j];
    //                 alpha_BA = (tdc_BA[i * Nstate + j] - sum(P_nuc, mass, gij)) / sum(P_nuc, mass);
    //                 nac_BAeff[i * Nstate + j] = gij + alpha_BA * P_nuc / mass;
    //                 nac[i * Nstate + j] = nac_BAeff[i * Nstate + j] * tdc_BA[i * Nstate + j] / sum(P_nuc, mass, nac_BAeff[i * Nstate + j]);
    //             }
    //             if (sum(abs(nac_old[i * Nstate + j])) > 1E-10 && sum(nac_old[i * Nstate + j], nac[i * Nstate + j]) < 0) {
    //                 nac[i * Nstate + j] = -nac[i * Nstate + j];
    //             }
    //         }
    //     }
    //     memcpy(nac_old, nac, Nstate * Nstate * Ndof1 * Ndof2 * sizeof(double complex));
    // }

}

// void cal_dvadia_state(){}

void cal_propagator_adia(int Nstate, double dt, double complex *U){
    int i, j;
    double complex H_eff[Nstate * Nstate], C[Nstate * Nstate];
    double E[Nstate], eig[Nstate * Nstate];
    double sineig[Nstate * Nstate], coseig[Nstate * Nstate];
    double real_pro[Nstate * Nstate], img_pro[Nstate * Nstate];
    double complex mean, rho[Nstate * Nstate], E_t[Nstate * Nstate];
    double x00[Nstate], p00[Nstate];
    double P_eff[Ndof1 * Ndof2];
    double E_avg;
    double complex eiet[Nstate * Nstate];
    double sum;
    double complex tempcm1[Nstate * Nstate], tempcm2[Nstate * Nstate];
    

    // for (i = 0; i < Nstate * Nstate; i++) {
    //     H_eff[i] = 0.0 + 0.0 * I;
    // }
    memset(H_eff,0,Nstate*Nstate*sizeof(double complex));

    for (i = 0; i < Nstate; i++) {
        for (j = 0; j < Nstate; j++) {
            sum=0;
            for (int k=0;k<Ndof1*Ndof2;k++){
                sum += P_nuc[k] * nac[i*Nstate*Ndof1*Ndof2+j*Ndof1*Ndof2+k]/mass[k];
            }
            H_eff[i * Nstate + j] = E_adia[i] * Kronecker_delta(i, j) - I * sum;
        }
    }

    dia_hermitemat(Nstate, H_eff, E, C);

  
    memset(sineig,0,Nstate * Nstate * sizeof(double));
    memset(coseig,0,Nstate * Nstate * sizeof(double));

    for (i = 0; i < Nstate; i++) {
        // eig[i * Nstate + i] = E[i];
        sineig[i * Nstate + i] = sin(E[i] * dt);
        coseig[i * Nstate + i] = cos(E[i] * dt);
    }

    // matmul(C, coseig - I * sineig, transpose_conjg(C), U);
    for (i = 0; i < Nstate * Nstate; i++){
        eiet[i] = coseig[i] - I * sineig[i];
    }
    diagger(C,tempcm1,Nstate);
    cc_matmul(eiet,tempcm1,tempcm2,Nstate,Nstate,Nstate);
    cc_matmul(C,tempcm2,U,Nstate,Nstate,Nstate);

    // int slavecore_id = athread_get_id(-1);
    // if(slavecore_id == 10){
    //     for(int i=0;i<Nstate*Nstate;i++){
    //         printf("temp=%18.8E %18.8E\n",creal(tempcm2[i]),cimag(tempcm2[i]));
    //     }
    //     for(int i=0;i<Nstate*Nstate;i++){
    //         printf("H=%18.8E %18.8E\n",creal(H_eff[i]),cimag(H_eff[i]));
    //     }
    //     for(int i=0;i<Nstate*Nstate;i++){
    //         printf("eiet=%18.8E %18.8E\n",creal(eiet[i]),cimag(eiet[i]));
    //     }
    //     for(int i=0;i<Nstate*Nstate;i++){
    //         printf("C=%18.8E %18.8E\n",creal(C[i]),cimag(C[i]));
    //     }
    //     for(int i=0;i<Nstate*Nstate;i++){
    //         printf("CD=%18.8E %18.8E\n",creal(tempcm1[i]),cimag(tempcm1[i]));
    //     }
    //     for(int i=0;i<Nstate*Nstate;i++){
    //         printf("U=%18.8E %18.8E\n",creal(U[i]),cimag(U[i]));
    //     }
    // }       


    if (type_prop_adia == 1) {
        // for (i = 0; i < Nstate * Nstate; i++) {
        //     eiet[i] = 0;
        // }
        memset(eiet,0,Nstate * Nstate * sizeof(double complex));
        for (i = 0; i < Nstate; i++) {
            eiet[i * Nstate + i] = cos(E_adia[i] * dt) - I * sin(E_adia[i] * dt);
        }

        switch (type_algorithm) {
            case 1:
                // matmul(eiet, overlap_adia, U);
                cd_matmul(eiet,overlap_adia,U,Nstate,Nstate,Nstate);
                break;
            // case 2:
            //     matmul(U_d2a, eiet, U);
            //     break;
            // case 3:
            //     matmul(eiet, overlap_adia, U);
            //     break;
            // case 4:
            //     if (type_prop_4cont == 1) {
            //         matmul(U_d2a, eiet, U);
            //     } else if (type_prop_4cont == 2) {
            //         matmul(eiet, transpose(U_d2a), U);
            //     }
            //     break;
        }
    }
}

// void cal_propagator_adia_unsmash(){}
// void cal_propagator_dia_unsmash(){}

