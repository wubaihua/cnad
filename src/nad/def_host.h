#ifndef DEF_HOST_H
#define DEF_HOST_H

#include <complex.h>
#include <stdint.h>
// #include "def.h"
#include "constant.h"
#include "gmath.h"
#include <stdio.h>
#include <math.h>
#include "cJSON.h"
// #include "msmodel.h"
#include <stdbool.h>


// extern int mpi_size, mpi_rank, mpi_ierr;

// // Dynamically allocated arrays
// extern double complex *mpi_den; // 3D complex array: [size1][size2][size3]
// extern double complex *mpi_cfall; // 5D complex array: [size1][size2][size3][size4][size5]
// extern double complex *mpi_cfeff; // 1D complex array: [size1]
// extern double *real_rho; // 3D double array: [size1][size2][size3]
// extern double *imag_rho; // 3D double array: [size1][size2][size3]
// extern double *real_cfall; // 5D double array: [size1][size2][size3][size4][size5]
// extern double *imag_cfall; // 5D double array: [size1][size2][size3][size4][size5]
// extern double *real_cfeff; // 1D double array: [size1]
// extern double *imag_cfeff; // 1D double array: [size1]
// extern double *mpi_population; // 2D double array: [size1][size2]
// extern double *mpi_pop_fb; // 3D double array: [size1][size2][size3]
// extern double *mpi_R_nuc_mean; // 3D double array: [size1][size2][size3]
// extern double *mpi_P_nuc_mean; // 3D double array: [size1][size2][size3]
// extern double *mpi_R2_nuc_mean; // 3D double array: [size1][size2][size3]
// extern double *mpi_P2_nuc_mean; // 3D double array: [size1][size2][size3]
// extern unsigned long long *mpi_N_nan_sum; // 1D int array: [size1]
// extern double *mpi_real_den;
// extern double *mpi_imag_den;
// extern double *mpi_real_cfeff;
// extern double *mpi_imag_cfeff;


// extern double complex *fi_den; // 3D complex array: [size1][size2][size3]
// extern double complex *fi_cfall; // 5D complex array: [size1][size2][size3][size4][size5]
// extern double complex *fi_cfeff; // 1D complex array: [size1]
// extern double *fi_population; // 2D double array: [size1][size2]
// extern double *fi_pop_fb; // 3D double array: [size1][size2][size3]
// extern double *fi_R_nuc_mean; // 3D double array: [size1][size2][size3]
// extern double *fi_P_nuc_mean; // 3D double array: [size1][size2][size3]
// extern double *fi_R2_nuc_mean; // 3D double array: [size1][size2][size3]
// extern double *fi_P2_nuc_mean; // 3D double array: [size1][size2][size3]
// extern unsigned long long *fi_N_nan_sum; // 1D int array: [size1]
// extern double *fi_real_den;
// extern double *fi_imag_den;
// extern double *fi_real_cfeff;
// extern double *fi_imag_cfeff;
// extern double *fi_time_grid;


// // int *count_st; // 2D int array: [size1][size2]
// // int *mpi_count_st; // 2D int array: [size1][size2]
// // int *count_pertraj; // 1D int array: [size1]

// // double *R_nuc_mean; // 3D double array: [size1][size2][size3]
// // double *P_nuc_mean; // 3D double array: [size1][size2][size3]
// // double *R2_nuc_mean; // 3D double array: [size1][size2][size3]
// // double *P2_nuc_mean; // 3D double array: [size1][size2][size3]

// extern char *filepath; // Path of Model input file
// extern char *workpath;

// extern int Ndof1, Ndof2;
// extern int Nstate;
// extern int init_occ_4read;

// // double *xe, *pe; // 1D double arrays: [size1] Meyer-Miller mapping variables
// // double *ye, *pxe, *pye; // 1D double arrays: [size1] Li-Miller mapping variables
// // double complex *ce; // 1D complex array: [size1] Electronic amplitude
// // double complex *correfun_t; // 2D complex array: [size1][size2] Electronic reduced density matrix for evolution
// // double complex *gamma_cv; // 2D complex array: [size1][size2] Commutator matrix or zero-point energy parameters
// // double complex *den_e, *den_e4nuc; // 2D complex arrays: [size1][size2]

// // double *xe_mb, *pe_mb; // 2D double arrays: [size1][size2]
// // int *state_unsmash; // 1D int array: [size1]
// extern int if_typemb;
// // double *gc_unsmash, *gp_unsmash; // 1D double arrays: [size1]

// // double *xe_cv, *pe_cv; // 1D double arrays: [size1]

// // More variables...

// extern int type_evo;
// extern double gamma_zpe;
// extern double sigma2_lsc;
// extern double gamma1_lsc, gamma2_lsc;
// // double identity_M1, identity_M2;
// extern int scheme_cvk, scheme_focus;
// extern double alpha_gdtwa, beta_gdtwa, delta_gdtwa, eps_gdtwa;
// extern int if_alpha;
// // int id_state, id_hop;
// // int id_state_old;
// // double *t_decoh, *t_coh, *L_red; // 1D double arrays: [size1]
// extern double w_dish;

// // double *xe_old, *pe_old; // 1D double arrays: [size1]
// // double *P_nuc_old, *R_nuc_old; // 2D double arrays: [size1][size2]
// // double *force_old; // 2D double array: [size1][size2]
// // double complex *gamma_cv_old, *den_e_old; // 2D complex arrays: [size1][size2]
// // double *nacv_old, *V_old, *dV_old, *E_adia_old, *dv_adia_old; // Various dimensional arrays


// // 动态分配的数组声明
// // double *nac_check_old; // 4D double array: [size1][size2][size3][size4]
// // double *U_ref_old; // 2D double array: [size1][size2]

// // double *P_nuc_old_traj; // 2D double array: [size1][size2]
// // double *R_nuc_old_traj; // 2D double array: [size1][size2]

// // double *deltaR_afssh; // 4D double array: [size1][size2][size3][size4]
// // double *deltaP_afssh; // 4D double array: [size1][size2][size3][size4]
// // double *deltaF_afssh; // 4D double array: [size1][size2][size3][size4]

// extern int index_t0;
// extern int index_t0_1, index_t0_2;
// // double complex correfun_0;

// // double complex *correfun_0_ms2; // 1D complex array: [size1]

// // double complex *cf0; // 2D complex array: [size1][size2]
// // double complex *cfall; // 5D complex array: [size1][size2][size3][size4][size5]
// // double complex *cfeff; // 1D complex array: [size1]
// // double *weight0; // 2D double array: [size1][size2]
// // double *weightt; // 2D double array: [size1][size2]
// extern int if_allcf, allcf_times;

// extern double complex correfun_0_pldm1, correfun_0_pldm2;
// // double complex *prop_pldm; // 2D complex array: [size1][size2]

// // double *U_d2a; // 2D double array: [size1][size2]
// // double *E_adia; // 1D double array: [size1]
// // double *nac; // 4D double array: [size1][size2][size3][size4]
// // double *dv_adia; // 4D double array: [size1][size2][size3][size4]
// // double *P_kin; // 2D double array: [size1][size2]
// // double *nac_check; // 4D double array: [size1][size2][size3][size4]
// // double *U_ref; // 2D double array: [size1][size2]
// // double *overlap_adia; // 2D double array: [size1][size2]
// extern int rep; // 0 for diabatic, 1 for adiabatic

// // double *U_d2a_old; // 2D double array: [size1][size2]

// extern int ifcv; // 0: no adjustment; -1: adjustment without evolution; 1: cv adjustment 
// extern int ifid;
// // double complex *den; // 3D complex array: [size1][size2][size3] (total electronic reduced density matrix)
// // double *population; // 2D double array: [size1][size2]
// // double *pop_fb; // 3D double array: [size1][size2][size3]

// // double complex *den_traj; // 3D complex array: [size1][size2][size3]
// // double complex *den2_traj; // 3D complex array: [size1][size2][size3]
// // double *population_traj; // 2D double array: [size1][size2]
// // double *population2_traj; // 2D double array: [size1][size2]
// // double *pop_fb_traj; // 3D double array: [size1][size2][size3]
// // double *pop_fb2_traj; // 3D double array: [size1][size2][size3]

// extern int if_st_fb;

// // double *V; // 2D double array: [size1][size2]
// // double *dV; // 4D double array: [size1][size2][size3][size4]
// // double *V_ref; // 3D double array: [size1][size2][size3]
// // double *dV_ref; // 5D double array: [size1][size2][size3][size4][size5]
// // double complex *propagator; // 2D complex array: [size1][size2]
// // double complex *propagator_path; // 2D complex array: [size1][size2]
// // double complex *propagator_switch; // 2D complex array: [size1][size2]

// // double *R_nuc; // 2D double array: [size1][size2]
// // double *P_nuc; // 2D double array: [size1][size2]
// // double *mass; // 2D double array: [size1][size2]

// // double *force; // 2D double array: [size1][size2]
// // double *force_nuc; // 2D double array: [size1][size2]
// // double *force_ref; // 3D double array: [size1][size2][size3]
// // double *force_nuc_ref; // 3D double array: [size1][size2][size3]
// extern int type_traj_sed;

// extern double beta, temperature;

// extern char method[20];

// extern double dt, ttot;
// // double *timegrid; // 1D double array: [size1]
// extern long long Nbreak, Ngrid;
// extern long long Ntraj;

// extern char unit_t[20];
// extern double unittrans_t;

// extern int outputtype;

// extern int calforcetype;

// extern int ifoutputmpi;

// extern int sampletype;

// extern int if_st_nan;
// // unsigned long long *N_nan_sum; // 1D int array: [size1]

// extern int if_traj;

// extern int type_phase;

// extern int type_ad_fssh;

// extern int if_ref;

// extern int if_1st;

// extern int if_inv_focus;

// extern int if_Pdis, s_N;
// extern double s_start, s_end;
// // double *s; // 1D double array: [size1]
// // double *real_expisP; // 1D double array: [size1]
// // double *imag_expisP; // 1D double array: [size1]
// // double complex *expisP; // 1D complex array: [size1]
// // double complex *mpi_expisP; // 1D complex array: [size1]
// // double *mpi_real_expisP;
// // double *mpi_imag_expisP;

// // double *R_nuc_ref; // 3D double array: [size1][size2][size3]
// // double *P_nuc_ref; // 3D double array: [size1][size2][size3]
// // double *xe_ref; // 2D double array: [size1][size2]
// // double *pe_ref; // 2D double array: [size1][size2]
// // double *ye_ref; // 2D double array: [size1][size2]
// // double *pxe_ref; // 2D double array: [size1][size2]
// // double *pye_ref; // 2D double array: [size1][size2]
// // double complex *gamma_cv_ref; // 3D complex array: [size1][size2][size3]
// // double complex *den_e_ref; // 3D complex array: [size1][size2][size3]
// // double complex *inverse_kernel; // 2D complex array: [size1][size2]

// extern int Nref;

// // double *pop0; // 1D double array: [size1]

// extern bool if_ad_nac;

// extern int mean_nuc;

// extern bool if_occ_rand;

// extern int if_engconsv;
// // double complex *engconsv_adjmat; // 2D complex array: [size1][size2]

// extern int if_RBC;

// // double measure_mash;
// extern int type_mash;
// // double complex rho0_mash[4], rhot_mash[4];

// // double U0_mash[4], mea_mat_mash[4];

// // double complex *rho0_unsmash; // 2D complex array: [size1][size2]
// // double complex *rhot_unsmash; // 2D complex array: [size1][size2]
// // double *U0_unsmash; // 2D double array: [size1][size2]
// // double *mea_mat_unsmash; // 2D double array: [size1][size2]

// extern int ifBA;
// // double *E_adia_old_traj; // 1D double array: [size1]
// // double *nac_BAeff; // 4D double array: [size1][size2][size3][size4]
// // double *tdc_BA; // 2D double array: [size1][size2]
// // double *nac_old; // 4D double array: [size1][size2][size3][size4]
// // double *P_nuc_BA_old; // 2D double array: [size1][size2]

// // int iflbd, occ_lbd;
// // double *rate_lbd; // 1D double array: [size1]
// // double rate_para_lbd;

// // int ifmfsh, occ_mfsh;
// // double thres;

// // int occ_rdmfocus;

// // int ifrw;
// // double rwfactor, beta_rw, temperature_rw;

// extern int ifmodprop;
// // double complex *propagator_ref; // 2D complex array: [size1][size2]

// extern int ifcorreden;

// // int memorylength;
// // int itime_save, i_re_save;

// extern int if_st_eng;
// // double *energy_est; // 1D double array: [size1]
// // double *mpi_energy_est; // 1D double array: [size1]

// extern int typeevo_ele;

// // double *A_jump; // 2D double array: [size1][size2]
// // double *lambda_jump; // 1D double array: [size1]
// // double *U_jump; // 2D double array: [size1][size2]
// // double *alpha_jump; // 1D double array: [size1]
// // int type_jump;

// // int if_switchcv, occ_switchcv, max_switchcv;

// // double Pn0_mf2, Wn0_mf2, energy_mf2;
// // int if_inv_evo, if_BO_mf2, id_max_mf2sh, type_eom_mf2sh, id_init_occ_mf2, if_bak_mf2cv;
// // double *da_mf2; // 1D double array: [size1]
// // double *R_nuc_old_mf2cv; // 2D double array: [size1][size2]
// // double *P_nuc_old_mf2cv; // 2D double array: [size1][size2]
// // double *xe_old_mf2cv; // 1D double array: [size1]
// // double *pe_old_mf2cv; // 1D double array: [size1]
// // double *dE_old_mf2cv; // 1D double array: [size1]
// // double *dE_mf2cv; // 1D double array: [size1]

// extern int type_prop_adia;

// // double complex *G_xpconfg; // 2D complex array: [size1][size2]
// // double *permutation_ms2; // 3D double array: [size1][size2][size3]
// // double complex *G_ms2; // 2D complex array: [size1][size2]
// // double thres_ms2;
// // int *index_ms2; // 1D int array: [size1]
// // int *index_old_ms2; // 1D int array: [size1]

// // double *R_nuc_state; // 2D double array: [size1][size2]
// // double *P_nuc_state; // 2D double array: [size1][size2]
// // double *V_state; // 2D double array: [size1][size2]
// // double *dv_state; // 4D double array: [size1][size2][size3][size4]
// // double *dv_adia_state; // 4D double array: [size1][size2][size3][size4]
// // double *U_d2a_state; // 2D double array: [size1][size2]
// // double *U_ref_state; // 2D double array: [size1][size2]
// // double *force_nuc_state; // 2D double array: [size1][size2]
// // int if_statetraj;

// // bool ifBC_BCMF;

// // int scheme_cvsh;

// extern int ifswitchforce, ifmashforce;

// extern int ifscaleenergy;

// extern double gamma_rescale;
// extern int ifscalegamma;

// // double E_conserve;

// // int ifmsbranch, type_traj_msbranch;
// // int itime_start_msbranch, itime_end_msbranch;
// // int i_re_start_msbranch, igrid_start_msbranch;
// // int iwrong_msbranch;
// // double time_start_msbranch, time_end_msbranch;
// // double *R_nuc_brapoint; // 2D double array: [size1][size2]
// // double *P_nuc_brapoint; // 2D double array: [size1][size2]
// // double *xe_brapoint; // 1D double array: [size1]
// // double *pe_brapoint; // 1D double array: [size1]
// // double complex *gamma_cv_brapoint; // 2D complex array: [size1][size2]
// // double complex *den_e_brapoint; // 2D complex array: [size1][size2]

// // double complex *correfun_t_oldtraj; // 3D complex array: [size1][size2][size3]
// // double *R_nuc_oldtraj; // 3D double array: [size1][size2][size3]
// // double *P_nuc_oldtraj; // 3D double array: [size1][size2][size3]

// extern int direc_padj;

// extern int ifcount;

// extern int ifreflp;

// extern int ifreflp_mash;

// extern int ifhardwall;

// // double *eig_cv; // 1D double array: [size1]
// // double *eig_cv_mat; // 2D double array: [size1][size2]
// // double complex *commu_vari; // 2D complex array: [size1][size2]

// extern int ifzpecorr;

// extern int iflangevin;
// extern double eta_langevin;

// extern double scale_sqc2;

// extern int type_algorithm;

// extern int type_prop_4cont;

// extern int scaleenergy_type;

// extern int n_step_algo5;

// extern int allow_hop;

// extern int if_traceless_force;


// extern int forcetype;
// extern char msmodelname[200];


struct set_host
{
    int mpi_size, mpi_rank, mpi_ierr;

    // Dynamically allocated arrays
    double complex *mpi_den; // 3D complex array: [size1][size2][size3]
    double complex *mpi_cfall; // 5D complex array: [size1][size2][size3][size4][size5]
    double complex *mpi_cfeff; // 1D complex array: [size1]
    double *real_rho; // 3D double array: [size1][size2][size3]
    double *imag_rho; // 3D double array: [size1][size2][size3]
    double *real_cfall; // 5D double array: [size1][size2][size3][size4][size5]
    double *imag_cfall; // 5D double array: [size1][size2][size3][size4][size5]
    double *real_cfeff; // 1D double array: [size1]
    double *imag_cfeff; // 1D double array: [size1]
    double *mpi_population; // 2D double array: [size1][size2]
    double *save_population;
    double *mpi_pop_fb; // 3D double array: [size1][size2][size3]
    double *save_pop_fb;
    double *mpi_R_nuc_mean; // 3D double array: [size1][size2][size3]
    double *mpi_P_nuc_mean; // 3D double array: [size1][size2][size3]
    double *mpi_R2_nuc_mean; // 3D double array: [size1][size2][size3]
    double *mpi_P2_nuc_mean; // 3D double array: [size1][size2][size3]
    double *save_R_nuc_mean; // 3D double array: [size1][size2][size3]
    double *save_P_nuc_mean; // 3D double array: [size1][size2][size3]
    double *save_R2_nuc_mean; // 3D double array: [size1][size2][size3]
    double *save_P2_nuc_mean; // 3D double array: [size1][size2][size3]
    unsigned long long *mpi_N_nan_sum; // 1D int array: [size1]
    unsigned long long *save_N_nan_sum; 
    double *mpi_real_den;
    double *mpi_imag_den;
    double *save_real_den;
    double *save_imag_den;
    double *mpi_real_cfeff;
    double *mpi_imag_cfeff;
    double *save_real_cfeff;
    double *save_imag_cfeff;




    double complex *fi_den; // 3D complex array: [size1][size2][size3]
    double complex *fi_cfall; // 5D complex array: [size1][size2][size3][size4][size5]
    double complex *fi_cfeff; // 1D complex array: [size1]
    double *fi_population; // 2D double array: [size1][size2]
    double *fi_pop_fb; // 3D double array: [size1][size2][size3]
    double *fi_R_nuc_mean; // 3D double array: [size1][size2][size3]
    double *fi_P_nuc_mean; // 3D double array: [size1][size2][size3]
    double *fi_R2_nuc_mean; // 3D double array: [size1][size2][size3]
    double *fi_P2_nuc_mean; // 3D double array: [size1][size2][size3]
    unsigned long long *fi_N_nan_sum; // 1D int array: [size1]
    double *fi_real_den;
    double *fi_imag_den;
    double *fi_real_cfeff;
    double *fi_imag_cfeff;
    double *fi_time_grid;


    // int *count_st; // 2D int array: [size1][size2]
    // int *mpi_count_st; // 2D int array: [size1][size2]
    // int *count_pertraj; // 1D int array: [size1]

    // double *R_nuc_mean; // 3D double array: [size1][size2][size3]
    // double *P_nuc_mean; // 3D double array: [size1][size2][size3]
    // double *R2_nuc_mean; // 3D double array: [size1][size2][size3]
    // double *P2_nuc_mean; // 3D double array: [size1][size2][size3]

    // char *filepath; // Path of Model input file
    char filepath[300];
    char *workpath;

    int Ndof1, Ndof2;
    int Nstate;
    int init_occ_4read;

    // double *xe, *pe; // 1D double arrays: [size1] Meyer-Miller mapping variables
    // double *ye, *pxe, *pye; // 1D double arrays: [size1] Li-Miller mapping variables
    // double complex *ce; // 1D complex array: [size1] Electronic amplitude
    // double complex *correfun_t; // 2D complex array: [size1][size2] Electronic reduced density matrix for evolution
    // double complex *gamma_cv; // 2D complex array: [size1][size2] Commutator matrix or zero-point energy parameters
    // double complex *den_e, *den_e4nuc; // 2D complex arrays: [size1][size2]

    // double *xe_mb, *pe_mb; // 2D double arrays: [size1][size2]
    // int *state_unsmash; // 1D int array: [size1]
    int if_typemb;
    // double *gc_unsmash, *gp_unsmash; // 1D double arrays: [size1]

    // double *xe_cv, *pe_cv; // 1D double arrays: [size1]

    // More variables...

    int type_evo;
    double gamma_zpe;
    double sigma2_lsc;
    double gamma1_lsc, gamma2_lsc;
    // double identity_M1, identity_M2;
    int scheme_cvk, scheme_focus;
    double alpha_gdtwa, beta_gdtwa, delta_gdtwa, eps_gdtwa;
    int if_alpha;
    // int id_state, id_hop;
    // int id_state_old;
    // double *t_decoh, *t_coh, *L_red; // 1D double arrays: [size1]
    double w_dish;

    // double *xe_old, *pe_old; // 1D double arrays: [size1]
    // double *P_nuc_old, *R_nuc_old; // 2D double arrays: [size1][size2]
    // double *force_old; // 2D double array: [size1][size2]
    // double complex *gamma_cv_old, *den_e_old; // 2D complex arrays: [size1][size2]
    // double *nacv_old, *V_old, *dV_old, *E_adia_old, *dv_adia_old; // Various dimensional arrays


    // 动态分配的数组声明
    // double *nac_check_old; // 4D double array: [size1][size2][size3][size4]
    // double *U_ref_old; // 2D double array: [size1][size2]

    // double *P_nuc_old_traj; // 2D double array: [size1][size2]
    // double *R_nuc_old_traj; // 2D double array: [size1][size2]

    // double *deltaR_afssh; // 4D double array: [size1][size2][size3][size4]
    // double *deltaP_afssh; // 4D double array: [size1][size2][size3][size4]
    // double *deltaF_afssh; // 4D double array: [size1][size2][size3][size4]

    int index_t0;
    int index_t0_1, index_t0_2;
    // double complex correfun_0;

    // double complex *correfun_0_ms2; // 1D complex array: [size1]

    // double complex *cf0; // 2D complex array: [size1][size2]
    // double complex *cfall; // 5D complex array: [size1][size2][size3][size4][size5]
    // double complex *cfeff; // 1D complex array: [size1]
    // double *weight0; // 2D double array: [size1][size2]
    // double *weightt; // 2D double array: [size1][size2]
    int if_allcf, allcf_times;

    double complex correfun_0_pldm1, correfun_0_pldm2;
    // double complex *prop_pldm; // 2D complex array: [size1][size2]

    // double *U_d2a; // 2D double array: [size1][size2]
    // double *E_adia; // 1D double array: [size1]
    // double *nac; // 4D double array: [size1][size2][size3][size4]
    // double *dv_adia; // 4D double array: [size1][size2][size3][size4]
    // double *P_kin; // 2D double array: [size1][size2]
    // double *nac_check; // 4D double array: [size1][size2][size3][size4]
    // double *U_ref; // 2D double array: [size1][size2]
    // double *overlap_adia; // 2D double array: [size1][size2]
    int rep; // 0 for diabatic, 1 for adiabatic

    // double *U_d2a_old; // 2D double array: [size1][size2]

    int ifcv; // 0: no adjustment; -1: adjustment without evolution; 1: cv adjustment 
    int ifid;
    // double complex *den; // 3D complex array: [size1][size2][size3] (total electronic reduced density matrix)
    // double *population; // 2D double array: [size1][size2]
    // double *pop_fb; // 3D double array: [size1][size2][size3]

    // double complex *den_traj; // 3D complex array: [size1][size2][size3]
    // double complex *den2_traj; // 3D complex array: [size1][size2][size3]
    // double *population_traj; // 2D double array: [size1][size2]
    // double *population2_traj; // 2D double array: [size1][size2]
    // double *pop_fb_traj; // 3D double array: [size1][size2][size3]
    // double *pop_fb2_traj; // 3D double array: [size1][size2][size3]

    int if_st_fb;

    // double *V; // 2D double array: [size1][size2]
    // double *dV; // 4D double array: [size1][size2][size3][size4]
    // double *V_ref; // 3D double array: [size1][size2][size3]
    // double *dV_ref; // 5D double array: [size1][size2][size3][size4][size5]
    // double complex *propagator; // 2D complex array: [size1][size2]
    // double complex *propagator_path; // 2D complex array: [size1][size2]
    // double complex *propagator_switch; // 2D complex array: [size1][size2]

    // double *R_nuc; // 2D double array: [size1][size2]
    // double *P_nuc; // 2D double array: [size1][size2]
    // double *mass; // 2D double array: [size1][size2]

    // double *force; // 2D double array: [size1][size2]
    // double *force_nuc; // 2D double array: [size1][size2]
    // double *force_ref; // 3D double array: [size1][size2][size3]
    // double *force_nuc_ref; // 3D double array: [size1][size2][size3]
    int type_traj_sed;

    double beta, temperature;

    char method[20];

    double dt, ttot;
    // double *timegrid; // 1D double array: [size1]
    long long Nbreak, Ngrid;
    long long Ntraj;

    char unit_t[20];
    double unittrans_t;

    int outputtype;

    int calforcetype;

    int ifoutputmpi;

    int sampletype;

    int if_st_nan;
    // unsigned long long *N_nan_sum; // 1D int array: [size1]

    int if_traj;

    int type_phase;

    int type_ad_fssh;

    int if_ref;

    int if_1st;

    int if_inv_focus;

    int if_Pdis, s_N;
    double s_start, s_end;
    // double *s; // 1D double array: [size1]
    // double *real_expisP; // 1D double array: [size1]
    // double *imag_expisP; // 1D double array: [size1]
    // double complex *expisP; // 1D complex array: [size1]
    // double complex *mpi_expisP; // 1D complex array: [size1]
    // double *mpi_real_expisP;
    // double *mpi_imag_expisP;

    // double *R_nuc_ref; // 3D double array: [size1][size2][size3]
    // double *P_nuc_ref; // 3D double array: [size1][size2][size3]
    // double *xe_ref; // 2D double array: [size1][size2]
    // double *pe_ref; // 2D double array: [size1][size2]
    // double *ye_ref; // 2D double array: [size1][size2]
    // double *pxe_ref; // 2D double array: [size1][size2]
    // double *pye_ref; // 2D double array: [size1][size2]
    // double complex *gamma_cv_ref; // 3D complex array: [size1][size2][size3]
    // double complex *den_e_ref; // 3D complex array: [size1][size2][size3]
    // double complex *inverse_kernel; // 2D complex array: [size1][size2]

    int Nref;

    // double *pop0; // 1D double array: [size1]

    // bool if_ad_nac;

    int mean_nuc;

    bool if_occ_rand;

    int if_engconsv;
    // double complex *engconsv_adjmat; // 2D complex array: [size1][size2]

    int if_RBC;

    // double measure_mash;
    int type_mash;
    // double complex rho0_mash[4], rhot_mash[4];

    // double U0_mash[4], mea_mat_mash[4];

    // double complex *rho0_unsmash; // 2D complex array: [size1][size2]
    // double complex *rhot_unsmash; // 2D complex array: [size1][size2]
    // double *U0_unsmash; // 2D double array: [size1][size2]
    // double *mea_mat_unsmash; // 2D double array: [size1][size2]

    int ifBA;
    // double *E_adia_old_traj; // 1D double array: [size1]
    // double *nac_BAeff; // 4D double array: [size1][size2][size3][size4]
    // double *tdc_BA; // 2D double array: [size1][size2]
    // double *nac_old; // 4D double array: [size1][size2][size3][size4]
    // double *P_nuc_BA_old; // 2D double array: [size1][size2]

    // int iflbd, occ_lbd;
    // double *rate_lbd; // 1D double array: [size1]
    // double rate_para_lbd;

    // int ifmfsh, occ_mfsh;
    // double thres;

    // int occ_rdmfocus;

    // int ifrw;
    // double rwfactor, beta_rw, temperature_rw;

    int ifmodprop;
    // double complex *propagator_ref; // 2D complex array: [size1][size2]

    int ifcorreden;

    // int memorylength;
    // int itime_save, i_re_save;

    int if_st_eng;
    // double *energy_est; // 1D double array: [size1]
    // double *mpi_energy_est; // 1D double array: [size1]

    int typeevo_ele;

    // double *A_jump; // 2D double array: [size1][size2]
    // double *lambda_jump; // 1D double array: [size1]
    // double *U_jump; // 2D double array: [size1][size2]
    // double *alpha_jump; // 1D double array: [size1]
    // int type_jump;

    // int if_switchcv, occ_switchcv, max_switchcv;

    // double Pn0_mf2, Wn0_mf2, energy_mf2;
    // int if_inv_evo, if_BO_mf2, id_max_mf2sh, type_eom_mf2sh, id_init_occ_mf2, if_bak_mf2cv;
    // double *da_mf2; // 1D double array: [size1]
    // double *R_nuc_old_mf2cv; // 2D double array: [size1][size2]
    // double *P_nuc_old_mf2cv; // 2D double array: [size1][size2]
    // double *xe_old_mf2cv; // 1D double array: [size1]
    // double *pe_old_mf2cv; // 1D double array: [size1]
    // double *dE_old_mf2cv; // 1D double array: [size1]
    // double *dE_mf2cv; // 1D double array: [size1]

    int type_prop_adia;

    // double complex *G_xpconfg; // 2D complex array: [size1][size2]
    // double *permutation_ms2; // 3D double array: [size1][size2][size3]
    // double complex *G_ms2; // 2D complex array: [size1][size2]
    // double thres_ms2;
    // int *index_ms2; // 1D int array: [size1]
    // int *index_old_ms2; // 1D int array: [size1]

    // double *R_nuc_state; // 2D double array: [size1][size2]
    // double *P_nuc_state; // 2D double array: [size1][size2]
    // double *V_state; // 2D double array: [size1][size2]
    // double *dv_state; // 4D double array: [size1][size2][size3][size4]
    // double *dv_adia_state; // 4D double array: [size1][size2][size3][size4]
    // double *U_d2a_state; // 2D double array: [size1][size2]
    // double *U_ref_state; // 2D double array: [size1][size2]
    // double *force_nuc_state; // 2D double array: [size1][size2]
    // int if_statetraj;

    // bool ifBC_BCMF;

    // int scheme_cvsh;

    int ifswitchforce, ifmashforce;

    int ifscaleenergy;

    double gamma_rescale;
    int ifscalegamma;


    int if_scale_sqc;
    // double E_conserve;

    // int ifmsbranch, type_traj_msbranch;
    // int itime_start_msbranch, itime_end_msbranch;
    // int i_re_start_msbranch, igrid_start_msbranch;
    // int iwrong_msbranch;
    // double time_start_msbranch, time_end_msbranch;
    // double *R_nuc_brapoint; // 2D double array: [size1][size2]
    // double *P_nuc_brapoint; // 2D double array: [size1][size2]
    // double *xe_brapoint; // 1D double array: [size1]
    // double *pe_brapoint; // 1D double array: [size1]
    // double complex *gamma_cv_brapoint; // 2D complex array: [size1][size2]
    // double complex *den_e_brapoint; // 2D complex array: [size1][size2]

    // double complex *correfun_t_oldtraj; // 3D complex array: [size1][size2][size3]
    // double *R_nuc_oldtraj; // 3D double array: [size1][size2][size3]
    // double *P_nuc_oldtraj; // 3D double array: [size1][size2][size3]

    int direc_padj;

    int ifcount;

    int ifreflp;

    int ifreflp_mash;

    int ifhardwall;

    // double *eig_cv; // 1D double array: [size1]
    // double *eig_cv_mat; // 2D double array: [size1][size2]
    // double complex *commu_vari; // 2D complex array: [size1][size2]

    int ifzpecorr;

    int iflangevin;
    double eta_langevin;

    // double scale_sqc2;

    int type_algorithm;

    int type_prop_4cont;

    int scaleenergy_type;

    int n_step_algo5;

    int allow_hop;

    int if_traceless_force;


    int nproc_sw;



    int type_hop;

    int if_eld;



    int if_default;


    int if_classical;



    // model parameters


        int forcetype;
        char msmodelname[200];


        // Spin-Boson Model parameters
        int N_bath_SBM, bathtype; // bathtype=1 for Ohmic; bathtype=2 for Debye
        double eps_SBM, delta_SBM, alpha_SBM, omega_c_SBM, lambda_SBM, s_SBM;
        double *c_SBM, *omega_SBM;


        // 3-state Morse potential models parameters
        int type_morse3;
        double A_morse3[9], beta_morse3[3], D_morse3[3], C_morse3[3], R0_morse3[3];
        double R00_morse3[9], alpha_morse3[9];
        double Re_morse3;
        double omega_morse3, mass_morse3;


        // site-exciton model parameters
        int N_bath_SEM, Nstate_SEM;
        double bias_SEM, delta_SEM, k_eff_SEM, omega_c_SEM, lambda_SEM;
        double *c_SEM, *omega_SEM, *H_ele_SEM;


        // atom-in-cavity model parameters
        int N_mode_aic, Nstate_aic;
        double L_aic;
        double *eps_aic, *miu_aic, *omega_aic, *lambda_aic;


        // Linear vibronic coupling model parameters
        int N_mode_lvcm, Nstate_lvcm;
        double *H_ele_lvcm, *c_lvcm, *omega_lvcm ;
        double *R0_lvcm, *P0_lvcm, *alpha_lvcm;


        // pyrazine LVCM
        int type_pyrazine;


        // crco5 LVCM
        int type_crco5;

        // DNA LVCM
        int type_dnalvcm;


        // Tully model parameters
        int type_tully;
        double A_tully, B_tully, C_tully, D_tully, E_tully, A2_tully, f_tully, Z_tully;
        double R0_tully, P0_tully, gamma_tully;


        // SEM with dipole moment parameters:
        int type_cal_SEMdp; // 1 for PT, 2 for NP
        int N_bath_SEMdp, Nstate_SEMdp;
        double bias_SEMdp, delta_SEMdp, k_eff_SEMdp, omega_c_SEMdp, lambda_SEMdp;
        double *c_SEMdp, *omega_SEMdp, *H_ele_SEMdp;
        double *dipole_SEMdp;



        // rubrene parameters:
        int Nstate_rubrene, N_mode_rubrene, Nmole_rubrene;
        double Vc_rubrene, factor_freq_rubrene;
        double *omega_rubrene, *lambda_rubrene, *g_rubrene, *mass_rubrene;


        // frozen nuclei limit
        int Nstate_frozen;
        double *Hele_frozen;
        


        // dual harmonic oscilators:
        int type_dualho;
        // double A1_dualho, A2_dualho, R01_dualho, R02_dualho, 

};



void initial_para(struct set_host *seth);
void readinp(struct set_host *seth);
void init_host(struct set_host *seth);

void print_info(struct set_host *seth);

void fileout(struct set_host *seth) ;

void fileout_mpi(int id, struct set_host *seth);

#endif // DEF_HOST_H