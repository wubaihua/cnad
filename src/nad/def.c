

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
#ifdef sunway
    #include "iofun_slave.h"
    #include <slave.h>
#endif

// int mpi_size, mpi_rank, mpi_ierr;

// Dynamically allocated arrays
// double complex *mpi_sets->den; // 3D complex array: [size1][size2][size3]
// double complex *mpi_cfall; // 5D complex array: [size1][size2][size3][size4][size5]
// double complex *mpi_sets->cfeff; // 1D complex array: [size1]
// double *real_rho; // 3D double array: [size1][size2][size3]
// double *imag_rho; // 3D double array: [size1][size2][size3]
// double *real_cfall; // 5D double array: [size1][size2][size3][size4][size5]
// double *imag_cfall; // 5D double array: [size1][size2][size3][size4][size5]
// double *real_sets->cfeff; // 1D double array: [size1]
// double *imag_sets->cfeff; // 1D double array: [size1]
// double *mpi_sets->population; // 2D double array: [size1][size2]
// double *mpi_sets->pop_fb; // 3D double array: [size1][size2][size3]
// double *mpi_sets->R_nuc_mean; // 3D double array: [size1][size2][size3]
// double *mpi_sets->P_nuc_mean; // 3D double array: [size1][size2][size3]
// double *mpi_sets->R2_nuc_mean; // 3D double array: [size1][size2][size3]
// double *mpi_sets->P2_nuc_mean; // 3D double array: [size1][size2][size3]
// unsigned long long *mpi_sets->N_nan_sum; // 1D int array: [size1]
// double *mpi_real_sets->den;
// double *mpi_imag_sets->den;
// double *mpi_real_sets->cfeff;
// double *mpi_imag_sets->cfeff;


// int *count_st; // 2D int array: [size1][size2]
// int *mpi_count_st; // 2D int array: [size1][size2]
// int *count_sets->pertraj; // 1D int array: [size1]

// double *sets->R_nuc_mean; // 3D double array: [size1][size2][size3]
// double *sets->P_nuc_mean; // 3D double array: [size1][size2][size3]
// double *sets->R2_nuc_mean; // 3D double array: [size1][size2][size3]
// double *sets->P2_nuc_mean; // 3D double array: [size1][size2][size3]

// // char *filepath; // Path of Model input file
// // char *workpath;

// int seth->Ndof1, seth->Ndof2;
// int seth->Nstate;
// int sets->init_occ;

// double *sets->xe, *sets->pe; // 1D double arrays: [size1] Meyer-Miller mapping variables
// double *ye, *psets->xe, *pye; // 1D double arrays: [size1] Li-Miller mapping variables
// double complex *ce; // 1D complex array: [size1] Electronic amplitude
// double complex *sets->correfun_t; // 2D complex array: [size1][size2] Electronic reduced sets->density matrix for evolution
// double complex *sets->gamma_cv; // 2D complex array: [size1][size2] Commutator matrix or zero-point energy parameters
// double complex *sets->den_e, *sets->den_e4nuc; // 2D complex arrays: [size1][size2]

// double *sets->xe_mb, *sets->pe_mb; // 2D double arrays: [size1][size2]
// int *state_unsmash; // 1D int array: [size1]
// // int if_tysets->pemb;
// double *gc_unsmash, *gp_unsmash; // 1D double arrays: [size1]

// double *sets->xe_cv, *sets->pe_cv; // 1D double arrays: [size1]

// // More variables...

// // int seth->type_evo;
// // double seth->gamma_zpe;
// // double seth->sigma2_lsc;
// // double gamma1_lsc, gamma2_lsc;
// double isets->dentity_M1, isets->dentity_M2;
// // int scheme_cvk, scheme_focus;
// // double alpha_gdtwa, seth->beta_gdtwa, delta_gdtwa, eps_gdtwa;
// // int if_alpha;
// int sets->id_state, id_hop;
// int sets->id_state_old;
// double *t_decoh, *t_coh, *L_red; // 1D double arrays: [size1]
// // double w_dish;

// double *sets->xe_old, *sets->pe_old; // 1D double arrays: [size1]
// double *sets->P_nuc_old, *sets->R_nuc_old; // 2D double arrays: [size1][size2]
// double *sets->force_old; // 2D double array: [size1][size2]
// double complex *sets->gamma_cv_old, *sets->den_e_old; // 2D complex arrays: [size1][size2]
// double *sets->nacv_old, *V_old, *sets->dV_old, *sets->E_adia_old, *sets->dv_adia_old; // Various dimensional arrays


// // 动态分配的数组声明
// double *sets->nac_check_old; // 4D double array: [size1][size2][size3][size4]
// double *sets->U_ref_old; // 2D double array: [size1][size2]

// double *sets->P_nuc_old_traj; // 2D double array: [size1][size2]
// double *sets->R_nuc_old_traj; // 2D double array: [size1][size2]

// double *deltaR_afssh; // 4D double array: [size1][size2][size3][size4]
// double *deltaP_afssh; // 4D double array: [size1][size2][size3][size4]
// double *deltaF_afssh; // 4D double array: [size1][size2][size3][size4]

// // int seth->index_t0;
// // int seth->index_t0_1, seth->index_t0_2;
// double complex sets->correfun_0;

// double complex *sets->correfun_0_ms2; // 1D complex array: [size1]

// double complex *sets->cf0; // 2D complex array: [size1][size2]
// double complex *cfall; // 5D complex array: [size1][size2][size3][size4][size5]
// double complex *sets->cfeff; // 1D complex array: [size1]
// double *weight0; // 2D double array: [size1][size2]
// double *weightt; // 2D double array: [size1][size2]
// // int seth->if_allcf, allcf_times;

// // double complex sets->correfun_0_pldm1, sets->correfun_0_pldm2;
// double complex *prop_pldm; // 2D complex array: [size1][size2]

// double *sets->U_d2a; // 2D double array: [size1][size2]
// double *sets->E_adia; // 1D double array: [size1]
// double *sets->nac; // 4D double array: [size1][size2][size3][size4]
// double *sets->dv_adia; // 4D double array: [size1][size2][size3][size4]
// double *P_kin; // 2D double array: [size1][size2]
// double *sets->nac_check; // 4D double array: [size1][size2][size3][size4]
// double *sets->U_ref; // 2D double array: [size1][size2]
// double *sets->overlap_adia; // 2D double array: [size1][size2]
// // int seth->rep; // 0 for diabatic, 1 for adiabatic

// double *sets->U_d2a_old; // 2D double array: [size1][size2]

// // int seth->ifcv; // 0: no adjustment; -1: adjustment without evolution; 1: cv adjustment 
// // int seth->ifid;
// double complex *sets->den; // 3D complex array: [size1][size2][size3] (total electronic reduced sets->density matrix)
// double *sets->population; // 2D double array: [size1][size2]
// double *sets->pop_fb; // 3D double array: [size1][size2][size3]

// double complex *sets->den_traj; // 3D complex array: [size1][size2][size3]
// double complex *sets->den2_traj; // 3D complex array: [size1][size2][size3]
// double *sets->population_traj; // 2D double array: [size1][size2]
// double *sets->population2_traj; // 2D double array: [size1][size2]
// double *sets->pop_fb_traj; // 3D double array: [size1][size2][size3]
// double *sets->pop_fb2_traj; // 3D double array: [size1][size2][size3]

// // int seth->if_st_fb;

// double *V; // 2D double array: [size1][size2]
// double *sets->dV; // 4D double array: [size1][size2][size3][size4]
// double *V_ref; // 3D double array: [size1][size2][size3]
// double *sets->dV_ref; // 5D double array: [size1][size2][size3][size4][size5]
// double complex *sets->propagator; // 2D complex array: [size1][size2]
// double complex *sets->propagator_path; // 2D complex array: [size1][size2]
// double complex *sets->propagator_switch; // 2D complex array: [size1][size2]

// double *sets->R_nuc; // 2D double array: [size1][size2]
// double *sets->P_nuc; // 2D double array: [size1][size2]
// double *sets->mass; // 2D double array: [size1][size2]

// double *sets->force; // 2D double array: [size1][size2]
// double *sets->force_nuc; // 2D double array: [size1][size2]
// double *sets->force_ref; // 3D double array: [size1][size2][size3]
// double *sets->force_nuc_ref; // 3D double array: [size1][size2][size3]
// // int tysets->pe_traj_sed;

// // double seth->beta, seth->temperature;

// // char seth->method[20];

// // double dt, seth->ttot,
// double sets->t_now;
// double *sets->timegrid; // 1D double array: [size1]
// // long long Nbreak, seth->Ngrid;
// // long long Ntraj;

// // char seth->unit_t[20];
// // double seth->unittrans_t;

// // int seth->outputtype;

// // int seth->calforcetype;

// // int ifoutputmpi;

// // int seth->sampletype;

// // int seth->if_st_nan;
// unsigned long long *sets->N_nan_sum; // 1D int array: [size1]

// // int if_traj;

// // int tysets->pe_phase;

// // int tysets->pe_ad_fssh;

// // int if_ref;

// // int if_1st;

// // int if_inv_focus;

// // int seth->if_Pdis, seth->s_N;
// // double s_start, s_end;
// double *s; // 1D double array: [size1]
// double *real_sets->expisP; // 1D double array: [size1]
// double *imag_sets->expisP; // 1D double array: [size1]
// double complex *sets->expisP; // 1D complex array: [size1]
// double complex *mpi_sets->expisP; // 1D complex array: [size1]
// double *mpi_real_sets->expisP;
// double *mpi_imag_sets->expisP;

// double *sets->R_nuc_ref; // 3D double array: [size1][size2][size3]
// double *sets->P_nuc_ref; // 3D double array: [size1][size2][size3]
// double *sets->xe_ref; // 2D double array: [size1][size2]
// double *sets->pe_ref; // 2D double array: [size1][size2]
// double *ye_ref; // 2D double array: [size1][size2]
// double *psets->xe_ref; // 2D double array: [size1][size2]
// double *pye_ref; // 2D double array: [size1][size2]
// double complex *sets->gamma_cv_ref; // 3D complex array: [size1][size2][size3]
// double complex *sets->den_e_ref; // 3D complex array: [size1][size2][size3]
// double complex *inverse_kernel; // 2D complex array: [size1][size2]

// // int Nref;

// double *pop0; // 1D double array: [size1]

// // bool seth->if_ad_nac;

// // int seth->mean_nuc;

// // bool seth->if_occ_rand;

// // int if_engconsv;
// double complex *engconsv_adjmat; // 2D complex array: [size1][size2]

// // int if_RBC;

// double sets->measure_mash;
// // int tysets->pe_mash;
// double complex sets->rho0_mash[4], rhot_mash[4];

// double sets->U0_mash[4], mea_mat_mash[4];

// double complex *rho0_unsmash; // 2D complex array: [size1][size2]
// double complex *rhot_unsmash; // 2D complex array: [size1][size2]
// double *sets->U0_unsmash; // 2D double array: [size1][size2]
// double *mea_mat_unsmash; // 2D double array: [size1][size2]

// // int ifBA;
// double *sets->E_adia_old_traj; // 1D double array: [size1]
// double *sets->nac_BAeff; // 4D double array: [size1][size2][size3][size4]
// double *tdc_BA; // 2D double array: [size1][size2]
// double *sets->nac_old; // 4D double array: [size1][size2][size3][size4]
// double *sets->P_nuc_BA_old; // 2D double array: [size1][size2]

// int iflbd, occ_lbd;
// double *rate_lbd; // 1D double array: [size1]
// double rate_para_lbd;

// int ifmfsh, occ_mfsh;
// double thres;

// int occ_rdmfocus;

// int ifrw;
// double rwfactor, seth->beta_rw, seth->temperature_rw;

// // int ifmodprop;
// double complex *sets->propagator_ref; // 2D complex array: [size1][size2]

// // int ifcorresets->den;

// int memorylength;
// int itime_save, i_re_save;

// // int if_st_eng;
// double *sets->energy_est; // 1D double array: [size1]
// double *mpi_sets->energy_est; // 1D double array: [size1]

// // int seth->type_evo_ele;

// double *A_jump; // 2D double array: [size1][size2]
// double *lambda_jump; // 1D double array: [size1]
// double *U_jump; // 2D double array: [size1][size2]
// double *alpha_jump; // 1D double array: [size1]
// int tysets->pe_jump;

// int if_switchcv, occ_switchcv, max_switchcv;

// double Pn0_mf2, Wn0_mf2, energy_mf2;
// int if_inv_evo, if_BO_mf2, id_max_mf2sh, tysets->pe_eom_mf2sh, id_sets->init_occ_mf2, if_bak_mf2cv;
// double *da_mf2; // 1D double array: [size1]
// double *sets->R_nuc_old_mf2cv; // 2D double array: [size1][size2]
// double *sets->P_nuc_old_mf2cv; // 2D double array: [size1][size2]
// double *sets->xe_old_mf2cv; // 1D double array: [size1]
// double *sets->pe_old_mf2cv; // 1D double array: [size1]
// double *dE_old_mf2cv; // 1D double array: [size1]
// double *dE_mf2cv; // 1D double array: [size1]

// // int seth->type_prop_adia;

// double complex *G_xpconfg; // 2D complex array: [size1][size2]
// double *sets->permutation_ms2; // 3D double array: [size1][size2][size3]
// double complex *G_ms2; // 2D complex array: [size1][size2]
// double thres_ms2;
// int *index_ms2; // 1D int array: [size1]
// int *index_old_ms2; // 1D int array: [size1]

// double *sets->R_nuc_state; // 2D double array: [size1][size2]
// double *sets->P_nuc_state; // 2D double array: [size1][size2]
// double *V_state; // 2D double array: [size1][size2]
// double *dv_state; // 4D double array: [size1][size2][size3][size4]
// double *sets->dv_adia_state; // 4D double array: [size1][size2][size3][size4]
// double *sets->U_d2a_state; // 2D double array: [size1][size2]
// double *sets->U_ref_state; // 2D double array: [size1][size2]
// double *sets->force_nuc_state; // 2D double array: [size1][size2]
// int if_statetraj;

// bool ifBC_BCMF;

// int scheme_cvsh;

// // int seth->ifswitchforce, ifmashsets->force;

// // int seth->ifscaleenergy;

// // double seth->gamma_rescale;
// // int seth->ifscalegamma;

// double E_conserve;

// int ifmsbranch, tysets->pe_traj_msbranch;
// int itime_start_msbranch, itime_end_msbranch;
// int i_re_start_msbranch, igrid_start_msbranch;
// int iwrong_msbranch;
// double time_start_msbranch, time_end_msbranch;
// double *sets->R_nuc_brapoint; // 2D double array: [size1][size2]
// double *sets->P_nuc_brapoint; // 2D double array: [size1][size2]
// double *sets->xe_brapoint; // 1D double array: [size1]
// double *sets->pe_brapoint; // 1D double array: [size1]
// double complex *sets->gamma_cv_brapoint; // 2D complex array: [size1][size2]
// double complex *sets->den_e_brapoint; // 2D complex array: [size1][size2]

// double complex *sets->correfun_t_oldtraj; // 3D complex array: [size1][size2][size3]
// double *sets->R_nuc_oldtraj; // 3D double array: [size1][size2][size3]
// double *sets->P_nuc_oldtraj; // 3D double array: [size1][size2][size3]

// // int direc_padj;

// // int ifcount;

// // int ifreflp;

// // int ifreflp_mash;

// // int ifhardwall;

// double *eig_cv; // 1D double array: [size1]
// double *eig_cv_mat; // 2D double array: [size1][size2]
// double complex *commu_vari; // 2D complex array: [size1][size2]

// // int ifzsets->pecorr;

// // int iflangevin;
// // double eta_langevin;

// double scale_sqc2;

// // int seth->type_algorithm;

// // int tysets->pe_prop_4cont;

// // int seth->scaleenergy_type;

// // int n_step_algo5;

// // int allow_hop;

// // int seth->if_traceless_force;






void initial_vari(struct set_slave *sets,struct set_host *seth) {
    int i;


    // printf("11111111111\n");
    // printf("N=%d %d %d\n",seth->Ndof1,seth->Ndof2,seth->Nstate);

    // 分配内存
    
   
    sets->R_nuc = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->P_nuc = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->mass = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->force = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    // printf("22222\n");

    if (seth->forcetype == 1) {
        sets->force_nuc = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    }

    
    sets->V = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    sets->dV = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));

    sets->R_nuc_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->P_nuc_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->force_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->V_old = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    sets->dV_old = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    
    sets->propagator = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
    // sets->propagator_switch = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
    // if (seth->ifmodprop == 1) {
    //     sets->propagator_ref = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
    // }

    // seth->Ngrid = (long long)(seth->ttot / seth->dt) / seth->Nbreak + 1;
    //  printf("33333\n");
     
    sets->timegrid = (double *)malloc(seth->Ngrid * sizeof(double));
    
    //  printf("345454545\n");
    if (seth->if_allcf == 0) {
        if (seth->outputtype >= 0) {
            sets->den = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double complex));
            memset(sets->den, 0, seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double complex));
           
        }
        if (seth->outputtype != 0) {
            sets->population = (double *)malloc(seth->Nstate * seth->Ngrid * sizeof(double));
            memset(sets->population, 0, seth->Nstate * seth->Ngrid * sizeof(double));
            
        }
    // } else if (seth->if_allcf == 1) {
    //     sets->cf0 = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
    //     memset(sets->cf0, 0, seth->Nstate * seth->Nstate * sizeof(double complex ));
    //     // cfall = (double *)malloc(seth->Nstate * seth->Nstate * seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double ));
    //     // memset(cfall, 0, seth->Nstate * seth->Nstate * seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
    //     seth->type_evo = 1;
    // } else if (seth->if_allcf >= 2) {
    //     sets->cf0 = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
    //     memset(sets->cf0, 0, seth->Nstate * seth->Nstate * sizeof(double complex ));
    //     sets->cfeff = (double  complex *)malloc(seth->Ngrid * sizeof(double complex ));
    //     memset(sets->cfeff, 0, seth->Ngrid * sizeof(double complex ));
    //     seth->type_evo = 1;
    //     sets->weight0 = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->weightt = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    // }
    
    // if (seth->if_st_eng == 1) {
    //     sets->energy_est = (double *)malloc(seth->Ngrid * sizeof(double));
    // }
    //   printf("44444\n");
    
    // if (iflbd > 0) {
    //     rate_lbd = (double *)malloc(seth->Nstate * sizeof(double));
    //     memset(rate_lbd, 0, seth->Nstate * sizeof(double));
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    //     occ_lbd = sets->init_occ;
    }

    // if (ifmfsh > 0) {
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    // }

    sets->correfun_t = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
    memset(sets->correfun_t, 0, seth->Nstate * seth->Nstate * sizeof(double complex));
    sets->gamma_cv = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
    memset(sets->gamma_cv, 0, seth->Nstate * seth->Nstate * sizeof(double complex));
    sets->gamma_cv_old = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
    
    // if (seth->ifcv == 3) {
    //     commu_vari = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     eig_cv = (double *)malloc(seth->Nstate * sizeof(double));
    //     eig_cv_mat = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    // }

    // if (strcmp(seth->method, "GDTWA") == 0 || strcmp(seth->method, "eGDTWA") == 0) {
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    //     if (seth->if_inv_focus == 1) {
    //         sets->inverse_kernel = (double  complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double  complex ));
    //     }
   
    // } else if (strcmp(seth->method, "test") == 0) {
    //     // seth->type_evo = 1;
    // } else if (strcmp(seth->method, "DISH") == 0 || strcmp(seth->method, "dish") == 0) {
    //     t_decoh = (double *)malloc(seth->Nstate * sizeof(double));
    //     t_coh = (double *)malloc(seth->Nstate * sizeof(double));
    //     L_red = (double *)malloc(seth->Nstate * sizeof(double));
    // } else if (strcmp(seth->method, "A-FSSH") == 0 || strcmp(seth->method, "a-fssh") == 0 ||
    //            strcmp(seth->method, "AFSSH") == 0 || strcmp(seth->method, "afssh") == 0) {
    //     deltaR_afssh = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     deltaP_afssh = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     deltaF_afssh = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    // } else if (strcmp(seth->method, "GFSH") == 0 || strcmp(seth->method, "gfsh") == 0 ||
    //            strcmp(seth->method, "SC-FSSH") == 0 || strcmp(seth->method, "sc-fssh") == 0 ||
    //            strcmp(seth->method, "CC-FSSH") == 0 || strcmp(seth->method, "cc-fssh") == 0) {
    //     sets->xe_old = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->pe_old = (double *)malloc(seth->Nstate * sizeof(double));
    // } else if (strcmp(seth->method, "rdmfocus") == 0) {
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    // } else if (strcmp(seth->method, "focus_jump") == 0) {
    //     A_jump = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     lambda_jump = (double *)malloc(seth->Nstate * sizeof(double));
    //     U_jump = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     alpha_jump = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->propagator_path = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    // } else if (strcmp(seth->method, "mf2-sh") == 0 || strcmp(seth->method, "MF2-SH") == 0 ||
    //            strcmp(seth->method, "mf2cv") == 0 || strcmp(seth->method, "MF2CV") == 0) {
    //     da_mf2 = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->R_nuc_old_mf2cv = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->P_nuc_old_mf2cv = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->xe_old_mf2cv = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->pe_old_mf2cv = (double *)malloc(seth->Nstate * sizeof(double));
    //     dE_old_mf2cv = (double *)malloc(seth->Nstate * sizeof(double));
    //     dE_mf2cv = (double *)malloc(seth->Nstate * sizeof(double));
    // } else if (strcmp(seth->method, "ms2") == 0 || strcmp(seth->method, "MS2") == 0 ||
    //            strcmp(seth->method, "ms3") == 0 || strcmp(seth->method, "MS3") == 0) {
    //     seth->type_evo = 4;
    //     sets->correfun_0_ms2 = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->permutation_ms2 = (double *)malloc(seth->Nstate * seth->Nstate * seth->Nstate * sizeof(double));
    //     for (i = 0; i < seth->Nstate; i++) {
    //         gen_sets->permutation(seth->Nstate, &sets->permutation_ms2[i * seth->Nstate * seth->Nstate], i);
    //     }
    //             index_ms2 = (double *)malloc(seth->Nstate * sizeof(double));
    //     index_old_ms2 = (double *)malloc(seth->Nstate * sizeof(double));
    // } else if (strcmp(seth->method, "mstraj") == 0 || strcmp(seth->method, "mstraj2") == 0) {
    //     sets->R_nuc_state = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->P_nuc_state = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     dv_state = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     V_state = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->dv_adia_state = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->U_d2a_state = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->U_ref_state = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->force_nuc_state = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    // } else if (strcmp(seth->method, "bcmf") == 0 || strcmp(seth->method, "BCMF") == 0) {
    //     sets->E_adia_old_traj = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->xe_old = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->pe_old = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->R_nuc_old_traj = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->P_nuc_old_traj = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    // } else if (strcmp(seth->method, "sed2") == 0 || strcmp(seth->method, "SED2") == 0 ||
    //            strcmp(seth->method, "sed3") == 0 || strcmp(seth->method, "SED3") == 0) {
    //     if (seth->outputtype >= 0) {
    //         sets->den_traj = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
    //         memset(sets->den_traj, 0, seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
    //         sets->den2_traj = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
    //         memset(sets->den2_traj, 0, seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
    //     }
    //     if (seth->outputtype != 0) {
    //         sets->population_traj = (double *)malloc(seth->Nstate * seth->Ngrid * sizeof(double));
    //         memset(sets->population_traj, 0, seth->Nstate * seth->Ngrid * sizeof(double));
    //         sets->population2_traj = (double *)malloc(seth->Nstate * seth->Ngrid * sizeof(double));
    //         memset(sets->population2_traj, 0, seth->Nstate * seth->Ngrid * sizeof(double));
    //     }
    //     if (seth->if_st_fb == 1) {
    //         sets->pop_fb_traj = (double *)malloc(seth->Nstate * seth->Ngrid * 2 * sizeof(double));
    //         memset(sets->pop_fb_traj, 0, seth->Nstate * seth->Ngrid * 2 * sizeof(double));
    //         sets->pop_fb2_traj = (double *)malloc(seth->Nstate * seth->Ngrid * 2 * sizeof(double));
    //         memset(sets->pop_fb2_traj, 0, seth->Nstate * seth->Ngrid * 2 * sizeof(double));
    //     }
    // } else if (strcmp(seth->method, "unsmash") == 0 || strcmp(seth->method, "UNSMASH") == 0 ||
    //            strcmp(seth->method, "unSMASH") == 0 || strcmp(seth->method, "unsmash-mf") == 0 ||
    //            strcmp(seth->method, "UNSMASH-MF") == 0 || strcmp(seth->method, "unSMASH-MF") == 0) {
    //     if_tysets->pemb = 1;
    //     seth->type_evo = 5;
    //     sets->xe_mb = (double *)malloc(2 * (seth->Nstate - 1) * sizeof(double));
    //     sets->pe_mb = (double *)malloc(2 * (seth->Nstate - 1) * sizeof(double));
    //     state_unsmash = (double *)malloc((seth->Nstate - 1) * sizeof(double));
    //     rho0_unsmash = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     rhot_unsmash = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->U0_unsmash = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     gp_unsmash = (double *)malloc(seth->Nstate * sizeof(double));
    //     gc_unsmash = (double *)malloc(seth->Nstate * sizeof(double));
    //     mea_mat_unsmash = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    // }


    // if (ifmsbranch > 0) {
    //     sets->R_nuc_state = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->P_nuc_state = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     dv_state = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     V_state = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->dv_adia_state = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->U_d2a_state = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->U_ref_state = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->force_nuc_state = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->R_nuc_brapoint = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->P_nuc_brapoint = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->xe_brapoint = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->pe_brapoint = (double *)malloc(seth->Nstate * sizeof(double));
    //     sets->gamma_cv_brapoint = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->correfun_t_oldtraj = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
    // }
    // if (ifcount == 1) {
    //     count_st = (int *)malloc(5 * seth->Ngrid * sizeof(int));
    //     count_sets->pertraj = (int *)malloc(5 * sizeof(int));
    //     memset(count_st, 0, 5 * seth->Ngrid * sizeof(int));
    //     memset(count_sets->pertraj, 0, 5 * sizeof(int));
    // }
    sets->init_occ = seth->init_occ_4read;
    if (sets->init_occ == 0) {
        seth->if_occ_rand = 1;
    }
    // if (if_ref == 1) {
    //     sets->R_nuc_ref = (double *)malloc(Nref * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->P_nuc_ref = (double *)malloc(Nref * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->dV_ref = (double *)malloc(Nref * seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     V_ref = (double *)malloc(Nref * seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->force_ref = (double *)malloc(Nref * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->force_nuc_ref = (double *)malloc(Nref * seth->Ndof1 * seth->Ndof2 * sizeof(double));
    //     sets->xe_ref = (double *)malloc(Nref * seth->Nstate * sizeof(double));
    //     sets->pe_ref = (double *)malloc(Nref * seth->Nstate * sizeof(double));
    //     sets->den_e_ref = (double *)malloc(Nref * seth->Nstate * seth->Nstate * sizeof(double));
    //     sets->gamma_cv_ref = (double *)malloc(Nref * seth->Nstate * seth->Nstate * sizeof(double));
    //     prop_pldm = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    //     memset(prop_pldm, 0, seth->Nstate * seth->Nstate * sizeof(double));
    //     for (i = 0; i < seth->Nstate; i++) {
    //         prop_pldm[i * seth->Nstate + i] = 1;
    //     }
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    // }
    switch (seth->type_evo) {
        case 0:
        case 5:
            sets->xe = (double *)malloc(seth->Nstate * sizeof(double));
            sets->pe = (double *)malloc(seth->Nstate * sizeof(double));
            sets->xe_old = (double *)malloc(seth->Nstate * sizeof(double));
            sets->pe_old = (double *)malloc(seth->Nstate * sizeof(double));
            break;
        case 1:
        case 2:
            sets->xe = (double *)malloc(seth->Nstate * sizeof(double));
            sets->pe = (double *)malloc(seth->Nstate * sizeof(double));
            sets->den_e = (double  complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double  complex ));
            sets->xe_old = (double *)malloc(seth->Nstate * sizeof(double));
            sets->pe_old = (double *)malloc(seth->Nstate * sizeof(double));
            sets->den_e_old = (double  complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double  complex ));
            break;
        // case 3:
        //     sets->xe = (double *)malloc(seth->Nstate * sizeof(double));
        //     sets->pe = (double *)malloc(seth->Nstate * sizeof(double));
        //     sets->den_e4nuc = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        //     sets->den_e = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        //     break;
        // case 4:
        //     sets->xe = (double *)malloc(seth->Nstate * sizeof(double));
        //     sets->pe = (double *)malloc(seth->Nstate * sizeof(double));
        //     G_xpconfg = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        //     break;
    }
    
    // sets->xe_cv = (double *)malloc(seth->Nstate * sizeof(double));
    // sets->pe_cv = (double *)malloc(seth->Nstate * sizeof(double));
    if (seth->sampletype == 2) {
        seth->rep = 1;
    }
    if (seth->rep == 1 || seth->sampletype == 3) {
        sets->U_d2a = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        sets->U_ref = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        sets->U_d2a_old = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        sets->U_ref_old = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        memset(sets->U_ref, 0, seth->Nstate * seth->Nstate * sizeof(double));
        sets->E_adia = (double *)malloc(seth->Nstate * sizeof(double));
        sets->dv_adia = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        sets->nac = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        sets->E_adia_old = (double *)malloc(seth->Nstate * sizeof(double));
        sets->dv_adia_old = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        sets->nac_old = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        // P_kin = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
        sets->nac_check = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        sets->nac_check_old = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        if (seth->type_prop_adia > 0) {
            sets->overlap_adia = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        }
        memset(sets->nac_check, 0, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        // if (ifBA == 1) {
        //     sets->E_adia_old_traj = (double *)malloc(seth->Nstate * sizeof(double));
        //     sets->nac_BAeff = (double *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));
        //     tdc_BA = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        //     sets->P_nuc_BA_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
        //     memset(tdc_BA, 0, seth->Nstate * seth->Nstate * sizeof(double));
        // }
    }
    
    if (seth->if_st_fb == 1) {
        sets->pop_fb = (double *)malloc(seth->Nstate * seth->Ngrid * 2 * sizeof(double));
        memset(sets->pop_fb, 0, seth->Nstate * seth->Ngrid * 2 * sizeof(double));
    }
    // if (seth->if_Pdis == 1) {
    //     sets->s = (double *)malloc(seth->s_N * sizeof(double));
    //     sets->expisP = (double  complex *)malloc(seth->s_N * sizeof(double complex ));
    //     memset(sets->expisP, 0, seth->s_N * sizeof(double complex ));
    //     for (i = 0; i < seth->s_N; i++) {
    //         sets->s[i] = seth->s_start + i * (seth->s_end - seth->s_start) / seth->s_N;
    //     }
    // }
   
    // if (seth->temperature != 0.0 && seth->beta == 0.0) {
    //     seth->beta = 1 / (kb * seth->temperature);
    //     if (seth->temperature < 1.0) seth->beta = 10000000;
    // }
    
    // if (ifrw > 0) {
    //     if (seth->beta_rw < 0) {
    //         seth->beta_rw = 1.0 / (kb * seth->temperature_rw);
    //         if (seth->temperature_rw < 1.0) seth->beta_rw = 10000000;
    //     }
    // }
    // if (strcmp(seth->unit_t, "au") == 0) {
    //     seth->unittrans_t = 1.0;
    // } else if (strcmp(seth->unit_t, "fs") == 0) {
    //     seth->unittrans_t = 1.0 / au_2_fs;
    // }
    
    // seth->ttot *= seth->unittrans_t;
    // seth->dt *= seth->unittrans_t;

    // if (seth->mean_nuc == 1) {
    //     sets->R_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     sets->P_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     sets->R2_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     sets->P2_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     memset(sets->R_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     memset(sets->P_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     memset(sets->R2_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     memset(sets->P2_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     // if (ifmsbranch > 0) {
    //     //     sets->R_nuc_oldtraj = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     //     sets->P_nuc_oldtraj = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
    //     // }
    // }

    // printf("111166\n");
    // if (if_engconsv == 1) {
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    //     engconsv_adjmat = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    // }
    // printf("111166\n");
    //  printf("seth->Ngrid=%d\n",seth->Ngrid);
    // unsigned long long Ntemp=seth->Ngrid;
    sets->N_nan_sum = (unsigned long long *)malloc(seth->Ngrid * sizeof(unsigned long long));
    // printf("22266\n");
    // mpi_N_nan_sum = (unsigned long long *)malloc(seth->Ngrid * sizeof(unsigned long long));
    // printf("333366\n");
    memset(sets->N_nan_sum, 0, seth->Ngrid * sizeof(unsigned long long));
    // printf("444466\n");
    if (seth->ifswitchforce > 0) {
        if (seth->rep == 0) {
            sets->U_d2a = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
            sets->E_adia = (double *)malloc(seth->Nstate * sizeof(double));
            sets->U_d2a_old = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
            sets->E_adia_old = (double *)malloc(seth->Nstate * sizeof(double));
        }
    }

    // printf("666666\n");

}



void sample_ele(struct set_slave *sets,struct set_host *seth) {
    double theta[seth->Nstate], action[seth->Nstate];
    double thetaref[seth->Nstate], actionref[seth->Nstate];
    double *R0 = NULL, *Es = NULL, *C = NULL;
    double complex CC[seth->Nstate * seth->Nstate];
    double E[seth->Nstate], E_diag[seth->Nstate * seth->Nstate];
    int i, j, k, l, nid, iref;
    double p1, x1, x2, ps1, ps2;
    double alpha_mash, beta_mash;
    double c_main[seth->Nstate], sumc_main[seth->Nstate];
    double matA[seth->Nstate * seth->Nstate], matB[seth->Nstate * seth->Nstate];
    double xe_save[seth->Nstate], pe_save[seth->Nstate];
    double complex gamma_cv_save[seth->Nstate * seth->Nstate], den_e_save[seth->Nstate * seth->Nstate];
    double tempdm1[seth->Nstate * seth->Nstate], tempdm2[seth->Nstate * seth->Nstate], tempdv1[seth->Nstate], tempdv2[seth->Nstate];
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate * seth->Nstate];

    // Generate random numbers for theta
    for (i = 0; i < seth->Nstate; i++) {
        theta[i] = ((double) rand() / RAND_MAX) * 2 * M_PI;
    }

    if (seth->if_occ_rand) {
        x2 = (double) rand() / RAND_MAX;
        sets->init_occ = 1 + floor(x2 * seth->Nstate);
    }

    if (strcmp(seth->method, "MFT") == 0 || strcmp(seth->method, "mft") == 0 ||
        strcmp(seth->method, "BCMF") == 0 || strcmp(seth->method, "bcmf") == 0) {
        
        if (seth->if_occ_rand) {
            x2 = (double) rand() / RAND_MAX;
            sets->init_occ = 1 + floor(x2 * seth->Nstate);
        }

        for (i = 0; i < seth->Nstate; i++) {
            action[i] = 0;
        }
        action[sets->init_occ-1] = 1;
    
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
            
        }
        
        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        sets->correfun_0 = 0.5 * (sets->xe[sets->init_occ-1] * sets->xe[sets->init_occ-1] + sets->pe[sets->init_occ-1] * sets->pe[sets->init_occ-1]);

        if (seth->type_evo >= 1) {
            for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->den_e[i] = 0.0 + 0.0 * I;
            }
            sets->den_e[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] = 1.0 + 0.0 * I;
        }

        if (seth->if_allcf != 0) {
            for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->cf0[i] = sets->den_e[i] - sets->gamma_cv[i];
            }
        }

        // if (if_ref == 1) {
        //     for (iref = 0; iref < Nref; iref++) {
        //         for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
        //             sets->gamma_cv_ref[iref * seth->Nstate * seth->Nstate + i] = 0.0 + 0.0 * I;
        //         }

        //         if (seth->type_evo >= 1) {
        //             for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
        //                 sets->den_e_ref[iref * seth->Nstate * seth->Nstate + i] = 0.0 + 0.0 * I;
        //             }
        //             sets->den_e_ref[iref * seth->Nstate * seth->Nstate + sets->init_occ * seth->Nstate + sets->init_occ] = 1.0 + 0.0 * I;
        //         }
        //     }
        // }

    } else if (strcmp(seth->method, "eCMM") == 0 || strcmp(seth->method, "ecmm") == 0 ||
        strcmp(seth->method, "CMM") == 0 || strcmp(seth->method, "cmm") == 0) {
        
        // 调用 random_prob 函数
        random_prob(seth->Nstate, action);

        // 计算 sets->xe 和 sets->pe
        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * (1 + seth->Nstate * seth->gamma_zpe) * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * (1 + seth->Nstate * seth->gamma_zpe) * action[i]) * sin(theta[i]);
        }

        // 初始化 sets->gamma_cv
        for (int i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
            }
        }

        // 计算 sets->correfun_0
        if (seth->index_t0 == 0) {
            sets->correfun_0 = seth->Nstate * (0.5 * (sets->xe[sets->init_occ-1] * sets->xe[sets->init_occ-1] + sets->pe[sets->init_occ-1] * sets->pe[sets->init_occ-1]) - seth->gamma_zpe);
        } else if (seth->index_t0 == 1) {
            sets->correfun_0 = seth->Nstate * 0.5 * (sets->xe[seth->index_t0_1-1] - I * sets->pe[seth->index_t0_1-1]) * (sets->xe[seth->index_t0_2-1] + I * sets->pe[seth->index_t0_2-1]);
        }

        // 计算 sets->den_e
        if (seth->type_evo >= 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

        // // 计算 sets->cf0
        // if (seth->if_allcf != 0) {
        //     for (int i = 0; i < seth->Nstate; i++) {
        //         for (int j = 0; j < seth->Nstate; j++) {
        //             sets->cf0[i * seth->Nstate + j] = seth->Nstate * (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
        //         }
        //     }
        // }
    } else if (strcmp(seth->method, "focus") == 0 || strcmp(seth->method, "langer") == 0 ||
        strcmp(seth->method, "deltamft") == 0 || strcmp(seth->method, "dmft") == 0) {
        
        // 调用 random_prob 函数
        // random_prob(seth->Nstate, action);
        for (i = 0; i < seth->Nstate; i++) {
            action[i] = seth->gamma_zpe;
        }
        action[sets->init_occ - 1] += 1.0;

        // 计算 sets->xe 和 sets->pe
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }

        // 初始化 sets->gamma_cv
        for (i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
            }
        }

        // 计算 sets->correfun_0
        if (seth->index_t0 == 0) {
            sets->correfun_0 = 1;
        } else if (seth->index_t0 == 1) {
            sets->correfun_0 = seth->Nstate * 0.5 * (sets->xe[seth->index_t0_1-1] - I * sets->pe[seth->index_t0_1-1]) * (sets->xe[seth->index_t0_2-1] + I * sets->pe[seth->index_t0_2-1]);
        }

        // 计算 sets->den_e
        if (seth->type_evo >= 1) {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

        // // 计算 sets->cf0
        // if (seth->if_allcf != 0) {
        //     for (int i = 0; i < seth->Nstate; i++) {
        //         for (int j = 0; j < seth->Nstate; j++) {
        //             sets->cf0[i * seth->Nstate + j] = seth->Nstate * (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
        //         }
        //     }
        // }
    } else if (strcmp(seth->method, "sqc") == 0 || strcmp(seth->method, "SQC") == 0) {
        for (i = 0; i < seth->Nstate; i++) {
            action[i] = ((double) rand() / RAND_MAX);
        }
        p1 = 100000;
        while (1.0 - action[sets->init_occ - 1] < p1){
            for (i = 0; i < seth->Nstate; i++) {
                action[i] = ((double) rand() / RAND_MAX);
            }
            p1 = ((double) rand() / RAND_MAX);
        }
        for (i = 0; i < seth->Nstate; i++) {
            if (i == sets->init_occ - 1) continue;
            action[i] *= (1.0 - action[sets->init_occ - 1]);
        }
        action[sets->init_occ - 1] += 1;
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }

        seth->gamma_zpe = 1.0 / 3.0 ;
        // 初始化 sets->gamma_cv
        for (i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
            }
        }

        // 计算 sets->correfun_0
        // if (seth->index_t0 == 0) {
            sets->correfun_0 = 1.0;
        // } else if (seth->index_t0 == 1) {
        //     sets->correfun_0 = seth->Nstate * 0.5 * (sets->xe[seth->index_t0_1-1] - I * sets->pe[seth->index_t0_1-1]) * (sets->xe[seth->index_t0_2-1] + I * sets->pe[seth->index_t0_2-1]);
        // }

        // 计算 sets->den_e
        if (seth->type_evo >= 1) {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }
    }

    if (seth->ifcv == -1 || seth->ifcv == 1) {
        for (int i = 0; i < seth->Nstate; i++) {
            if (i == sets->init_occ - 1) {
                sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - 1;
                if (seth->ifscalegamma == 1) {
                    sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - 1;
                }
            } else {
                sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
                if (seth->ifscalegamma == 1) {
                    sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe);
                }
            }
        }    
    }

    if (seth->ifcv == -2 || seth->ifcv == 2) {
        for (int i = 0; i < seth->Nstate; i++) {
            if (i == sets->init_occ - 1) {
                sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - 1;
                if (seth->ifscalegamma == 1) {
                    sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - 1;
                }
            } else {
                sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
                if (seth->ifscalegamma == 1) {
                    sets->gamma_cv[i * seth->Nstate + i] = 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe);
                }
            }
        }
        for (int i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                if (i == j || i == sets->init_occ || j == sets->init_occ) continue;
                sets->gamma_cv[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                if (seth->ifscalegamma == 1) {
                    sets->gamma_cv[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe);
                }
            }
        }
    }



    seth->if_ad_nac = 0;
    if (seth->sampletype == 2){
        
        dV_msmodel(sets->R_nuc, sets->dV,seth);
        V_msmodel(sets->R_nuc, sets->V, 0.0,seth);
        cal_NACV(sets,seth);

        // sets->xe = matmul(transpose(sets->U_d2a), sets->xe)
        transpose(sets->U_d2a, tempdm1, seth->Nstate);
        memcpy(tempdv1, sets->xe, seth->Nstate * sizeof(double));
        dd_matmul(tempdm1, tempdv1, sets->xe, seth->Nstate, seth->Nstate, 1);
        memcpy(tempdv1, sets->pe, seth->Nstate * sizeof(double));
        dd_matmul(tempdm1, tempdv1, sets->pe, seth->Nstate, seth->Nstate, 1);
        // sets->pe = matmul(transpose(sets->U_d2a), sets->pe)
        dc_matmul(tempdm1,sets->gamma_cv,tempcm1, seth->Nstate, seth->Nstate, seth->Nstate);
        cd_matmul(tempcm1, sets->U_d2a, sets->gamma_cv, seth->Nstate, seth->Nstate, seth->Nstate);
        // sets->gamma_cv = matmul(matmul(transpose(sets->U_d2a), sets->gamma_cv), sets->U_d2a)
        
        // 这里需要实现矩阵乘法和转置操作

        if (seth->type_evo == 1 || seth->type_evo == 3) {
            // sets->den_e = matmul(matmul(transpose(sets->U_d2a), sets->den_e), sets->U_d2a)
            dc_matmul(tempdm1,sets->den_e,tempcm1, seth->Nstate, seth->Nstate, seth->Nstate);
            cd_matmul(tempcm1, sets->U_d2a, sets->den_e, seth->Nstate, seth->Nstate, seth->Nstate);
        }
        // if (sets->den_e4nuc != NULL) {
        //     // sets->den_e4nuc = matmul(matmul(transpose(sets->U_d2a), sets->den_e4nuc), sets->U_d2a)
        // }
        // if (seth->type_evo == 4) {
        //     // G_xpconfg = matmul(transpose(sets->U_d2a), G_xpconfg)
        // }

    //     if (strcmp(seth->method, "fssh") == 0 || strcmp(seth->method, "FSSH") == 0 || strcmp(seth->method, "SC-FSSH") == 0 || strcmp(seth->method, "sc-fssh") == 0 ||
    //         strcmp(seth->method, "CC-FSSH") == 0 || strcmp(seth->method, "cc-fssh") == 0 || strcmp(seth->method, "fsshswitch") == 0 || strcmp(seth->method, "pcsh") == 0 ||
    //         strcmp(seth->method, "PCSH") == 0 || strcmp(seth->method, "PCSH-NAF") == 0 || strcmp(seth->method, "pcsh-naf") == 0 || strcmp(seth->method, "BCSH") == 0 ||
    //         strcmp(seth->method, "bcsh") == 0 || strcmp(seth->method, "BCSH-NAF") == 0 || strcmp(seth->method, "bcsh-naf") == 0) {

    //         double x2 = (double)rand() / RAND_MAX;
    //         double ps1, ps2;
    //         for (int i = 0; i < seth->Nstate; i++) {
    //             sets->id_state = i;
                
    //             if (i == 0) {
    //                 ps1 = 0;
    //                 ps2 = sets->U_d2a[sets->init_occ * seth->Nstate + 0] * sets->U_d2a[sets->init_occ * seth->Nstate + 0];
    //             } else {
    //                 ps1 = 0;
    //                 for (int k = 0; k < i - 1; k++) {
    //                     ps1 += sets->U_d2a[sets->init_occ * seth->Nstate + k] * sets->U_d2a[sets->init_occ * seth->Nstate + k];
    //                 }
    //                 ps2 = ps1 + sets->U_d2a[sets->init_occ * seth->Nstate + i] * sets->U_d2a[sets->init_occ * seth->Nstate + i];
    //             }
    //             if (x2 >= ps1 && x2 < ps2) {
    //                 break;
    //             }
    //         }
    //     } else if (strcmp(seth->method, "mash") == 0 || strcmp(seth->method, "MASH") == 0) {
            
    //         if (0.5 * (sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0]) >= 0.5) {
    //             sets->id_state = 0;
    //         } else {
    //             sets->id_state = 1;
    //         }
    //         sets->measure_mash = fabs(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0] - sets->xe[1] * sets->xe[1] - sets->pe[1] * sets->pe[1]);
           
    //         memcpy(sets->U0_mash, sets->U_d2a, seth->Nstate * seth-> Nstate * sizeof(double));
    //         memset(sets->rho0_mash, 0, seth->Nstate * seth-> Nstate * sizeof(double));
    //         if (0.5 * (sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0]) >= 0.5) {
    //             sets->rho0_mash[0] = 1;
    //         } else {
    //             sets->rho0_mash[3] = 1;
    //         }
    //         sets->rho0_mash[1] = 0.5 * (sets->xe[0] + I * sets->pe[0]) * (sets->xe[1] - I * sets->pe[1]);
    //         sets->rho0_mash[2] = 0.5 * (sets->xe[0] - I * sets->pe[0]) * (sets->xe[1] + I * sets->pe[1]);
    //     // } else if (strcmp(seth->method, "unsmash") == 0 || strcmp(seth->method, "UNSMASH") == 0 || strcmp(seth->method, "unSMASH") == 0) {
    //     //     memcpy(sets->U0_unsmash, sets->U_d2a, sizeof(sets->U0_unsmash));
    //     } else if (strcmp(seth->method, "mash-mf") == 0 || strcmp(seth->method, "MASH-MF") == 0) {
            
    //         if (0.5 * (sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0]) >= 0.5) {
    //             sets->id_state = 0;
    //         } else {
    //             sets->id_state = 1;
    //         }
    //     } else if (strcmp(seth->method, "ms-mash") == 0 || strcmp(seth->method, "MS-MASH") == 0 || strcmp(seth->method, "ms-mash-focus") == 0 || strcmp(seth->method, "MS-MASH-focus") == 0 ||
    //                strcmp(seth->method, "mash-mf3") == 0 || strcmp(seth->method, "MASH-MF3") == 0) {
           
    //         // double max_val = sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0];
    //         // for (int i = 1; i < seth->Nstate; i++) {
    //         //     double val = sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i];
    //         //     if (val > max_val) {
    //         //         max_val = val;
    //         //         sets->id_state = i;
    //         //     }
    //         // }
    //         for (int i = 0; i < seth->Nstate; i++){
    //             tempdv1[i] = sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i];
    //         }
    //         sets->id_state = maxloc(tempdv1, seth->Nstate);

    //     // } else if (strcmp(seth->method, "mf2-sh") == 0 || strcmp(seth->method, "MF2-SH") == 0) {
    //     //     int sets->id_state = 0;
    //     //     double max_val = sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0];
    //     //     for (int i = 1; i < seth->Nstate; i++) {
    //     //         double val = sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i];
    //     //         if (val > max_val) {
    //     //             max_val = val;
    //     //             sets->id_state = i;
    //     //         }
    //     //     }
    //     // } else if (strcmp(seth->method, "cvsh") == 0 || strcmp(seth->method, "CVSH") == 0) {
    //     //     double action[seth->Nstate];
    //     //     for (int i = 0; i < seth->Nstate; i++) {
    //     //         if (seth->type_evo == 1) {
    //     //             action[i] = sets->den_e[i * seth->Nstate + i];
    //     //         } else {
    //     //             action[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) / 2 - sets->gamma_cv[i * seth->Nstate + i];
    //     //         }
    //     //     }
    //     //     int sets->id_state = 0;
    //     //     double max_val = action[0];
    //     //     for (int i = 1; i < seth->Nstate; i++) {
    //     //         if (action[i] > max_val) {
    //     //             max_val = action[i];
    //     //             sets->id_state = i;
    //     //         }
    //     //     }
    //     }
    }

    // if(seth->ifswitchforce>0 && seth->ifswitchforce<4){
    //     if(seth->rep == 0){
    //         V_msmodel(sets->R_nuc, sets->V, 0.0,seth);
    //         dia_symmat(seth->Nstate,sets->V,sets->E_adia,sets->U_d2a);
    //         if(seth->type_evo == 0){
    //             memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
    //             memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
    //             transpose(U_d2a,tempdm1,seth->Nstate);
    //             dd_matmul(tempdm1,xe_save,xe,seth->Nstate,seth->Nstate,1);
    //             dd_matmul(tempdm1,pe_save,pe,seth->Nstate,seth->Nstate,1);
    //         } else if (seth->type_evo == 1) {
    //             memcpy(den_e_save,den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
    //             transpose(U_d2a,tempdm1,seth->Nstate);

    //         }
            
    //     }
    // }
    if (seth->ifswitchforce > 0) {
        if (seth->rep == 0) {
            
            V_msmodel(sets->R_nuc, sets->V, 0.0, seth);
            #ifdef sunway
            int slavecore_id = athread_get_id(-1);
            #endif
          
            // printf("%d %18.8E %18.8E %18.8E %18.8E\n",slavecore_id,sets->V[0],sets->V[1],sets->V[2],sets->V[3]);
                
           
            
            dia_symmat(seth->Nstate, sets->V, sets->E_adia, sets->U_d2a);

            // printf("%d\n",slavecore_id);

            if (seth->type_evo == 0) {
                memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                transpose(sets->U_d2a,tempdm1,seth->Nstate);
                dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
            } else if (seth->type_evo == 1) {
                memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                transpose(sets->U_d2a,tempdm1,seth->Nstate);
                dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
                cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            }
            memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
            dc_matmul(tempdm1,gamma_cv_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
            cd_matmul(tempcm1,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        }

        for (int i = 0; i < seth->Nstate; i++) {
            if (seth->type_evo == 0) {
                c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) / 2 - creal(sets->gamma_cv[i * seth->Nstate + i]);
            } else if (seth->type_evo == 1) {
                c_main[i] = creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]);
            }
        }

        if (seth->ifscalegamma == 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                if (seth->type_evo == 0) {
                    c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) / 2 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]);
                } else if (seth->type_evo == 1) {
                    c_main[i] = creal(sets->den_e[i * seth->Nstate + i] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + i]);
                }
            }
        }

        sets->id_state = maxloc(c_main, seth->Nstate);

        if (seth->rep == 0) {
            if (seth->type_evo == 0) {
                memcpy(sets->xe,xe_save,seth->Nstate*sizeof(double));
                memcpy(sets->pe,pe_save,seth->Nstate*sizeof(double));
            } else if (seth->type_evo == 1) {
                memcpy(sets->den_e,den_e_save,seth->Nstate * seth->Nstate * sizeof(double complex));
            }
            memcpy(sets->gamma_cv,gamma_cv_save,seth->Nstate * seth->Nstate * sizeof(double complex));
        }
    }


}



void cal_correfun(struct set_slave *sets,struct set_host *seth) {
    int i, j, k, l;
    double x1, x2;
    double ct_sqc;
    double x0[seth->Nstate], p0[seth->Nstate], pop_test[seth->Nstate], rad_test;
    double theta_cycle[2], theta_sw, varphi_sw, psi_gp;
    double basisv1[3], basisv2[3], pvec_sw[3];
    double weight_switch[2], A[2 * 2], b[2], E_vec[2], E_mat[2 * 2], U[2 * 2];
    double complex c_st[seth->Nstate], sum_cst[seth->Nstate];
    double action[seth->Nstate], theta[seth->Nstate];
    int i_st;
    double rdm_st;
    double alpha_mash, beta_mash;
    double tempv[seth->Nstate];
    double complex tempcm[seth->Nstate*seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate],tempdm2[seth->Nstate*seth->Nstate];
    double complex tempcm2[seth->Nstate*seth->Nstate];
    #ifdef sunway
    int slavecore_id = athread_get_id(-1);
    #endif
    
    
    // tempmat1[seth->Nstate*seth->Nstate],tempmat2[seth->Nstate*seth->Nstate],tempmat3[seth->Nstate*seth->Nstate]



    if (seth->sampletype == 2) {
        // sets->xe = matmul(sets->U_d2a, sets->xe)
        // matmul(sets->U_d2a, sets->xe, seth->Nstate, seth->Nstate, 1);
        
        memcpy(tempv,sets->xe,seth->Nstate * sizeof(double));
        dd_matmul(sets->U_d2a,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
        // sets->pe = matmul(sets->U_d2a, sets->pe)
        // matmul(sets->U_d2a, sets->pe, seth->Nstate, seth->Nstate, 1);
        memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
        dd_matmul(sets->U_d2a,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
        // sets->gamma_cv = matmul(matmul(sets->U_d2a, sets->gamma_cv), transpose(sets->U_d2a))
        // matmul(sets->U_d2a, sets->gamma_cv, seth->Nstate, seth->Nstate, seth->Nstate);
        
        memcpy(tempcm,sets->gamma_cv,seth->Nstate * seth->Nstate*sizeof(double complex));
        
        transpose(sets->U_d2a,tempdm,seth->Nstate);
        cd_matmul(sets->gamma_cv,tempdm,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
        dc_matmul(sets->U_d2a,tempcm2,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);

        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            memcpy(tempcm,sets->den_e,seth->Nstate*seth->Nstate*sizeof(double complex));
            cd_matmul(sets->den_e,tempdm,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            dc_matmul(sets->U_d2a,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
        }
        // if (seth->type_evo == 4) {
        //     // G_xpconfg = matmul(sets->U_d2a, G_xpconfg)
        //     matmul(sets->U_d2a, G_xpconfg, seth->Nstate, seth->Nstate, 1);
        // }
    }

    
    if (strcmp(seth->method, "MFT") == 0 || strcmp(seth->method, "mft") == 0 ||
        strcmp(seth->method, "BCMF") == 0 || strcmp(seth->method, "bcmf") == 0) {
        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->correfun_t[i] = sets->den_e[i];
            }
        } else {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->correfun_t[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

    } else if (strcmp(seth->method, "eCMM") == 0 || strcmp(seth->method, "ecmm") == 0 ||
        strcmp(seth->method, "CMM") == 0 || strcmp(seth->method, "cmm") == 0) {
        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->correfun_t[i] = sets->den_e[i] * (1.0 + seth->Nstate) / pow(1 + seth->Nstate * seth->gamma_zpe, 2);
            }
            for (int i = 0; i < seth->Nstate; i++) {
                sets->correfun_t[i * seth->Nstate + i] -= (1.0 - seth->gamma_zpe) / (1 + seth->Nstate * seth->gamma_zpe);
            }
        } else {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->correfun_t[i * seth->Nstate + j] = (1.0 + seth->Nstate) / (2 * pow(1 + seth->Nstate * seth->gamma_zpe, 2)) * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - (1.0 - seth->gamma_zpe) / (1 + seth->Nstate * seth->gamma_zpe) * (i == j ? 1 : 0);
                }
            }
        }

    } else if (strcmp(seth->method, "focus") == 0 || strcmp(seth->method, "langer") == 0 ||
        strcmp(seth->method, "deltamft") == 0 || strcmp(seth->method, "dmft") == 0) {
        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->correfun_t[i] = sets->den_e[i] ;
            }
            for (i = 0; i < seth->Nstate; i++) {
                sets->correfun_t[i * seth->Nstate + i] -= seth->gamma_zpe;
            }
        } else {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->correfun_t[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - seth->gamma_zpe * (i == j ? 1 : 0);
                }
            }
        }

    } else if (strcmp(seth->method, "sqc") == 0 || strcmp(seth->method, "SQC") == 0) {
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (i = 0; i < seth->Nstate ; i++) {
                sets->correfun_t[i * seth->Nstate + i] = 1 ;
                for (j = 0; j < seth->Nstate; j++){
                    if( (j == i && creal(sets->den_e[j * seth->Nstate + j]) < 1) || 
                        (j != i && creal(sets->den_e[j * seth->Nstate + j]) >= 1)){
                        sets->correfun_t[i * seth->Nstate + i] = 0 ;
                    }
                    if(i != j) sets->correfun_t[i*seth->Nstate+j] = sets->den_e[i*seth->Nstate+j];
                }
            }
        } else {
            for (i = 0; i < seth->Nstate ; i++) {
                sets->correfun_t[i * seth->Nstate + i] = 1 ;
                for (j = 0; j < seth->Nstate; j++){
                    if( (j == i && 0.5 * (sets->xe[j] * sets->xe[j] + sets->pe[j] * sets->pe[j]) < 1) || 
                        (j != i && 0.5 * (sets->xe[j] * sets->xe[j] + sets->pe[j] * sets->pe[j]) >= 1)){
                        sets->correfun_t[i * seth->Nstate + i] = 0 ;
                    }
                    if(i != j) sets->correfun_t[i*seth->Nstate+j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }
    }





    if (seth->sampletype == 2) {
        // sets->xe = matmul(sets->U_d2a, sets->xe)
        // matmul(sets->U_d2a, sets->xe, seth->Nstate, seth->Nstate, 1);
        transpose(sets->U_d2a,tempdm,seth->Nstate);
        memcpy(tempv,sets->xe,seth->Nstate*sizeof(double));
        dd_matmul(tempdm,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
        // sets->pe = matmul(sets->U_d2a, sets->pe)
        // matmul(sets->U_d2a, sets->pe, seth->Nstate, seth->Nstate, 1);
        memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
        dd_matmul(tempdm,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
        // sets->gamma_cv = matmul(matmul(sets->U_d2a, sets->gamma_cv), transpose(sets->U_d2a))
        // matmul(sets->U_d2a, sets->gamma_cv, seth->Nstate, seth->Nstate, seth->Nstate);
       
        memcpy(tempcm,sets->gamma_cv,seth->Nstate*seth->Nstate*sizeof(double complex));
        cd_matmul(sets->gamma_cv,sets->U_d2a,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
        dc_matmul(tempdm,tempcm2,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);

        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            memcpy(tempcm,sets->den_e,seth->Nstate*seth->Nstate*sizeof(double complex));
            cd_matmul(sets->den_e,sets->U_d2a,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            dc_matmul(tempdm,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
        }
        // if (seth->type_evo == 4) {
        //     // G_xpconfg = matmul(sets->U_d2a, G_xpconfg)
        //     matmul(sets->U_d2a, G_xpconfg, seth->Nstate, seth->Nstate, 1);
        // }
    }


}


void cal_propagator(int Nstate, double *H, double dt, double complex *U,struct set_slave *sets,struct set_host *seth) {
    int i, j;
    double E[Nstate], C[Nstate * Nstate];
    double sineig[Nstate * Nstate], coseig[Nstate * Nstate];
    double real_pro[Nstate * Nstate], img_pro[Nstate * Nstate];
    
    // printf("%18.8E %18.8E %18.8E %18.8E \n",H[0],H[1],H[2],H[3]);
    // 调用dia_symmat函数
    dia_symmat(Nstate, H, E, C);
    //  printf("z  1111111 \n");
    // 初始化sineig和coseig
    // for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
    //     sineig[i] = 0.0;
    //     coseig[i] = 0.0;
    // }
    // printf("223333222\n");
    memset(sineig,0,Nstate*Nstate*sizeof(double));
    memset(coseig,0,Nstate*Nstate*sizeof(double));
    // printf("2222222222222\n");
    // 计算sineig和coseig
    for (i = 0; i < Nstate; i++) {
        sineig[i * Nstate + i] = -1.0*sin(E[i] * dt);
        coseig[i * Nstate + i] = cos(E[i] * dt);
    }
    //  printf("z 22222222  \n");
    // 计算real_pro和img_pro
    double tempdm1[Nstate*Nstate],tempdm2[Nstate*Nstate];
    
    transpose(C, tempdm1, Nstate);
    dd_matmul(coseig,tempdm1,tempdm2,Nstate,Nstate,Nstate);
    dd_matmul(C,tempdm2,real_pro,Nstate,Nstate,Nstate);
    dd_matmul(sineig,tempdm1,tempdm2,Nstate,Nstate,Nstate);
    dd_matmul(C,tempdm2,img_pro,Nstate,Nstate,Nstate);
    
    // printf("z  3333333 \n");
    
    // 计算U
    for (i = 0; i < Nstate * Nstate; i++) {
        U[i] = real_pro[i] + I * img_pro[i];
    }
    //  printf("z  4444444444 \n");
    // 检查sets->E_adia是否已分配
    // printf("%d %d\n",seth->ifswitchforce, seth->rep);
    if (seth->ifswitchforce > 0 && seth->rep == 0) {
        // printf("11111\n");
        memcpy(sets->U_d2a,C,Nstate*Nstate*sizeof(double));
        memcpy(sets->E_adia,E,Nstate*sizeof(double));
        // for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
        //     sets->U_d2a[i] = C[i];
        //     if(i<seth->Nstate) {sets->E_adia[i] = E[i]};
        // }
        // for (i = 0; i < seth->Nstate; i++) {
        //     sets->E_adia[i] = E[i];
        // }
    }
    //  printf("z 55555  \n");
}
        
    
void evo_traj_ele(double deltat,struct set_slave *sets,struct set_host *seth) {
    double x0[seth->Nstate], p0[seth->Nstate];
    double x0_mb[2 * (seth->Nstate - 1)], p0_mb[2 * (seth->Nstate - 1)];
    double complex propagator_unsmash[2 * 2 * (seth->Nstate - 1)];
    int i;
    double complex tempv1[seth->Nstate],tempv2[seth->Nstate];
    double complex tempcm1[seth->Nstate*seth->Nstate],tempcm2[seth->Nstate*seth->Nstate];
    // if (sets->U_d2a_old != NULL) memcpy(sets->U_d2a_old, sets->U_d2a, sizeof(sets->U_d2a));
    // printf("11111\n");
    switch (seth->type_evo) {
        case 0:
            // printf("y1111111111\n");
            if (seth->rep == 0) cal_propagator(seth->Nstate, sets->V, deltat, sets->propagator,sets,seth);
            //  printf("y2222222222\n");
            if (seth->rep == 1) cal_propagator_adia(seth->Nstate, deltat, sets->propagator,sets,seth);
            memcpy(x0, sets->xe, seth->Nstate * sizeof(double));
            memcpy(p0, sets->pe, seth->Nstate * sizeof(double));
            //  printf("y3333333333333\n");
            // matmul_real_imag(seth->Nstate, sets->propagator, x0, p0, sets->xe, sets->pe);
            // printf("22222\n");
            
            cd_matmul(sets->propagator,x0,tempv1,seth->Nstate,seth->Nstate,1);
            cd_matmul(sets->propagator,p0,tempv2,seth->Nstate,seth->Nstate,1);
            //  printf("y4444444444444444444\n");
            for(i=0;i<seth->Nstate;i++){
                sets->xe[i]=creal(tempv1[i])-cimag(tempv2[i]);
                sets->pe[i]=creal(tempv2[i])+cimag(tempv1[i]);
            }
            // printf("y5555555555555\n");
            break;
        
        case 1:
            if (seth->rep == 0) cal_propagator(seth->Nstate, sets->V, deltat, sets->propagator,sets,seth);
            if (seth->rep == 1) cal_propagator_adia(seth->Nstate, deltat, sets->propagator,sets,seth);
            // matmul_complex(seth->Nstate, sets->propagator, sets->den_e, sets->den_e);
            
            diagger(sets->propagator,tempcm1,seth->Nstate);
            cc_matmul(sets->den_e,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(sets->propagator,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            // if (inverse_kernel != NULL) {
            //     matmul_complex(seth->Nstate, sets->propagator, inverse_kernel, inverse_kernel);
            // }
            // break;
        
            
        // case 2:
        //     if (seth->rep == 0) cal_propagator(seth->Nstate, V, deltat, sets->propagator);
        //     if (seth->rep == 1 && seth->type_evo_ele == 0) cal_propagator_adia(seth->Nstate, deltat, sets->propagator);
        //     memcpy(x0, sets->xe, seth->Nstate * sizeof(double));
        //     memcpy(p0, sets->pe, seth->Nstate * sizeof(double));
        //     matmul_real_imag(seth->Nstate, sets->propagator, x0, p0, sets->xe, sets->pe);
        //     matmul_complex(seth->Nstate, sets->propagator, sets->den_e, sets->den_e);
        //     break;
        // case 3:
        //     if (seth->rep == 0) cal_propagator(seth->Nstate, V, deltat, sets->propagator);
        //     if (seth->rep == 1 && seth->type_evo_ele == 0) cal_propagator_adia(seth->Nstate, deltat, sets->propagator);
        //     matmul_complex(seth->Nstate, sets->propagator, sets->den_e, sets->den_e);
        //     matmul_complex(seth->Nstate, sets->propagator, sets->den_e4nuc, sets->den_e4nuc);
        //     break;
        // case 4:
        //     if (seth->rep == 0) cal_propagator(seth->Nstate, V, deltat, sets->propagator);
        //     if (seth->rep == 1 && seth->type_evo_ele == 0) cal_propagator_adia(seth->Nstate, deltat, sets->propagator);
        //     memcpy(x0, sets->xe, seth->Nstate * sizeof(double));
        //     memcpy(p0, sets->pe, seth->Nstate * sizeof(double));
        //     matmul_real_imag(seth->Nstate, sets->propagator, x0, p0, sets->xe, sets->pe);
        //     matmul_complex(seth->Nstate, sets->propagator, G_xpconfg, G_xpconfg);
        //     break;
        // case 5:
        //     if (seth->rep == 0) cal_propagator_dia_unsmash(seth->Nstate, V, sets->id_state, state_unsmash, deltat, sets->propagator_unsmash);
        //     if (seth->rep == 1 && seth->type_evo_ele == 0) cal_propagator_adia_unsmash(seth->Nstate, sets->id_state, state_unsmash, deltat, sets->propagator_unsmash);
        //     memcpy(x0_mb, sets->xe_mb, 2 * (seth->Nstate - 1) * sizeof(double));
        //     memcpy(p0_mb, sets->pe_mb, 2 * (seth->Nstate - 1) * sizeof(double));
        //     for (i = 0; i < seth->Nstate - 1; i++) {
        //         matmul_real_imag(2, &sets->propagator_unsmash[4 * i], &x0_mb[2 * i], &p0_mb[2 * i], &sets->xe_mb[2 * i], &sets->pe_mb[2 * i]);
        //     }
        //     break;
    }

    if (seth->ifcv == 1 || seth->ifcv == 2) {
        // matmul_complex(seth->Nstate, sets->propagator, sets->gamma_cv, sets->gamma_cv);
        diagger(sets->propagator,tempcm1,seth->Nstate);
        cc_matmul(sets->gamma_cv,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
        cc_matmul(sets->propagator,tempcm2,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
    }
}



// void evo_traj_ele_prop(double deltat){

// TODO !!!!!!!!!!!!!!!!!!!!!!!!!!

// }

void evo_traj_nucP(double deltat,struct set_slave *sets,struct set_host *seth) {
    // 假设 sets->P_nuc 和 sets->force 是一维数组
   
    for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
        sets->P_nuc[i] += sets->force[i] * deltat;
    }

    // if (seth->ifscaleenergy == 3) energy_conserve_naf_3(deltat);
}

void evo_traj_nucR(double deltat,struct set_slave *sets,struct set_host *seth) {
    double x1, x2;

    // if (iflangevin == 1) {
        // for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
        //     sets->R_nuc[i] += sets->P_nuc[i] * 0.5 * deltat / sets->mass[i];
        // }

        // for (int i = 0; i < seth->Ndof1; i++) {
        //     for (int j = 0; j < seth->Ndof2; j++) {
        //         box_muller(&x2, &x1, 1.0, 0.0);
        //         int index = i * seth->Ndof2 + j;
        //         sets->P_nuc[index] = exp(-eta_langevin * deltat) * sets->P_nuc[index] + sqrt((1 - exp(-2 * eta_langevin * deltat)) * sets->mass[index] / seth->beta) * x2;
        //     }
        // }

        // for (int i = 0; i < Ndof; i++) {
        //     sets->R_nuc[i] += sets->P_nuc[i] * 0.5 * deltat / sets->mass[i];
        // }
    // } else {
        // if (seth->ifswitchforce == 3 && sum_abs_sets->dv_adia(sets->id_state, sets->id_state) > 1e-30) {
        //     x2 = 0.0;
        //     for (int i = 0; i < seth->Nstate; i++) {
        //         for (int j = 0; j < seth->Nstate; j++) {
        //             if (i == j) continue;
        //             x2 += (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - sets->gamma_cv[i * seth->Nstate + j]) * (sets->E_adia[j] - sets->E_adia[i]) * sum_sets->nac_sets->P_nuc_sets->mass(i, j);
        //         }
        //     }
        //     for (int i = 0; i < Ndof; i++) {
        //         sets->R_nuc[i] += sets->P_nuc[i] * deltat / sets->mass[i] + x2 / sets->dv_adia[sets->id_state * Ndof + sets->id_state] * deltat;
        //     }
        // } else {
            for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
                sets->R_nuc[i] += sets->P_nuc[i] * deltat / sets->mass[i];
            }
        // }
    // }
}


void evo_traj_calProp(int igrid_cal,struct set_slave *sets,struct set_host *seth) {
    int i, j, icfall;
    double x2;
   
    cal_correfun(sets,seth);

    // if (ifmsbranch > 0) {
    //     if (if_statetraj == -1) {
    //         for (i = 0; i < seth->Nstate; i++) {
    //             for (j = 0; j < seth->Nstate; j++) {
    //                 sets->correfun_t[i * seth->Nstate + j] -= sets->correfun_t_oldtraj[i * seth->Nstate + j + igrid_cal * seth->Nstate * seth->Nstate];
    //                 sets->correfun_t_oldtraj[i * seth->Nstate + j + igrid_cal * seth->Nstate * seth->Nstate] = sets->correfun_t[i * seth->Nstate + j] + sets->correfun_t_oldtraj[i * seth->Nstate + j + igrid_cal * seth->Nstate * seth->Nstate];
    //             }
    //         }
    //     } else {
    //         memcpy(&sets->correfun_t_oldtraj[igrid_cal * seth->Nstate * seth->Nstate], sets->correfun_t, seth->Nstate * seth->Nstate * sizeof(double));
    //     }
    // }

    if (any_isnan(sets->correfun_t, seth->Nstate * seth->Nstate)) {
        sets->N_nan_sum[igrid_cal]++;
        if (seth->if_st_nan == 1) {
            memset(sets->correfun_t, 0, seth->Nstate * seth->Nstate * sizeof(double complex));
        }
    } else {
        sets->N_nan_sum[igrid_cal] = sets->N_nan_sum[igrid_cal];
    }

    // printf("ii=%d\n",igrid_cal);

    if (seth->outputtype >= 0) {
        // if (strcmp(seth->method, "gauss") == 0 || strcmp(seth->method, "genLSC") == 0 || strcmp(seth->method, "genlsc") == 0) {
        //     if (seth->ifid == 0) {
        //         for (i = 0; i < seth->Nstate; i++) {
        //             for (j = 0; j < seth->Nstate; j++) {
        //                 sets->den[i * seth->Nstate * seth->Ngrid + j * seth->Ngrid  + igrid_cal] += sets->correfun_0 * sets->correfun_t[i * seth->Nstate + j];
        //             }
        //         }
        //     } else if (seth->ifid == 1) {
        //         for (i = 0; i < seth->Nstate; i++) {
        //             for (j = 0; j < seth->Nstate; j++) {
        //                 sets->den[i * seth->Nstate * seth->Ngrid + j * seth->Ngrid  + igrid_cal] += sets->correfun_0 * sets->correfun_t[i * seth->Nstate + j];
        //             }
        //         }
        //         for (j = 0; j < seth->Nstate; j++) {
        //             sets->den[j * seth->Nstate * seth->Ngrid + j * seth->Ngrid + igrid_cal] += 1.0 / seth->Nstate - 1.0 / (seth->sigma2_lsc * seth->sigma2_lsc * seth->Nstate * seth->Nstate);
        //         }
        //     }
        // } else {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->den[i * seth->Nstate * seth->Ngrid + j * seth->Ngrid + igrid_cal] += sets->correfun_0 * sets->correfun_t[i * seth->Nstate + j];
                }
            }
        // }
    }

   

    if (seth->outputtype != 0) {
        for (i = 0; i < seth->Nstate; i++) {
        //     if (strcmp(seth->method, "gauss") == 0 || strcmp(seth->method, "genLSC") == 0 || strcmp(seth->method, "genlsc") == 0) {
        //         if (seth->ifid == 0) {
        //             sets->population[i * seth->Ngrid  + igrid_cal] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
        //             if (seth->if_st_fb == 1) {
        //                 if (sets->P_nuc[0] > 0) {
        //                     sets->pop_fb[i * seth->Ngrid *2  + igrid_cal*2] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
        //                 } else {
        //                     sets->pop_fb[i * seth->Ngrid *2  + igrid_cal*2 + 1] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
        //                 }
        //             }
        //         } else if (seth->ifid == 1) {
        //             sets->population[i * seth->Ngrid  + igrid_cal] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]) + 1.0 / seth->Nstate - 1.0 / (seth->sigma2_lsc * seth->sigma2_lsc * seth->Nstate * seth->Nstate);
        //             if (seth->if_st_fb == 1) {
        //                 if (sets->P_nuc[0] > 0) {
        //                     sets->pop_fb[i * seth->Ngrid *2  + igrid_cal*2] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]) + 1.0 / seth->Nstate - 1.0 / (seth->sigma2_lsc * seth->sigma2_lsc * seth->Nstate * seth->Nstate);
        //                 } else {
        //                     sets->pop_fb[i * seth->Ngrid *2  + igrid_cal*2 + 1] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]) + 1.0 / seth->Nstate - 1.0 / (seth->sigma2_lsc * seth->sigma2_lsc * seth->Nstate * seth->Nstate);
        //                 }
        //             }
        //         }
        //     } else {
                sets->population[i * seth->Ngrid  + igrid_cal] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                if (seth->if_st_fb == 1) {
                    if (sets->P_nuc[0] > 0) {
                        sets->pop_fb[i * seth->Ngrid *2  + igrid_cal * 2 + 0] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                    } else {
                        sets->pop_fb[i * seth->Ngrid *2  + igrid_cal * 2 + 1] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                    }
                }
        //     }
        }
    }

    

    

    // if (allocated(cfall)) {
    //     for (i = 0; i < seth->Nstate; i++) {
    //         for (j = 0; j < seth->Nstate; j++) {
    //             for (int k = 0; k < seth->Ndof1; k++) {
    //                 for (int l = 0; l < seth->Ndof2; l++) {
    //                     cfall[i * seth->Nstate + j + k * seth->Nstate * seth->Nstate + l * seth->Nstate * seth->Nstate * seth->Ndof1 + igrid_cal * seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2] += sets->cf0[i * seth->Nstate + j] * sets->correfun_t[k * seth->Nstate + l];
    //                 }
    //             }
    //         }
    //     }
    // }

    // if (sets->cfeff!=NULL) {
    //     if (seth->if_allcf == 2) {
    //         sets->cfeff[igrid_cal] += sum(weight0, seth->Nstate) * sum(weightt, seth->Nstate) * sum(sets->correfun_t, seth->Nstate * seth->Nstate);
    //     } else if (seth->if_allcf == 3) {
    //         for (icfall = 0; icfall < allcf_times; icfall++) {
    //             cfweight_msmodel(weight0, weightt, seth->beta, icfall);
    //             sets->cfeff[igrid_cal] += sum(weight0, seth->Nstate) * sum(weightt, seth->Nstate) * sum(sets->correfun_t, seth->Nstate * seth->Nstate);
    //         }
    //     } else if (seth->if_allcf == 4) {
    //         for (icfall = 0; icfall < allcf_times; icfall++) {
    //             cfweight_msmodel(weight0, weightt, seth->beta, icfall);
    //             sets->cfeff[igrid_cal] += sum(weight0, seth->Nstate) * sum(weightt, seth->Nstate) * sum(sets->correfun_t, seth->Nstate * seth->Nstate);
    //         }
    //     }
    // }

    // if (sets->R_nuc_mean!=NULL) {
    //     if (strcmp(seth->method, "sqc") == 0 
    //     ||  strcmp(seth->method, "SQC") == 0 
    //     ||  strcmp(seth->method, "mf3") == 0 
    //     ||  strcmp(seth->method, "MF3") == 0 
    //     ||  strcmp(seth->method, "sqc2") == 0 
    //     ||  strcmp(seth->method, "SQC2") == 0 
    //     ||  strcmp(seth->method, "sqc3") == 0 
    //     ||  strcmp(seth->method, "SQC3") == 0) {
    //         x2 = 0;
    //         for (i = 0; i < seth->Nstate; i++) {
    //             x2 += sets->correfun_t[i * seth->Nstate + i];
    //         }
    //         for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
    //             sets->R_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * creal(sets->correfun_0) * x2;
    //             sets->P_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * creal(sets->correfun_0) * x2;
    //             sets->R2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0) * x2;
    //             sets->P2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0) * x2;
    //         }
    //     } else if (strcmp(seth->method, "mash") == 0 
    //     || strcmp(seth->method, "MASH") == 0) {
    //         for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
    //             sets->R_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
    //             sets->P_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
    //             sets->R2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
    //             sets->P2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
    //         }
    //     // } else if (strcmp(seth->method, "unsmash") == 0 || strcmp(seth->method, "UNSMASH") == 0 || strcmp(seth->method, "unSMASH") == 0) {
    //     //     for (i = 0; i < Ndof; i++) {
    //     //         sets->R_nuc_mean[i + igrid_cal * Ndof] += sets->R_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
    //     //         sets->P_nuc_mean[i + igrid_cal * Ndof] += sets->P_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
    //     //         sets->R2_nuc_mean[i + igrid_cal * Ndof] += sets->R_nuc[i] * sets->R_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
    //     //         sets->P2_nuc_mean[i + igrid_cal * Ndof] += sets->P_nuc[i] * sets->P_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
    //     //     }
    //     } else if (strcmp(seth->method, "mash-mf") == 0 || strcmp(seth->method, "MASH-MF") == 0) {
    //         for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
    //             sets->R_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
    //             sets->P_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
    //             sets->R2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
    //             sets->P2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
    //         }
    //     } else {
    //         for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
    //             sets->R_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * creal(sets->correfun_0);
    //             sets->P_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * creal(sets->correfun_0);
    //             sets->R2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0);
    //             sets->P2_nuc_mean[i + igrid_cal * seth->Ndof1*seth->Ndof2] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0);
    //         }
    //     }
    // }

    

    // if (ifmsbranch > 0) {
    //     if (if_statetraj == -1) {
    //         for (i = 0; i < Ndof; i++) {
    //             sets->R_nuc_mean[i + igrid_cal * Ndof] -= sets->R_nuc_oldtraj[i + igrid_cal * Ndof];
    //             sets->P_nuc_mean[i + igrid_cal * Ndof] -= sets->P_nuc_oldtraj[i + igrid_cal * Ndof];
    //         }
    //     }
    //     memcpy(&sets->R_nuc_oldtraj[igrid_cal * Ndof], sets->R_nuc, Ndof * sizeof(double));
    //     memcpy(&sets->P_nuc_oldtraj[igrid_cal * Ndof], sets->P_nuc, Ndof * sizeof(double));
    // }

    
    // if (sets->energy_est != NULL) {
    //     double x2 = 0.0;
    //     double Ekin=0.0;
    //     // V_msmodel(sets->R_nuc, V, sets->t_now); // Uncomment if needed
    //     for(int i=0;i<seth->Ndof1*seth->Ndof2;i++){
    //         Ekin+=0.5*sets->P_nuc[i]*sets->P_nuc[i]/sets->mass[i];
    //     }
    //     if (seth->rep == 1 && seth->sampletype == 0) {
    //         // cal_NACV(); // Uncomment if needed
    //         for (int i = 0; i < seth->Nstate; i++) {
    //             x2 += (sets->E_adia[i] + Ekin) * sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i];
    //         }
    //     } else {
    //         for (int i = 0; i < seth->Nstate; i++) {
    //             for (int j = 0; j < seth->Nstate; j++) {
    //                 if (i == j) {
    //                     x2 += Ekin * sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i];
    //                 }
    //                 x2 += sets->V[i * seth->Nstate + j] * sets->correfun_0 * sets->correfun_t[j * seth->Nstate + i];
    //             }
    //         }
    //     }

    //     if (seth->ifswitchforce >= 1) {
    //         // cal_NACV(); // Uncomment if needed
    //         sets->energy_est[igrid_cal] += Ekin + sets->E_adia[sets->id_state];
    //     } else {
    //         if (strcmp(seth->method, "FSSH") == 0 || strcmp(seth->method, "fssh") == 0 ||
    //             strcmp(seth->method, "mash") == 0 || strcmp(seth->method, "MASH") == 0 ||
    //             strcmp(seth->method, "ms-mash") == 0 || strcmp(seth->method, "MS-MASH") == 0 ||
    //             strcmp(seth->method, "unsmash") == 0 || strcmp(seth->method, "UNSMASH") == 0 ||
    //             strcmp(seth->method, "unSMASH") == 0) {
    //             sets->energy_est[igrid_cal] += Ekin + sets->E_adia[sets->id_state];
    //         } else {
    //             sets->energy_est[igrid_cal] += x2;
    //         }
    //     }
    // }

    

    // if (count_st != NULL) {
    //     for (int i = 0; i < seth->Nstate; i++) {
    //         count_st[i * igrid_cal] += count_sets->pertraj;
    //     }
    //     count_sets->pertraj = 0;
    

}


void energy_conserve_naf_1(double E0, double *dE_naf,struct set_slave *sets,struct set_host *seth) {
    *dE_naf = E0 - sets->E_adia[sets->id_state];
    if (*dE_naf >= 0) {
        double P_nuc_sum = 0.0;
        for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
            P_nuc_sum += sets->P_nuc[i] * sets->P_nuc[i] / sets->mass[i];
        }
        double factor = sqrt(*dE_naf / (0.5 * P_nuc_sum));
        for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
            sets->P_nuc[i] *= factor;
        }
    }
}


// void energy_conserve_naf_3(double deltat,struct set_slave *sets,struct set_host *seth) {
//     double dE_naf, x2 = 0.0;
//     for (int i = 0; i < seth->Nstate; i++) {
//         for (int j = 0; j < seth->Nstate; j++) {
//             if (i == j) continue;
//             double sum_nac = 0.0;
//             for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
//                 sum_nac += sets->nac[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k] * sets->P_nuc[k] / sets->mass[k];
//             }
//             x2 += (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i*seth->Nstate+j])) * (sets->E_adia[j] - sets->E_adia[i]) * sum_nac;
//         }
//     }
//     double P_nuc_sum = 0.0;
//     for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
//         P_nuc_sum += sets->P_nuc[i] * sets->P_nuc[i] / sets->mass[i];
//     }
//     double factor = x2 / P_nuc_sum * deltat;
//     for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
//         sets->P_nuc[i] += sets->P_nuc[i] * factor;
//     }
// }

void evo_traj_algorithm1(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    
    evo_traj_nucP(deltat / 2,sets,seth);
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    evo_traj_ele(deltat,sets,seth);
    cal_force(sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}



void evo_traj_savetraj(struct set_slave *sets,struct set_host *seth) {
    // 假设所有变量已经在def.h中声明
    memcpy(sets->R_nuc_old, sets->R_nuc, seth->Ndof1 * seth->Ndof2 * sizeof(double));
    memcpy(sets->P_nuc_old, sets->P_nuc, seth->Ndof1 * seth->Ndof2 * sizeof(double));
    memcpy(sets->xe_old, sets->xe, seth->Nstate * sizeof(double));
    memcpy(sets->pe_old, sets->pe, seth->Nstate * sizeof(double));
    // printf("x11111\n");
    if (seth->type_evo == 1) memcpy(sets->den_e_old, sets->den_e, seth->Nstate * seth->Nstate * sizeof(double complex));
    // // printf("x122222\n");
    memcpy(sets->gamma_cv_old, sets->gamma_cv, seth->Nstate * seth->Nstate * sizeof(double complex));
    memcpy(sets->V_old, sets->V, seth->Nstate * seth->Nstate * sizeof(double));
    memcpy(sets->dV_old, sets->dV, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
    // // printf("xxxxxxxxxxx\n");
    
    
    if (seth->rep == 1){
        if (sets->E_adia != NULL) memcpy(sets->E_adia_old, sets->E_adia, seth->Nstate * sizeof(double));
        if (sets->dv_adia != NULL) memcpy(sets->dv_adia_old, sets->dv_adia, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
        if (sets->nac != NULL) memcpy(sets->nac_old, sets->nac, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
        if (sets->U_d2a != NULL) memcpy(sets->U_d2a_old, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double));
        if (sets->U_ref != NULL) memcpy(sets->U_ref_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double));
        if (sets->nac_check != NULL) memcpy(sets->nac_check_old, sets->nac_check, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
    }

    if (seth->ifswitchforce > 0 && seth->rep == 0){
        if (sets->E_adia != NULL) memcpy(sets->E_adia_old, sets->E_adia, seth->Nstate * sizeof(double));
        if (sets->U_d2a != NULL) memcpy(sets->U_d2a_old, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double));
    }
    
    sets->id_state_old = sets->id_state;

    memcpy(sets->force_old, sets->force, seth->Ndof1 * seth->Ndof2 *sizeof(double));

    // if (strcmp(trim(adjustl(msmodelname)), "mole") == 0) {
    //     if (file_exists("imomap.dat")) {
    //         system("cp imomap.dat imomap_save.dat");
    //     }
    // }
}

void evo_traj_back(struct set_slave *sets,struct set_host *seth) {
    // 假设所有变量已经在def.h中声明
    memcpy(sets->R_nuc, sets->R_nuc_old, seth->Ndof1 * seth->Ndof2 * sizeof(double));
    memcpy(sets->P_nuc, sets->P_nuc_old, seth->Ndof1 * seth->Ndof2 * sizeof(double));
    memcpy(sets->xe, sets->xe_old, seth->Nstate * sizeof(double));
    memcpy(sets->pe, sets->pe_old, seth->Nstate * sizeof(double));
    if (seth->type_evo == 1) memcpy(sets->den_e, sets->den_e_old, seth->Nstate * seth->Nstate * sizeof(double complex));
    memcpy(sets->gamma_cv, sets->gamma_cv_old, seth->Nstate * seth->Nstate * sizeof(double complex));
    memcpy(sets->V, sets->V_old, seth->Nstate * seth->Nstate * sizeof(double));
    memcpy(sets->dV, sets->dV_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
    
    if (seth->rep == 1){
        memcpy(sets->E_adia, sets->E_adia_old, seth->Nstate * sizeof(double));
        memcpy(sets->dv_adia, sets->dv_adia_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
        memcpy(sets->nac, sets->nac_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
        memcpy(sets->U_d2a, sets->U_d2a_old, seth->Nstate * seth->Nstate * sizeof(double));
        memcpy(sets->U_ref, sets->U_ref_old, seth->Nstate * seth->Nstate * sizeof(double));
        memcpy(sets->nac_check, sets->nac_check_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double));
    }


    if (seth->ifswitchforce > 0 && seth->rep == 0){
        memcpy(sets->E_adia, sets->E_adia_old, seth->Nstate * sizeof(double));
        memcpy(sets->U_d2a, sets->U_d2a_old, seth->Nstate * seth->Nstate * sizeof(double));
    }
    
    sets->id_state = sets->id_state_old;


    memcpy(sets->force, sets->force_old, seth->Ndof1 * seth->Ndof2 *sizeof(double));

    // if (strcmp(trim(adjustl(msmodelname)), "mole") == 0) {
    //     if (file_exists("imomap_save.dat")) {
    //         system("cp imomap_save.dat imomap.dat");
    //     }
    // }
}


void evo_traj_new(int itraj,struct set_slave *sets,struct set_host *seth) {
    int i_re, igrid;
    int nstep, itime, icfall, iref;
    int i, j, k;//, i_lbd, hop_lbd;
    double x0[seth->Nstate], p0[seth->Nstate], bound, deltaE, deltaE_all[seth->Nstate];
    double tt1, tt2, x1, x2, sumpop, diagden[seth->Nstate];
    // double E_diag[seth->Nstate * seth->Nstate];
    // double P_s[seth->Ndof1 * seth->Ndof2], P_s_all[seth->Nstate * seth->Ndof1 * seth->Ndof2], P_s_main[seth->Ndof1 * seth->Ndof2], R_s[seth->Ndof1 * seth->Ndof2], f_s[seth->Ndof1 * seth->Ndof2], sets->pex, pt, proj[seth->Nstate];
    // double Etot, Ekin, Epot, dE;
    // double vmp[seth->Nstate * seth->Nstate];
    // double Rinit, Pinit;
    // double complex rho_save[seth->Nstate * seth->Nstate], Aforsave[seth->Nstate * seth->Nstate];
    // double sets->P_nuc_BA_old[seth->Ndof1 * seth->Ndof2];
    // double sets->P_nuc_old[seth->Ndof1 * seth->Ndof2 * memorylength], sets->R_nuc_old[seth->Ndof1 * seth->Ndof2 * memorylength], sets->E_adia_old[seth->Nstate * memorylength], sets->nac_old[seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * memorylength], sets->force_old[seth->Ndof1 * seth->Ndof2 * memorylength], deltaE_all_old[seth->Nstate * memorylength];
    // int loc_bak, id_memory, jstate;
    bool if_bak;
    // long long *approachtimes;
    // double action_switchcv[seth->Nstate], off_sets->gamma_cv, Vij;
    // double complex c_main[seth->Nstate];
    // int i_max, npack, index_state[seth->Nstate];
    // double *depack;
    // int *index_pack;
    double dt_evo;
    // int nstep_small, istep_small;
    // double sets->t_now_small;
    // bool alive;
    int slavecore_id;
    #ifdef sunway
    slavecore_id=athread_get_id(-1);
    #endif
    
    // if_bak = false;
    // itime_save = 0;

    // count_sets->pertraj = 0;

    sets->t_now = 0;
    nstep = (int)(seth->ttot / seth->dt) + 1;

    V_msmodel(sets->R_nuc, sets->V, 0.0,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    //debug
    // int slavecore_id;
    // slavecore_id=athread_get_id(-1);
    // if(slavecore_id == 0){
        // for (i=0; i<seth->Nstate; i++){
            // printf("sets->E_adia(%d)=%18.8E\n",i,sets->E_adia[i]);
            // for (j=0;j<seth->Nstate;j++){
                // printf("U(%d,%d)=%18.8E\n",i,j,sets->U_d2a[i*seth->Nstate+j]);
    //             for (k=0;k<1;k++){
    //                 for (int l=0;l<1;l++){
    //                     printf("sets->dv_adia(%d,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
    //                     printf("sets->nac(%d,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->nac[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
    //                 }
    //             }
            
    //         }
    //     }
    // }
    
   

    i_re = seth->Nbreak;
    igrid = 0;

    if (seth->ifscaleenergy > 0) {
        sets->E_conserve = 0.0;
        for(i = 0; i < seth->Ndof1*seth->Ndof2; i++){
            sets->E_conserve += sets->P_nuc[i] * sets->P_nuc[i] / sets->mass[i];
        }
        sets->E_conserve *= 0.5;
        sets->E_conserve += sets->E_adia[sets->id_state];
    }

    if (seth->ifscaleenergy > 0) {
        switch (seth->ifscaleenergy) {
            case 1:
            case 4:
            case 5:
                seth->scaleenergy_type = 1;
                break;
            case 3:
                seth->scaleenergy_type = 3;
                break;
        }
    }

    

    cal_force(sets,seth);

    // debug
    // int slavecore_id;
    // slavecore_id=athread_get_id(-1);
    // if(slavecore_id == 0){
    //     for (i=0; i<seth->Ndof1; i++){
            
    //         for (j=0;j<seth->Ndof2;j++){
    //             printf("sets->force(%d,%d)=%18.8E\n",i,j,sets->force[i*seth->Ndof1+j]);
    //         }
    //     }
    // }
    //  exit(1);


    itime = 0;


    //debug
    // int slavecore_id;
    // slavecore_id=athread_get_id(-1);
    // if(slavecore_id == 0) printf("itime=%d,sets->force=%18.8E\n",itime,sets->force[299]);

   
    
    while (itime <= nstep) {
        if (i_re >= seth->Nbreak && igrid < seth->Ngrid) {

            
            evo_traj_calProp(igrid,sets,seth);
             
            // if(slavecore_id == 10) printf("%18.8E %18.8E %18.8E\n", sets->t_now, sets->R_nuc[0], sets->P_nuc[0]);
             
            sets->timegrid[igrid] = sets->t_now;
            igrid++;
            i_re = 0;
            
            //debug
            //  if(slavecore_id == 10) printf("%18.8E %18.8E %18.8E %18.8E %18.8E \n", sets->R_nuc[0],sets->P_nuc[0],sets->xe[0],sets->pe[0],sets->force[0]);
                
            
            // if(slavecore_id == 10) printf("%18.8E %18.8E\n", sets->force[0],sets->force[299]);
            

        }

        // printf("%d %d %f %f %f %f\n",slavecore_id, itime, sets->R_nuc[0],sets->P_nuc[0],sets->xe[0],sets->pe[0]);

     

        evo_traj_savetraj(sets,seth);

        dt_evo = seth->dt;


        

        switch (seth->type_algorithm) {
            case 1:
            
                evo_traj_algorithm1(dt_evo,sets,seth);
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
        // if(slavecore_id == 0) printf("itime=%d,sets->force=%18.8E\n",itime,sets->force[299]);

  
       
        if (seth->ifscaleenergy > 0) {
            if (seth->scaleenergy_type == 1) energy_conserve_naf_1(sets->E_conserve, &deltaE, sets, seth);
        //     if (seth->ifscaleenergy == 4 && seth->scaleenergy_type == 1 && deltaE < 0) {
        //         nstep_small = 2;
        //         dt_evo /= 2;
        //         sets->t_now_small = 0;
        //         evo_traj_back();
        //         for (istep_small = 1; istep_small <= nstep_small; istep_small++) {
        //             switch (seth->type_algorithm) {
        //                 case 1:
        //                     evo_traj_algorithm1(dt_evo);
        //                     break;
        //                 case 2:
        //                     evo_traj_algorithm2(dt_evo);
        //                     break;
        //             }
        //             if (seth->scaleenergy_type == 1) energy_conserve_naf_1(E_conserve, deltaE);
        //             if (seth->scaleenergy_type == 1 && deltaE < 0) {
        //                 if (dt_evo > dt / 1024) {
        //                     nstep_small *= 2;
        //                     dt_evo /= 2;
        //                     evo_traj_back();
        //                 }
        //             }
        //         }
        //         if (seth->scaleenergy_type == 3) seth->scaleenergy_type = 1;
        //     }
        }

        sets->t_now = (itime + 1) * seth->dt;
        i_re++;
        itime++;
        // if (ifzsets->pecorr > 0) zsets->pecorr_msmodel(sets->P_nuc, sets->R_nuc, ifzsets->pecorr);

        if (strcmp(seth->msmodelname, "morse3") == 0 || strcmp(seth->msmodelname, "Morse3") == 0) {
            if (seth->ifhardwall == 1) {
                if (sets->P_nuc[0] < 0 && sets->R_nuc[0] < 0) sets->P_nuc[0] = -sets->P_nuc[0];
            }
        }

        
        
    }

    // if (seth->if_Pdis == 1) {
    //     if (strcmp(seth->method, "mash") == 0 || strcmp(seth->method, "MASH") == 0) {
    
    //         for(i=0;i<seth->s_N;i++){
    //             // sets->expisP[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * 2 * sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
    //             sets->expisP[i] += sets->correfun_0 *(cos(sets->P_nuc[0] * sets->s[i]) + I * sin(sets->P_nuc[0] * sets->s[i])) * 2 * sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
    //         }
    //             // break;
    //     } else if (strcmp(seth->method, "unsmash") == 0 ||
    //               strcmp(seth->method, "UNSMASH") == 0 ||
    //               strcmp(seth->method, "unSMASH") == 0 ){
    //         for(i=0;i<seth->s_N;i++){
    //             // sets->expisP[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * seth->Nstate * rho0_unsmash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
    //             sets->expisP[i] += sets->correfun_0 * (cos(sets->P_nuc[0] * sets->s[i]) + I * sin(sets->P_nuc[0] * sets->s[i])) * seth->Nstate * sets->rho0_unsmash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
    //         }
    //             // break;
    //     } else if (strcmp(seth->method, "mash-mf") == 0 ||
    //                 strcmp(seth->method, "MASH-MF") == 0 ){
    //         for(i=0;i<seth->s_N;i++){
    //             // sets->expisP[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * 2 * sets->measure_mash;
    //             sets->expisP[i] += sets->correfun_0 * (cos(sets->P_nuc[0] * sets->s[i]) + I* sin(sets->P_nuc[0] * sets->s[i])) * 2 * sets->measure_mash;
    //         }
    //             // break;
    //     } else if (strcmp(seth->method, "sqc") == 0 ||
    //                strcmp(seth->method, "SQC") == 0 ||                 
    //                strcmp(seth->method, "mf3") == 0 ||
    //                strcmp(seth->method, "MF3") == 0 ||
    //                strcmp(seth->method, "sqc2") == 0 ||
    //                strcmp(seth->method, "SQC2") == 0 ||
    //                strcmp(seth->method, "sqc3") == 0 ||
    //                strcmp(seth->method, "SQC3") == 0 ){
    //         x2 = 0;
    //         for (i = 0; i < seth->Nstate; i++) {
    //             x2 += creal(sets->correfun_t[i * seth->Nstate + i]);
    //         }
    //         for(i=0;i<seth->s_N;i++){
    //             // sets->expisP[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * x2;
    //             sets->expisP[i] += sets->correfun_0 * ( cos( sets->P_nuc[0] * sets->s[i]) + I * sin( sets->P_nuc[0] * sets->s[i]) ) * x2;
    //         }
    //             // break;

    //     } else {
    //             for(i=0;i<seth->s_N;i++){
    //                 // sets->expisP[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]);
    //                 sets->expisP[i] += sets->correfun_0 * ( cos( sets->P_nuc[0] * sets->s[i]) + I * sin( sets->P_nuc[0] * sets->s[i]) );
    //             }
    //             // break;
    //     }
    // }
    
}


void cal_force(struct set_slave *sets,struct set_host *seth) {
    int iref, i, j;
    double frdm, x2;

    // sets->force = 0.0;
    
    memset(sets->force,0,seth->Ndof1*seth->Ndof2*sizeof(double));

    if (strcmp(seth->method, "FSSH") == 0 || strcmp(seth->method, "fssh") == 0) {
        // cal_force_fssh();
        // sets->force = -sets->dv_adia[sets->id_state][sets->id_state];
    } else {
        // if (if_ref == 1) {
        //     cal_force_mf_ref();
        // } else if (if_1st > 0) {
        //     cal_force_mf_1st();
        // } else if (ifmsbranch > 0) {
        //     cal_force_msbranch();
        if (seth->ifswitchforce == 1) {
            cal_force_switch(sets,seth);
        // } else if (seth->ifswitchforce == 2) {
        //     cal_force_switch2();
        // } else if (seth->ifswitchforce == 3) {
        //     cal_force_switch();
        // } else if (seth->ifswitchforce == 4) {
        //     cal_force_switch4();
        // } else if (seth->ifswitchforce == 5) {
        //     cal_force_switch5();
        // } else if (seth->ifswitchforce == 6) {
        //     cal_force_switch6();
        // } else if (seth->ifswitchforce == 7) {
        //     cal_force_switch7();
        // } else if (seth->ifswitchforce == 8) {
        //     cal_force_switch8();
        // } else if (ifmashsets->force > 0) {
        //     cal_force_mashsets->force();
        } else {
            cal_force_mf(sets,seth);
        }
    }
    
    if (seth->forcetype == 1) {
        nucforce_msmodel(sets->R_nuc, sets->force_nuc, seth);
        for(i=0 ; i<seth->Ndof1*seth->Ndof2; i++){
            sets->force[i] -= sets->force_nuc[i];
        }
        // if (if_ref == 1) {
        //     for (iref = 1; iref <= Nref; iref++) {
        //         nucsets->force_msmodel(sets->R_nuc_ref[iref], sets->force_nuc_ref[iref]);
        //         sets->force_ref[iref] -= sets->force_nuc_ref[iref];
        //     }
        // }
    }
    
}



void cal_force_mf(struct set_slave *sets,struct set_host *seth) {
    int i, j, k;
    double force_trace[seth->Ndof1 * seth->Ndof2];

    // if (seth->if_traceless_force == 1) {
    //     // for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
    //     //     sets->force_trace[i] = 0;
    //     // }
    //     memset(force_trace,0,seth->Ndof1*seth->Ndof2*sizeof(double));
    //     for (i = 0; i < seth->Nstate; i++) {
    //         for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    //             force_trace[j] += sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] / seth->Nstate;
    //         }
    //     }
    //     for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    //         for (i = 0; i < seth->Nstate; i++) {
    //             sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] -= force_trace[j];
    //         }
    //         sets->force[j] = -force_trace[j];
    //     }
    //     // for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    //     //     sets->force[j] = -sets->force_trace[j];
    //     // }
    // }

     
   

    switch (seth->type_evo) {
        case 0:
        case 2:
            if (seth->rep == 0) {
                if (seth->calforcetype == 1) {
                    for (i = 0; i < seth->Nstate; i++) {
                        for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
                            sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + i]));
                        }
                    }
                } else {
                    for (i = 0; i < seth->Nstate; i++) {
                        for (j = 0; j < seth->Nstate; j++) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * ((sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + j]));
                            }
                        }
                    }
                }
            } else if (seth->rep == 1) {
                for (i = 0; i < seth->Nstate; i++) {
                    for (j = 0; j < seth->Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - creal(sets->gamma_cv[i * seth->Nstate + i])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k];
                                // printf("sets->force(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
                            }
                        } else {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                                // sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                            }
                        }
                    }
                }
            }
            break;
        case 1:
            if (seth->rep == 0) {
                if (seth->calforcetype == 1) {
                    for (i = 0; i < seth->Nstate; i++) {
                        for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
                            sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]);
                        }
                    }
                } else {
                    for (i = 0; i < seth->Nstate; i++) {
                        for (j = 0; j < seth->Nstate; j++) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * creal(sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
                            }
                        }
                    }
                }
            } else if (seth->rep == 1) {
                for (i = 0; i < seth->Nstate; i++) {
                    for (j = 0; j < seth->Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k];
                            }
                        } else {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal(sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) * (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                            }
                        }
                    }
                }
            }
            break;
        // case 3:
        //     if (seth->rep == 0) {
        //         if (seth->calforcetype == 1) {
        //             for (i = 0; i < seth->Nstate; i++) {
        //                 for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
        //                     sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * creal(sets->den_e4nuc[i * seth->Nstate + i]);
        //                 }
        //             }
        //         } else {
        //             for (i = 0; i < seth->Nstate; i++) {
        //                 for (j = 0; j < seth->Nstate; j++) {
        //                     for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
        //                         sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * creal(sets->den_e4nuc[i * seth->Nstate + j]);
        //                     }
        //                 }
        //             }
        //         }
        //     } else if (seth->rep == 1) {
        //         for (i = 0; i < seth->Nstate; i++) {
        //             for (j = 0; j < seth->Nstate; j++) {
        //                 if (i == j) {
        //                     for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
        //                         sets->force[k] -= creal(sets->den_e4nuc[i * seth->Nstate + i]) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k];
        //                     }
        //                 } else {
        //                     for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
        //                         sets->force[k] -= creal(sets->den_e4nuc[i * seth->Nstate + j]) * (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
        //                     }
        //                 }
        //             }
        //         }
        //     }
        //     break;
    }
}

void cal_force_switch(struct set_slave *sets, struct set_host *seth) {
    double c_main[seth->Nstate], sumc_main[seth->Nstate];
    double xe_save[seth->Nstate], pe_save[seth->Nstate];
    double complex gamma_cv_save[seth->Nstate * seth->Nstate], den_e_save[seth->Nstate * seth->Nstate];
    double tempdm1[seth->Nstate * seth->Nstate], tempdm2[seth->Nstate * seth->Nstate], tempdv1[seth->Nstate], tempdv2[seth->Nstate];
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate * seth->Nstate];
    double deltavector[seth->Ndof1 * seth->Ndof2],P_para[seth->Ndof1 * seth->Ndof2],P_ver[seth->Ndof1 * seth->Ndof2];
    double sum;
    double Q_dia[seth->Nstate * seth->Nstate];
    double deltaE_mash;

    if (seth->rep == 0) {
        if (seth->type_evo == 0) {
            memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
            memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
            transpose(sets->U_d2a,tempdm1,seth->Nstate);
            dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
            dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
        } else if (seth->type_evo == 1) {
            memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
            transpose(sets->U_d2a,tempdm1,seth->Nstate);
            dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
            cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
        }
        memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
        dc_matmul(tempdm1,gamma_cv_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        cd_matmul(tempcm1,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        
    }

    for (int i = 0; i < seth->Nstate; i++) {
        if (seth->type_evo == 0) {
            c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) / 2 - creal(sets->gamma_cv[i * seth->Nstate + i]);
        } else if (seth->type_evo == 1) {
            c_main[i] = creal(sets->den_e[i * seth->Nstate + i]) - creal(sets->gamma_cv[i * seth->Nstate + i]);
        }
    }

    if (seth->ifscalegamma == 1) {
        for (int i = 0; i < seth->Nstate; i++) {
            if (seth->type_evo == 0) {
                c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) / 2 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]);
            } else if (seth->type_evo == 1) {
                c_main[i] = creal(sets->den_e[i * seth->Nstate + i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]);
            }
        }
    }

    int id_switch = maxloc(c_main, seth->Nstate);

    if (sets->id_state != id_switch) {
        switch (seth->direc_padj) {
            case 0:
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sets->P_nuc[i] /= sqrt(sets->mass[i]);
                    deltavector[i] = 0;
                }
                for (int i = 0; i < seth->Nstate; i++) {
                    for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
                        deltavector[j] += 1.0 / sqrt(sets->mass[j]) 
                        * creal(0.5 * (sets->xe[i] - I * sets->pe[i]) * (sets->xe[sets->id_state] + I * sets->pe[sets->id_state]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + j] 
                        - 0.5 * (sets->xe[i] - I * sets->pe[i]) * (sets->xe[id_switch] + I * sets->pe[id_switch]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + id_switch * seth->Ndof1 * seth->Ndof2 + j]);
                    }
                }
                sum = 0.0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sum += deltavector[i] * deltavector[i];
                }
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltavector[i] /= sqrt(sum);
                }
                sum = 0.0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sum += deltavector[i] * sets->P_nuc[i];
                }
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    P_para[i] = deltavector[i] * sum;
                    P_ver[i] = sets->P_nuc[i] - P_para[i];
                }
                deltaE_mash = 0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltaE_mash += P_para[i] * P_para[i];
                }
                deltaE_mash += 2 * (sets->E_adia[sets->id_state] - sets->E_adia[id_switch]);
                if (deltaE_mash >= 0) {
                    sum=0;
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sum += P_para[i] * P_para[i];
                    }
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        P_para[i] *= sqrt(deltaE_mash / sum);
                        sets->P_nuc[i] = P_para[i] + P_ver[i];
                        sets->P_nuc[i] *= sqrt(sets->mass[i]);
                    }
                    sets->id_state = id_switch;
                    // if (seth->count_pertraj) seth->count_pertraj[1] = 1;
                } else {
                    if (seth->ifreflp == 0) {
                        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                            P_para[i] = -P_para[i];
                        }
                    }
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sets->P_nuc[i] =P_para[i] + P_ver[i];
                        sets->P_nuc[i] *= sqrt(sets->mass[i]);
                    }
                    // if (seth->count_pertraj) seth->count_pertraj[2] = 1;
                    // seth->type_traj_sed = 1;
                }
                break;
            case 1:
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltavector[i] = sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + id_switch * seth->Ndof1 * seth->Ndof2 + i];
                }
                sum = 0.0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sum += deltavector[i] * deltavector[i];
                }
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltavector[i] /= sqrt(sum);
                }
                sum = 0.0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sum += deltavector[i] * sets->P_nuc[i];
                }
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    P_para[i] = deltavector[i] * sum;
                    P_ver[i] = sets->P_nuc[i] - P_para[i];
                }
                deltaE_mash = 0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltaE_mash += 0.5 * P_para[i] * P_para[i] / sets->mass[i];
                }
                deltaE_mash += (sets->E_adia[sets->id_state] - sets->E_adia[id_switch]);
                if (deltaE_mash >= 0) {
                    sum=0;
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sum += 0.5 * P_para[i] * P_para[i] / sets->mass[i];
                    }
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        P_para[i] *= sqrt(deltaE_mash / sum);
                        sets->P_nuc[i] = P_para[i] + P_ver[i];
                    }
                    sets->id_state = id_switch;
                    // if (seth->count_pertraj) seth->count_pertraj[1] = 1;
                } else {
                    if (seth->ifreflp == 0) {
                        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                            P_para[i] = -P_para[i];
                        }
                    }
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sets->P_nuc[i] =P_para[i] + P_ver[i];
                    }
                    // if (seth->count_pertraj) seth->count_pertraj[2] = 1;
                    // seth->type_traj_sed = 1;
                }
                break;

            case 2:
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltavector[i] = sets->P_nuc[i];
                }
                sum = 0.0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sum += deltavector[i] * deltavector[i];
                }
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltavector[i] /= sqrt(sum);
                }
                sum = 0.0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    sum += deltavector[i] * sets->P_nuc[i];
                }
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    P_para[i] = deltavector[i] * sum;
                    P_ver[i] = sets->P_nuc[i] - P_para[i];
                }
                deltaE_mash = 0;
                for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                    deltaE_mash += 0.5 * P_para[i] * P_para[i] / sets->mass[i];
                }
                deltaE_mash += (sets->E_adia[sets->id_state] - sets->E_adia[id_switch]);
                if (deltaE_mash >= 0) {
                    sum=0;
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sum += 0.5 * P_para[i] * P_para[i] / sets->mass[i];
                    }
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        P_para[i] *= sqrt(deltaE_mash / sum);
                        sets->P_nuc[i] = P_para[i] + P_ver[i];
                    }
                    sets->id_state = id_switch;
                    // if (seth->count_pertraj) seth->count_pertraj[1] = 1;
                } else {
                    if (seth->ifreflp == 0) {
                        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                            P_para[i] = -P_para[i];
                        }
                    }
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sets->P_nuc[i] =P_para[i] + P_ver[i];
                    }
                    // if (seth->count_pertraj) seth->count_pertraj[2] = 1;
                    // seth->type_traj_sed = 1;
                }
                break;
        }
    }

    for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
        Q_dia[i] = 0;
    }
    for (int i = 0; i < seth->Nstate; i++) {
        if (seth->type_evo == 0) {
            if (seth->ifscalegamma == 0) {
                Q_dia[i * seth->Nstate + i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + i]);
            } else {
                Q_dia[i * seth->Nstate + i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]);
            }
        } else if (seth->type_evo == 1) {
            if (seth->ifscalegamma == 0) {
                Q_dia[i * seth->Nstate + i] = creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]);
            } else {
                Q_dia[i * seth->Nstate + i] = creal(sets->den_e[i * seth->Nstate + i] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + i]);
            }
        }
    }
    for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
        Q_dia[i] = -Q_dia[i];
    }
    Q_dia[sets->id_state * seth->Nstate + sets->id_state] += 1;

    // matmul(sets->U_d2a, Q_dia, seth->Nstate, Q_dia);
    // transpose(sets->U_d2a, seth->Nstate, sets->U_d2a);
    // matmul(Q_dia, sets->U_d2a, seth->Nstate, Q_dia);
    transpose(sets->U_d2a,tempdm1,seth->Nstate);
    dd_matmul(sets->U_d2a,Q_dia,tempdm2,seth->Nstate,seth->Nstate,seth->Nstate);
    dd_matmul(tempdm2,tempdm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);

    if (seth->rep == 0) {
        if (seth->type_evo == 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                sets->xe[i] = xe_save[i];
                sets->pe[i] = pe_save[i];
            }
        } else if (seth->type_evo == 1) {
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->den_e[i] = den_e_save[i];
            }
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
            sets->gamma_cv[i] = gamma_cv_save[i];
        }
    }

    if (seth->rep == 0) {
        if (seth->calforcetype == 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    if (seth->type_evo == 0) {
                        if (seth->ifscalegamma == 0) {
                            sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]);
                        } else {
                            sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]);
                        }
                    } else if (seth->type_evo == 1) {
                        if (seth->ifscalegamma == 0) {
                            sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        (creal(sets->den_e[i * seth->Nstate + i]) - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]);
                        } else {
                            sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        (creal(sets->den_e[i * seth->Nstate + i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]);
                        }
                    }
                }
                
            }
        } else {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++){
                        if (seth->type_evo == 0) {
                            if (seth->ifscalegamma == 0) {
                                sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            ((sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + j]) + Q_dia[i * seth->Nstate + j]);
                            } else {
                                sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            ((sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) * 0.5 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + j]) + Q_dia[i * seth->Nstate + j]);
                            }
                        } else if (seth->type_evo == 1) {
                            if (seth->ifscalegamma == 0) {
                                sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            (creal(sets->den_e[i * seth->Nstate + j]) - creal(sets->gamma_cv[i * seth->Nstate + j]) + Q_dia[i * seth->Nstate + j]);
                            } else {
                                sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            (creal(sets->den_e[i * seth->Nstate + j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + j]) + Q_dia[i * seth->Nstate + j]);
                            }
                        }
                    }
                }
            }
        }
    } else if (seth->rep == 1) {
        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
            sets->force[i] = -sets->dv_adia[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + i];
        }
        for (int i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                if (i == j) continue;
                if (seth->type_evo == 0) {
                    if (seth->ifscalegamma == 0) {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) *
                                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                        }
                    } else {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + j])) *
                                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                        }
                    }
                } else if (seth->type_evo == 1) {
                    if (seth->ifscalegamma == 0) {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= (creal(sets->den_e[i * seth->Nstate + j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) *
                                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                        }
                    } else {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= (creal(sets->den_e[i * seth->Nstate + j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + j])) *
                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                        }
                    }
                }
            }
        }
    }

   


}

// void cal_force_fssh(){}


// void cal_force_eld(){}


void cal_NACV(struct set_slave *sets,struct set_host *seth){
    int i, j, imax;
    double dotp, norm_nac1, norm_nac2, cosphase;
    double overlap[seth->Nstate * seth->Nstate], overlap2[seth->Nstate * seth->Nstate], vmax;
    double gij[seth->Ndof1 * seth->Ndof2], alpha_BA;
    int id_max[seth->Nstate], idloc, id1, id2;
    double complex eiet[seth->Nstate * seth->Nstate];

    double tempdm1[seth->Nstate*seth->Nstate],tempdm2[seth->Nstate*seth->Nstate],tempdm3[seth->Nstate*seth->Nstate],tempdm4[seth->Nstate*seth->Nstate], tempdv1[seth->Nstate];

    dia_symmat(seth->Nstate, sets->V, sets->E_adia, sets->U_d2a);

    


    if (seth->if_ad_nac) {
        transpose(sets->U_d2a,tempdm1,seth->Nstate);
        dd_matmul(tempdm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);
        // matmul(transpose(sets->U_d2a), sets->U_ref, overlap);
        memset(overlap2, 0, seth->Nstate * seth->Nstate * sizeof(double));
       
        for (i = 0; i < seth->Nstate * seth->Nstate; i++){
            tempdm1[i] = fabs(overlap[i]);
        }

        for (i = 0; i < seth->Nstate; i++) {
            idloc=maxloc(tempdm1, seth->Nstate * seth->Nstate);
            id1 = idloc / seth->Nstate; 
            id2 = idloc % seth->Nstate; 
            overlap2[idloc] = (overlap[idloc] >= 0.0) ? 1.0 : -1.0;
            for (j = 0; j < seth->Nstate; j++) {
                tempdm1[id1 * seth->Nstate + j] = 0;
                tempdm1[j * seth->Nstate + id2] = 0;
            }
        }

        dd_matmul(sets->U_d2a, overlap2, tempdm1, seth->Nstate, seth->Nstate, seth->Nstate);
        memcpy(sets->U_d2a,tempdm1,seth->Nstate * seth->Nstate * sizeof(double));

        for (i = 0; i < seth->Nstate * seth->Nstate; i++){
            tempdm2[i] = fabs(overlap2[i]);
        }
        dd_matmul(sets->E_adia, tempdm2, tempdv1, 1, seth->Nstate, seth->Nstate);
        memcpy(sets->E_adia, tempdv1, seth->Nstate * sizeof(double));
    }

    if (seth->type_prop_adia > 0) {
        transpose(sets->U_d2a,tempdm1,seth->Nstate);
        dd_matmul(tempdm1, sets->U_ref, sets->overlap_adia, seth->Nstate, seth->Nstate, seth->Nstate);
    }

    memcpy(sets->U_d2a_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double));
    memcpy(sets->U_ref, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double));
    seth->if_ad_nac = 1;


    transpose(sets->U_d2a,tempdm1,seth->Nstate);
    for (i = 0; i < seth->Ndof1; i++) {
        for (j = 0; j < seth->Ndof2; j++) {
            for (int k = 0; k < seth->Nstate * seth->Nstate; k++){
                tempdm2[k] = sets->dV[k * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof2 + j];
            }
            dd_matmul(tempdm1,tempdm2,tempdm3,seth->Nstate,seth->Nstate,seth->Nstate);
            dd_matmul(tempdm3,sets->U_d2a,tempdm4,seth->Nstate,seth->Nstate,seth->Nstate);
            for (int k = 0; k < seth->Nstate * seth->Nstate; k++){
                sets->dv_adia[k * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof2 + j] = tempdm4[k];
            }
        }
    }

    memset(sets->nac, 0, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double));

    for (i = 0; i < seth->Nstate; i++) {
        for (j = 0; j < seth->Nstate; j++) {
            if (i == j) continue;
            for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++){
                sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] = sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] / (sets->E_adia[j] - sets->E_adia[i]);
            } 
        }
    }

}


// void cal_dvadia_state(){}


void cal_propagator_adia(int Nstate, double dt, double complex *U, struct set_slave *sets, struct set_host *seth){
    int i, j;
    double complex H_eff[seth->Nstate * seth->Nstate], C[seth->Nstate * seth->Nstate];
    double E[seth->Nstate], eig[seth->Nstate * seth->Nstate];
    double sineig[seth->Nstate * seth->Nstate], coseig[seth->Nstate * seth->Nstate];
    double real_pro[seth->Nstate * seth->Nstate], img_pro[seth->Nstate * seth->Nstate];
    double complex mean, rho[seth->Nstate * seth->Nstate], E_t[seth->Nstate * seth->Nstate];
    double x00[seth->Nstate], p00[seth->Nstate];
    double P_eff[seth->Ndof1 * seth->Ndof2];
    double E_avg;
    double complex eiet[seth->Nstate * seth->Nstate];
    double sum;
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate * seth->Nstate];


    // for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
    //     H_eff[i] = 0.0 + 0.0 * I;
    // }
    memset(H_eff,0,seth->Nstate * seth->Nstate * sizeof(double complex));

    for (i = 0; i < seth->Nstate; i++) {
        for (j = 0; j < seth->Nstate; j++) {
            sum=0;
            for (int k=0;k<seth->Ndof1*seth->Ndof2;k++){
                sum += sets->P_nuc[k] * sets->nac[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k]/sets->mass[k];
            }
            H_eff[i * seth->Nstate + j] = sets->E_adia[i] * Kronecker_delta(i, j) - I * sum;
        }
    }

    dia_hermitemat(seth->Nstate, H_eff, E, C);


    // int slavecore_id = athread_get_id(-1);
    // if(slavecore_id == 0) printf("V=%18.8E %18.8E %18.8E %18.8E\n",sets->V[0],sets->V[1],sets->V[2],sets->V[3]);
    // if(slavecore_id == 0) printf("E=%18.8E %18.8E\n",sets->E_adia[0],sets->E_adia[1]);
    // if(slavecore_id == 0) printf("H=%18.8E %18.8E %18.8E %18.8E\n",creal(H_eff[0]),creal(H_eff[1]),creal(H_eff[2]),creal(H_eff[3])) ;

  
    memset(sineig,0,seth->Nstate * seth->Nstate * sizeof(double));
    memset(coseig,0,seth->Nstate * seth->Nstate * sizeof(double));

    for (i = 0; i < seth->Nstate; i++) {
        // eig[i * seth->Nstate + i] = E[i];
        sineig[i * seth->Nstate + i] = sin(E[i] * dt);
        coseig[i * seth->Nstate + i] = cos(E[i] * dt);
    }

    // matmul(C, coseig - I * sineig, transpose_conjg(C), U);
    for (i = 0; i < seth->Nstate * seth->Nstate; i++){
        eiet[i] = coseig[i] - I * sineig[i];
    }
    diagger(C,tempcm1,seth->Nstate);
    cc_matmul(eiet,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
    cc_matmul(C,tempcm2,U,seth->Nstate,seth->Nstate,seth->Nstate);

    // int slavecore_id = athread_get_id(-1);
    // if(slavecore_id == 10){
    //     for(int i=0;i<seth->Nstate*seth->Nstate;i++){
    //         printf("temp=%18.8E %18.8E\n",creal(tempcm2[i]),cimag(tempcm2[i]));
    //     }
    //     for(int i=0;i<seth->Nstate*seth->Nstate;i++){
    //         printf("H=%18.8E %18.8E\n",creal(H_eff[i]),cimag(H_eff[i]));
    //     }
    //     for(int i=0;i<seth->Nstate*seth->Nstate;i++){
    //         printf("eiet=%18.8E %18.8E\n",creal(eiet[i]),cimag(eiet[i]));
    //     }
    //     for(int i=0;i<seth->Nstate*seth->Nstate;i++){
    //         printf("C=%18.8E %18.8E\n",creal(C[i]),cimag(C[i]));
    //     }
    //     for(int i=0;i<seth->Nstate*seth->Nstate;i++){
    //         printf("CD=%18.8E %18.8E\n",creal(tempcm1[i]),cimag(tempcm1[i]));
    //     }
    //     for(int i=0;i<seth->Nstate*seth->Nstate;i++){
    //         printf("U=%18.8E %18.8E\n",creal(U[i]),cimag(U[i]));
    //     }
    // }       


    if (seth->type_prop_adia == 1) {
        // for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
        //     eiet[i] = 0;
        // }
        memset(eiet,0,seth->Nstate * seth->Nstate * sizeof(double complex));
        for (i = 0; i < seth->Nstate; i++) {
            eiet[i * seth->Nstate + i] = cos(sets->E_adia[i] * dt) - I * sin(sets->E_adia[i] * dt);
        }

        switch (seth->type_algorithm) {
            case 1:
                // matmul(eiet, sets->overlap_adia, U);
                cd_matmul(eiet,sets->overlap_adia,U,seth->Nstate,seth->Nstate,seth->Nstate);
                break;
            // case 2:
            //     matmul(sets->U_d2a, eiet, U);
            //     break;
            // case 3:
            //     matmul(eiet, sets->overlap_adia, U);
            //     break;
            // case 4:
            //     if (tysets->pe_prop_4cont == 1) {
            //         matmul(sets->U_d2a, eiet, U);
            //     } else if (tysets->pe_prop_4cont == 2) {
            //         matmul(eiet, transpose(sets->U_d2a), U);
            //     }
            //     break;
        }
    }

    // int slavecore_id = athread_get_id(-1);
    // if(slavecore_id == 0) printf("U=%18.8E %18.8E %18.8E %18.8E\n",creal(U[0]),creal(U[1]),creal(U[2]),creal(U[3]));
}

// void cal_propagator_adia_unsmash(){}
// void cal_propagator_dia_unsmash(){}


void free_vari(struct set_slave *sets, struct set_host *seth) {
    free(sets->R_nuc);
    free(sets->P_nuc);
    free(sets->mass);
    free(sets->force);

    if (seth->forcetype == 1) {
        free(sets->force_nuc);
    }

    free(sets->V);
    free(sets->dV);

    free(sets->R_nuc_old);
    free(sets->P_nuc_old);
    free(sets->force_old);
    free(sets->V_old);
    free(sets->dV_old);

    free(sets->propagator);

    free(sets->timegrid);

    if (seth->if_allcf == 0) {
        if (seth->outputtype >= 0) {
            free(sets->den);
        }
        if (seth->outputtype != 0) {
            free(sets->population);
            if(seth->if_st_fb == 1){
                free(sets->pop_fb);
            }
        }
    }

    free(sets->correfun_t);
    free(sets->gamma_cv);
    free(sets->gamma_cv_old);

    switch (seth->type_evo) {
        case 0:
        case 5:
            free(sets->xe);
            free(sets->pe);
            free(sets->xe_old);
            free(sets->pe_old);
            break;
        case 1:
        case 2:
            free(sets->xe);
            free(sets->pe);
            free(sets->den_e);
            free(sets->xe_old);
            free(sets->pe_old);
            free(sets->den_e_old);
            break;
    }

     if (seth->rep == 1 || seth->sampletype == 3) {
        free(sets->U_d2a);
        free(sets->U_ref);
        free(sets->U_d2a_old);
        free(sets->U_ref_old);
        free(sets->E_adia);
        free(sets->dv_adia);
        free(sets->nac);
        free(sets->E_adia_old);
        free(sets->dv_adia_old);
        free(sets->nac_old);
        free(sets->nac_check);
        free(sets->nac_check_old);
        if (seth->type_prop_adia > 0) {
            free(sets->overlap_adia);
        }
    }

    free(sets->N_nan_sum);

}
