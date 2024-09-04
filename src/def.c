

#include <complex.h>
#include <stdint.h>
#include "def.h"
#include "constant.h"
#include "gmath.h"
#include <stdio.h>
#include <math.h>


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
double *mpi_pop_fb; // 3D double array: [size1][size2][size3]
double *mpi_R_nuc_mean; // 3D double array: [size1][size2][size3]
double *mpi_P_nuc_mean; // 3D double array: [size1][size2][size3]
double *mpi_R2_nuc_mean; // 3D double array: [size1][size2][size3]
double *mpi_P2_nuc_mean; // 3D double array: [size1][size2][size3]
int *mpi_N_nan_sum; // 1D int array: [size1]

int *count_st; // 2D int array: [size1][size2]
int *mpi_count_st; // 2D int array: [size1][size2]
int *count_pertraj; // 1D int array: [size1]

double *R_nuc_mean; // 3D double array: [size1][size2][size3]
double *P_nuc_mean; // 3D double array: [size1][size2][size3]
double *R2_nuc_mean; // 3D double array: [size1][size2][size3]
double *P2_nuc_mean; // 3D double array: [size1][size2][size3]

char filepath[200]; // Path of Model input file
char workpath[200];

int Ndof1, Ndof2;
int Nstate;
int init_occ, init_occ_adia;

double *xe, *pe; // 1D double arrays: [size1] Meyer-Miller mapping variables
double *ye, *pxe, *pye; // 1D double arrays: [size1] Li-Miller mapping variables
double complex *ce; // 1D complex array: [size1] Electronic amplitude
double complex *correfun_t; // 2D complex array: [size1][size2] Electronic reduced density matrix for evolution
double complex *gamma_cv; // 2D complex array: [size1][size2] Commutator matrix or zero-point energy parameters
double complex *den_e, *den_e4nuc; // 2D complex arrays: [size1][size2]

double *xe_mb, *pe_mb; // 2D double arrays: [size1][size2]
int *state_unsmash; // 1D int array: [size1]
int if_typemb;
double *gc_unsmash, *gp_unsmash; // 1D double arrays: [size1]

double *xe_cv, *pe_cv; // 1D double arrays: [size1]

// More variables...

int type_evo;
double gamma_zpe;
double sigma2_lsc;
double gamma1_lsc, gamma2_lsc;
double identity_M1, identity_M2;
int scheme_cvk, scheme_focus;
double alpha_gdtwa, beta_gdtwa, delta_gdtwa, eps_gdtwa;
int if_alpha;
int id_state, id_hop;
int id_state_old;
double *t_decoh, *t_coh, *L_red; // 1D double arrays: [size1]
double w_dish;

double *xe_old, *pe_old; // 1D double arrays: [size1]
double *P_nuc_old, *R_nuc_old; // 2D double arrays: [size1][size2]
double *force_old; // 2D double array: [size1][size2]
double complex *gamma_cv_old, *den_e_old; // 2D complex arrays: [size1][size2]
double *nacv_old, *V_old, *dV_old, *E_adia_old, *dv_adia_old; // Various dimensional arrays

#include <complex.h>
#include <stdbool.h>

// 动态分配的数组声明
double *nac_check_old; // 4D double array: [size1][size2][size3][size4]
double *U_ref_old; // 2D double array: [size1][size2]

double *P_nuc_old_traj; // 2D double array: [size1][size2]
double *R_nuc_old_traj; // 2D double array: [size1][size2]

double *deltaR_afssh; // 4D double array: [size1][size2][size3][size4]
double *deltaP_afssh; // 4D double array: [size1][size2][size3][size4]
double *deltaF_afssh; // 4D double array: [size1][size2][size3][size4]

int index_t0;
int index_t0_1, index_t0_2;
double complex correfun_0;

double complex *correfun_0_ms2; // 1D complex array: [size1]

double complex *cf0; // 2D complex array: [size1][size2]
double complex *cfall; // 5D complex array: [size1][size2][size3][size4][size5]
double complex *cfeff; // 1D complex array: [size1]
double *weight0; // 2D double array: [size1][size2]
double *weightt; // 2D double array: [size1][size2]
int if_allcf, allcf_times;

double complex correfun_0_pldm1, correfun_0_pldm2;
double complex *prop_pldm; // 2D complex array: [size1][size2]

double *U_d2a; // 2D double array: [size1][size2]
double *E_adia; // 1D double array: [size1]
double *nac; // 4D double array: [size1][size2][size3][size4]
double *dv_adia; // 4D double array: [size1][size2][size3][size4]
double *P_kin; // 2D double array: [size1][size2]
double *nac_check; // 4D double array: [size1][size2][size3][size4]
double *U_ref; // 2D double array: [size1][size2]
double *overlap_adia; // 2D double array: [size1][size2]
int rep; // 0 for diabatic, 1 for adiabatic

double *U_d2a_old; // 2D double array: [size1][size2]

int ifcv; // 0: no adjustment; -1: adjustment without evolution; 1: cv adjustment 
int ifid;
double complex *den; // 3D complex array: [size1][size2][size3] (total electronic reduced density matrix)
double *population; // 2D double array: [size1][size2]
double *pop_fb; // 3D double array: [size1][size2][size3]

double complex *den_traj; // 3D complex array: [size1][size2][size3]
double complex *den2_traj; // 3D complex array: [size1][size2][size3]
double *population_traj; // 2D double array: [size1][size2]
double *population2_traj; // 2D double array: [size1][size2]
double *pop_fb_traj; // 3D double array: [size1][size2][size3]
double *pop_fb2_traj; // 3D double array: [size1][size2][size3]

int if_st_fb;

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
int type_traj_sed;

double beta, temperature;

char method[20];

double dt, ttot, t_now;
double *timegrid; // 1D double array: [size1]
int Nbreak, Ngrid;
long long Ntraj;

char unit_t[20];
double unittrans_t;

int outputtype;

int calforcetype;

int ifoutputmpi;

int sampletype;

int if_st_nan;
int *N_nan_sum; // 1D int array: [size1]

int if_traj;

int type_phase;

int type_ad_fssh;

int if_ref;

int if_1st;

int if_inv_focus;

int if_Pdis, s_N;
double s_start, s_end;
double *s; // 1D double array: [size1]
double *real_expisP; // 1D double array: [size1]
double *img_expisP; // 1D double array: [size1]
double complex *expisP; // 1D complex array: [size1]
double complex *mpi_expisP; // 1D complex array: [size1]

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

int Nref;

double *pop0; // 1D double array: [size1]

bool if_ad_nac;

int mean_nuc;

bool if_occ_rand;

int if_engconsv;
double complex *engconsv_adjmat; // 2D complex array: [size1][size2]

int if_RBC;

double measure_mash;
int type_mash;
double complex rho0_mash[4], rhot_mash[4];

double U0_mash[4], mea_mat_mash[4];

double complex *rho0_unsmash; // 2D complex array: [size1][size2]
double complex *rhot_unsmash; // 2D complex array: [size1][size2]
double *U0_unsmash; // 2D double array: [size1][size2]
double *mea_mat_unsmash; // 2D double array: [size1][size2]

int ifBA;
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

int ifmodprop;
double complex *propagator_ref; // 2D complex array: [size1][size2]

int ifcorreden;

int memorylength;
int itime_save, i_re_save;

int if_st_eng;
double *energy_est; // 1D double array: [size1]
double *mpi_energy_est; // 1D double array: [size1]

int typeevo_ele;

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

int type_prop_adia;

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

int ifswitchforce, ifmashforce;

int ifscaleenergy;

double gamma_rescale;
int ifscalegamma;

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

int direc_padj;

int ifcount;

int ifreflp;

int ifreflp_mash;

int ifhardwall;

double *eig_cv; // 1D double array: [size1]
double *eig_cv_mat; // 2D double array: [size1][size2]
double complex *commu_vari; // 2D complex array: [size1][size2]

int ifzpecorr;

int iflangevin;
double eta_langevin;

double scale_sqc2;

int type_algorithm;

int type_prop_4cont;

int scaleenergy_type;

int n_step_algo5;

int allow_hop;

int if_traceless_force;

int forcetype;// should be removed!!!!!!!!!!!!!!!!!!!!
char msmodelname[200]; // should be removed!!!!!!!!!!!!!!!!!!!!


void initial_para() {
    Nbreak = 1;
    ifcv = 0;

    strcpy(unit_t, "au");
    unittrans_t = 1.0;

    temperature = 0.0;
    beta = 0.0;

    gamma_zpe = 0.0;

    index_t0 = 0;

    outputtype = 0;

    calforcetype = 0;

    ifoutputmpi = 0;

    sampletype = 0;

    // forcetype = 0;

    ifid = 0;

    if_st_nan = 0;

    scheme_focus = 0;

    rep = 0;

    if_st_fb = 0;

    if_traj = 0;

    type_evo = 0;

    type_phase = 0;

    if_alpha = 0;

    type_ad_fssh = 0;

    if_ref = 0;

    if_Pdis = 0;

    if_allcf = 0;

    allcf_times = 3;

    if_occ_rand = 0;

    mean_nuc = 0;

    if_engconsv = 0;

    if_inv_focus = 0;

    if_1st = 0;

    Nref = 1;

    if_RBC = 1;

    type_mash = 1;

    ifBA = 0;

    iflbd = 0;

    ifmfsh = 0;

    thres = 1e-8;

    ifrw = 0;

    beta_rw = -1.0;

    ifmodprop = 0;

    ifcorreden = 0;

    memorylength = 1;

    if_st_eng = 0;

    typeevo_ele = 0;

    if_switchcv = 0;

    if_inv_evo = 0;

    type_prop_adia = 0;

    thres_ms2 = 1e-30;

    ifmsbranch = 0;

    ifswitchforce = 0;

    ifmashforce = 0;

    ifscaleenergy = 0;

    ifcount = 0;

    ifreflp = 1;

    ifscalegamma = 0;

    direc_padj = 0;

    ifreflp_mash = 0;

    ifhardwall = 1;

    ifzpecorr = 0;

    iflangevin = 0;

    type_algorithm = 1;

    n_step_algo5 = 4;

    if_typemb = 0;

    allow_hop = 0;

    if_traceless_force = 0;
}



void initial_vari() {
    int i, Ngrid;

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
    // propagator_switch = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
    // if (ifmodprop == 1) {
    //     propagator_ref = (double complex *)malloc(Nstate * Nstate * sizeof(double complex ));
    // }

    Ngrid = (int)(ttot / dt) / Nbreak + 1;
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
    memset(correfun_t, 0, Nstate * Nstate * sizeof(double));
    gamma_cv = (double complex *)malloc(Nstate * Nstate * sizeof(double));
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
    // if (if_engconsv == 1) {
    //     if (type_evo != 3) type_evo = 1;
    //     engconsv_adjmat = (double *)malloc(Nstate * Nstate * sizeof(double));
    // }
    N_nan_sum = (int *)malloc(Ngrid * sizeof(int));
    mpi_N_nan_sum = (int *)malloc(Ngrid * sizeof(int));
    memset(N_nan_sum, 0, Ngrid * sizeof(int));
    if (ifswitchforce > 0) {
        if (rep == 0) {
            U_d2a = (double *)malloc(Nstate * Nstate * sizeof(double));
            E_adia = (double *)malloc(Nstate * sizeof(double));
            U_d2a_old = (double *)malloc(Nstate * Nstate * sizeof(double));
            E_adia_old = (double *)malloc(Nstate * sizeof(double));
        }
    }

}

void print_info(){
    
    printf("===========================Model details============================\n");
    printf("Model: %s\n", msmodelname);

    switch (rep) {
        case 0:
            printf("Diabatic representation will be used\n");
            break;
        case 1:
            printf("Adiabatic representation will be used\n");
            break;
    }

    if (index_t0 == 0) {
        printf("Initial occupied: %d\n", init_occ);
    } else if (index_t0 == 1) {
        printf("Detected index_t0=1\n");
        printf("Using off-diagonal term for initial condition\n");
        printf("index_t0_1=%d, index_t0_2=%d\n", index_t0_1, index_t0_2);
    }

    switch (sampletype) {
        case 0:
            break;
        case 1:
            printf("sampletype=1: initial occupied on corresponding (system) adiabatic state\n");
            break;
        case 2:
            printf("sampletype=2: initial occupied on diabatic state, but evolution in adiabatic representation\n");
            break;
        case 3:
            printf("sampletype=3: initial occupied on adiabatic state, but evolution in diabatic representation\n");
            break;
    }

    printf("=========================Simulation details=========================\n");

    if (strcmp(method, "MFT") == 0 || strcmp(method, "mft") == 0) {
        printf("Method: Ehrenfest mean-field trajectory (MFT)\n");
        printf("Related Publication: J. Chem. Phys. 2005, 123(8), 084106\n");
    } else if (strcmp(method, "BCMF") == 0 || strcmp(method, "bcmf") == 0) {
        printf("Method: branching corrected mean field (BCMF)\n");
        printf("Related Publication: J. Phys. Chem. Lett. 2020, 11, 8283-8291\n");
    } else if (strcmp(method, "FSSH") == 0 || strcmp(method, "fssh") == 0) {
        printf("Method: Fewest Switches Surface Hopping (FSSH)\n");
        printf("Related Publication: J. Chem. Phys. 1990, 93, 1061\n");
    } else if (strcmp(method, "fsshswitch") == 0) {
        printf("Method: Fewest Switches Surface Hopping (FSSH) with switch direction for adjustment P: %d\n", direc_padj);
    } else if (strcmp(method, "PCSH-NAF") == 0 || strcmp(method, "pcsh-naf") == 0) {
        printf("Method: PCSH with NAF\n");
        printf("Direction for adjustment P: %d\n", direc_padj);
    } else if (strcmp(method, "DISH") == 0 || strcmp(method, "dish") == 0) {
        printf("Method: Decoherence-induced Surface Hopping (DISH)\n");
        printf("Related Publication: J. Chem. Phys. 2012, 137, 22A545\n");
    } else if (strcmp(method, "PCSH") == 0 || strcmp(method, "pcsh") == 0) {
        printf("Method: Phase correction Surface Hopping (PCSH)\n");
        printf("Related Publication: J. Chem. Phys. 2011, 135, 024101\n");
    } else if (strcmp(method, "BCSH") == 0 || strcmp(method, "bcsh") == 0) {
        printf("Method: Branching Corrected Surface Hopping (BCSH)\n");
        printf("Related Publication: J. Chem. Phys. 2019, 150, 164101\n");
    } else if (strcmp(method, "BCSH-NAF") == 0 || strcmp(method, "bcsh-naf") == 0) {
        printf("Method: BCSH with NAF\n");
    } else if (strcmp(method, "A-FSSH") == 0 || strcmp(method, "a-fssh") == 0 || strcmp(method, "AFSSH") == 0 || strcmp(method, "afssh") == 0) {
        printf("Method: Augmented Fewest Switches Surface Hopping (A-FSSH)\n");
        printf("Related Publication: J. Chem. Phys. 2011, 134, 024105\n");
    } else if (strcmp(method, "GFSH") == 0 || strcmp(method, "gfsh") == 0) {
        printf("Method: Global Flux Surface Hopping (GFSH)\n");
        printf("Related Publication: J. Chem. Theory Comput. 2014, 10, 3598-3605\n");
    } else if (strcmp(method, "SC-FSSH") == 0 || strcmp(method, "sc-fssh") == 0) {
        printf("Method: Self-Consistent Fewest Switches Surface Hopping (SC-FSSH)\n");
        printf("Related Publication: J. Phys. Chem. Lett. 2014, 5, 713-719\n");
    } else if (strcmp(method, "CC-FSSH") == 0 || strcmp(method, "cc-fssh") == 0) {
        printf("Method: Crossing Corrected Fewest Switches Surface Hopping (CC-FSSH)\n");
        printf("Related Publication: J. Phys. Chem. Lett. 2018, 9, 4319-4325\n");
    } else if (strcmp(method, "eCMM") == 0 || strcmp(method, "ecmm") == 0) {
        printf("Method: extended Classical Mapping Model (eCMM)\n");
        printf("Related Publication: J. Chem. Phys. 2019, 151, 024105\n");
        printf("                    J. Phys. Chem. Lett. 2021, 12, 2496-2501\n");
        printf("ZPE gamma parameter: %f\n", gamma_zpe);
    } else if (strcmp(method, "focus") == 0 || strcmp(method, "efocus") == 0 || strcmp(method, "langer") == 0 || strcmp(method, "deltaMFT") == 0 || strcmp(method, "dmft") == 0) {
        printf("Method: Langer/spin-LSC focus/delta MFT\n");
        printf("Related Publication: J. Chem. Phys. 1979, 70, 3214\n");
        printf("                    J. Chem. Phys. 2020, 153, 194110\n");
        printf("                    J. Phys. Chem. Lett. in preparation\n");
        printf("ZPE gamma parameter: %f\n", gamma_zpe);
    } else if (strcmp(method, "sed2") == 0 || strcmp(method, "SED2") == 0) {
        printf("Method: SED2\n");
        printf("ZPE gamma parameter: %f\n", gamma_zpe);
    } else if (strcmp(method, "sed3") == 0 || strcmp(method, "SED3") == 0) {
        printf("Method: SED3\n");
        printf("ZPE gamma parameter: %f\n", gamma_zpe);
    } else if (strcmp(method, "sqc") == 0 || strcmp(method, "SQC") == 0) {
        printf("Method: Symmetric Quasi-Classical (SQC) approaches\n");
        printf("Related Publication: J. Chem. Phys. 2016, 145, 144108\n");
        printf("                    J. Chem. Phys. 2019, 150, 104101\n");
        printf("Using the default ZPE gamma parameter: 1/3\n");
    } else if (strcmp(method, "sqc2") == 0 || strcmp(method, "SQC2") == 0) {
        printf("Method: Symmetric Quasi-Classical approaches 2 (SQC2)\n");
        printf("ZPE gamma parameter: %f\n", gamma_zpe);
    } else if (strcmp(method, "sqc3") == 0 || strcmp(method, "SQC3") == 0) {
        printf("Method: Symmetric Quasi-Classical approaches 3 (SQC3)\n");
    } else if (strcmp(method, "gauss") == 0 || strcmp(method, "genLSC") == 0 || strcmp(method, "genlsc") == 0) {
        printf("Method: Generalized LSC (genLSC)/Gaussian sampling approaches\n");
        printf("Related Publication: in preparation\n");
        printf("sigma2=%f\n", sigma2_lsc);
        printf("gamma1=%f\n", gamma1_lsc);
        printf("gamma2=%f\n", gamma2_lsc);
        printf("default ZPE gamma parameter: 1/2\n");
    } else if (strcmp(method, "scmm") == 0) {
        printf("Method: symmetric cmm\n");
    } else if (strcmp(method, "cvk") == 0 || strcmp(method, "CVK") == 0) {
        printf("Method: commutator variables kernel approach (unpublished)\n");
    

    }else if (strcmp(method, "GDTWA") == 0 || strcmp(method, "eGDTWA") == 0) {
        printf("Method: (extend) Generalized discrete truncated Wigner approximation (GDTWA)\n");
        printf("Related Pulication: J. Chem. Phys. 2021, 155, 024111\n");
        if (if_inv_focus == 1) {
            printf("using unsymmetric kernels\n");
            printf("kernel: alpha=%f\n", alpha_gdtwa);
            printf("kernel: beta=%f\n", beta_gdtwa);
            printf("inverse: delta=%f\n", delta_gdtwa);
            printf("inverse: epsilon=%f\n", eps_gdtwa);
            printf("2*alpha*delta+(F-2)*beta*epsilon=%f\n", 2 * alpha_gdtwa * delta_gdtwa + (Nstate - 2) * beta_gdtwa * eps_gdtwa);
        }
    } else if (strcmp(method, "mash") == 0 || strcmp(method, "MASH") == 0) {
        printf("Method: Mapping Approach to Surface Hopping (MASH)\n");
        printf("Related Pulication: J. Chem. Phys. 2023, in press\n");
        if (ifreflp_mash == 1) printf("ifreflp_mash=%d\n", ifreflp_mash);
    } else if (strcmp(method, "ms-mash") == 0 || strcmp(method, "MS-MASH") == 0) {
        printf("Method: Multi-State Mapping Approach to Surface Hopping (MS-MASH)\n");
        printf("Related Pulication: arXiv:2305.08835\n");
        if (ifreflp_mash == 1) printf("ifreflp_mash=%d\n", ifreflp_mash);
    } else if (strcmp(method, "unsmash") == 0 || strcmp(method, "UNSMASH") == 0 || strcmp(method, "unSMASH") == 0) {
        printf("Method: uncoupled spheres Mapping Approach to Surface Hopping (unSMASH)\n");
        printf("Related Pulication: arXiv:2403.10627\n");
        // if (ifreflp_mash == 1) printf("ifreflp_mash=%d\n", ifreflp_mash);
    } else if (strcmp(method, "unsmash-mf") == 0 || strcmp(method, "UNSMASH-MF") == 0 || strcmp(method, "unSMASH-MF") == 0) {
        printf("Method: uncoupled spheres Mapping Approach to Surface Hopping with mean field (unSMASH-MF)\n");
    } else if (strcmp(method, "ms-mash-focus") == 0 || strcmp(method, "MS-MASH-focus") == 0) {
        printf("Method: Multi-State Mapping Approach to Surface Hopping with focus initial condition (MS-MASH-focus)\n");
        printf("Related Pulication: arXiv:2305.08835\n");
    } else if (strcmp(method, "ms-mash-mf") == 0 || strcmp(method, "MS-MASH-MF") == 0) {
        printf("Method: Multi-State Mapping Approach to Surface Hopping with mean field trajectories (MS-MASH-MF)\n");
        printf("Related Pulication: arXiv:2305.08835\n");
    } else if (strcmp(method, "ms-mash-mf2") == 0 || strcmp(method, "MS-MASH-MF2") == 0) {
        printf("Method: Multi-State Mapping Approach to Surface Hopping with mean field trajectories version 2 (MS-MASH-MF2)\n");
    } else if (strcmp(method, "mf2-sh") == 0 || strcmp(method, "MF2-SH") == 0) {
        printf("Method: MF2-SH\n");
    } else if (strcmp(method, "mf2cv") == 0 || strcmp(method, "MF2CV") == 0) {
        printf("Method: MF2CV\n");
    } else if (strcmp(method, "dmf2") == 0 || strcmp(method, "DMF2") == 0) {
        printf("Method: delta MF2\n");
        printf("gamma_zpe=%f\n", gamma_zpe);
    } else if (strcmp(method, "mf3") == 0 || strcmp(method, "MF3") == 0) {
        printf("Method: MF3\n");
        printf("gamma_zpe=%f\n", gamma_zpe);
    } else if (strcmp(method, "wsqc") == 0 || strcmp(method, "WSQC") == 0) {
        printf("Method: wSQC\n");
        printf("gamma_zpe=%f\n", gamma_zpe);
    } else if (strcmp(method, "ms2") == 0 || strcmp(method, "MS2") == 0) {
        printf("Method: MS2\n");
        printf("thres_ms2=%f\n", thres_ms2);
    } else if (strcmp(method, "mstraj") == 0) {
        printf("Method: mstraj\n");
        printf("thres_ms2=%f\n", thres_ms2);
    } else if (strcmp(method, "mstraj2") == 0) {
        printf("Method: mstraj2\n");
        printf("thres_ms2=%f\n", thres_ms2);
    } else if (strcmp(method, "ms3") == 0 || strcmp(method, "MS3") == 0) {
        printf("Method: MS3\n");
    } else if (strcmp(method, "mash-mf") == 0 || strcmp(method, "MASH-MF") == 0) {
        printf("Method: Mapping Approach to Surface Hopping with mean field trajectories (MASH-MF)\n");
        printf("Related Pulication: J. Chem. Phys. 2023, in press\n");
    } else if (strcmp(method, "mash-mf3") == 0 || strcmp(method, "MASH-MF3") == 0) {
        printf("Method: Mapping Approach to Surface Hopping with mean field trajectories version 3 (MASH-MF3)\n");
        printf("Related Pulication: J. Chem. Phys. 2023, in press\n");
    } else if (strcmp(method, "cvsh") == 0 || strcmp(method, "CVSH") == 0) {
        printf("Method: commutator variable surface hopping\n");
        printf("scheme_cvsh=%d\n", scheme_cvsh);
        if (scheme_cvsh == 1 || scheme_cvsh == 2 || scheme_cvsh == 3 || scheme_cvsh == 4) printf("ZPE gamma parameter: %f\n", gamma_zpe);
        if (scheme_cvsh == 9 || scheme_cvsh == 10 || scheme_cvsh == 11 || scheme_cvsh == 12) printf("alpha=%f, beta=%f\n", alpha_gdtwa, beta_gdtwa);
        if (scheme_cvsh == 13 || scheme_cvsh == 14 || scheme_cvsh == 15 || scheme_cvsh == 16) printf("scheme_cvk=%d\n", scheme_cvk);
    } else if (strcmp(method, "nmsse") == 0) {
        printf("Method: NMSSE\n");
    } else if (strcmp(method, "rdmfocus") == 0) {
        printf("Method: random focus\n");
    } else if (strcmp(method, "switch") == 0) {
        printf("Method: switch (unfinished)\n");
    } else if (strcmp(method, "focus_jump") == 0) {
        printf("Method: focus_jump, type_jump=%d\n", type_jump);
    } else if (strcmp(method, "test") == 0) {
        printf("!@TestTestTestTestTestTest===test module===TestTestTestTestTestTest@!\n");
    }


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
    
    // tempmat1[Nstate*Nstate],tempmat2[Nstate*Nstate],tempmat3[Nstate*Nstate]

    if (sampletype == 2 && typeevo_ele == 0) {
        // xe = matmul(U_d2a, xe)
        // matmul(U_d2a, xe, Nstate, Nstate, 1);
        double tempv[Nstate];
        memcpy(tempv,xe,Nstate*sizeof(double));
        dd_matmul(U_d2a,tempv,xe,Nstate,Nstate,1);
        // pe = matmul(U_d2a, pe)
        // matmul(U_d2a, pe, Nstate, Nstate, 1);
        memcpy(tempv,pe,Nstate*sizeof(double));
        dd_matmul(U_d2a,tempv,pe,Nstate,Nstate,1);
        // gamma_cv = matmul(matmul(U_d2a, gamma_cv), transpose(U_d2a))
        // matmul(U_d2a, gamma_cv, Nstate, Nstate, Nstate);
        double complex tempcm[Nstate*Nstate];
        memcpy(tempcm,gamma_cv,Nstate*Nstate*sizeof(double complex));
        double tempdm[Nstate*Nstate];
        transpose(U_d2a,tempdm,Nstate);
        double complex tempcm2[Nstate*Nstate];
        cd_matmul(gamma_cv,tempdm,tempcm2,Nstate,Nstate,Nstate);
        dc_matmul(U_d2a,tempcm2,gamma_cv,Nstate,Nstate,Nstate);

        
        // if (type_evo == 1 || type_evo == 3) {
            // den_e = matmul(matmul(U_d2a, den_e), transpose(U_d2a))
            // matmul(U_d2a, den_e, Nstate, Nstate, Nstate);
            // transpose(U_d2a, Nstate, Nstate);
            // matmul(den_e, U_d2a, Nstate, Nstate, Nstate);
        // }
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
    double complex tempcm1[Nstate*Nstate],tempcm2[Nstate*Nstate];
    // if (U_d2a_old != NULL) memcpy(U_d2a_old, U_d2a, sizeof(U_d2a));

    switch (type_evo) {
        case 0:{
            if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
            // if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
            memcpy(x0, xe, Nstate * sizeof(double));
            memcpy(p0, pe, Nstate * sizeof(double));
            // matmul_real_imag(Nstate, propagator, x0, p0, xe, pe);
            double complex tempv1[Nstate],tempv2[Nstate];
            cd_matmul(propagator,x0,tempv1,Nstate,Nstate,1);
            cd_matmul(propagator,p0,tempv2,Nstate,Nstate,1);
            for(i=0;i<Nstate;i++){
                xe[i]=creal(tempv1[i])-cimag(tempv2[i]);
                pe[i]=creal(tempv2[i])+cimag(tempv1[i]);
            }
            break;
        }
        case 1:{
            if (rep == 0) cal_propagator(Nstate, V, deltat, propagator);
            // if (rep == 1 && typeevo_ele == 0) cal_propagator_adia(Nstate, deltat, propagator);
            // matmul_complex(Nstate, propagator, den_e, den_e);
            
            diagger(propagator,tempcm1,Nstate);
            cc_matmul(den_e,tempcm1,tempcm2,Nstate,Nstate,Nstate);
            cc_matmul(propagator,tempcm2,den_e,Nstate,Nstate,Nstate);
            // if (inverse_kernel != NULL) {
            //     matmul_complex(Nstate, propagator, inverse_kernel, inverse_kernel);
            // }
            break;
        }
            
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

    if (den != NULL) {
        if (strcmp(trim(adjustl(method)), "gauss") == 0 || strcmp(trim(adjustl(method)), "genLSC") == 0 || strcmp(trim(adjustl(method)), "genlsc") == 0) {
            if (ifid == 0) {
                for (i = 0; i < Nstate; i++) {
                    for (j = 0; j < Nstate; j++) {
                        den[i * Nstate + j + igrid_cal * Nstate * Nstate] += correfun_0 * correfun_t[i * Nstate + j];
                    }
                }
            } else if (ifid == 1) {
                for (i = 0; i < Nstate; i++) {
                    for (j = 0; j < Nstate; j++) {
                        den[i * Nstate + j + igrid_cal * Nstate * Nstate] += correfun_0 * correfun_t[i * Nstate + j];
                    }
                }
                for (j = 0; j < Nstate; j++) {
                    den[j * Nstate + j + igrid_cal * Nstate * Nstate] += 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                }
            }
        } else {
            for (i = 0; i < Nstate; i++) {
                for (j = 0; j < Nstate; j++) {
                    den[i * Nstate + j + igrid_cal * Nstate * Nstate] += correfun_0 * correfun_t[i * Nstate + j];
                }
            }
        }
    }

    if (population!=NULL) {
        for (i = 0; i < Nstate; i++) {
            if (strcmp(trim(adjustl(method)), "gauss") == 0 || strcmp(trim(adjustl(method)), "genLSC") == 0 || strcmp(trim(adjustl(method)), "genlsc") == 0) {
                if (ifid == 0) {
                    population[i + igrid_cal * Nstate] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                    if (if_st_fb == 1) {
                        if (P_nuc[0] > 0) {
                            pop_fb[i * 2 + igrid_cal * Nstate * 2] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                        } else {
                            pop_fb[i * 2 + 1 + igrid_cal * Nstate * 2] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                        }
                    }
                } else if (ifid == 1) {
                    population[i + igrid_cal * Nstate] += creal(correfun_0 * correfun_t[i * Nstate + i]) + 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                    if (if_st_fb == 1) {
                        if (P_nuc[0] > 0) {
                            pop_fb[i * 2 + igrid_cal * Nstate * 2] += creal(correfun_0 * correfun_t[i * Nstate + i]) + 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                        } else {
                            pop_fb[i * 2 + 1 + igrid_cal * Nstate * 2] += creal(correfun_0 * correfun_t[i * Nstate + i]) + 1.0 / Nstate - 1.0 / (sigma2_lsc * sigma2_lsc * Nstate * Nstate);
                        }
                    }
                }
            } else {
                population[i + igrid_cal * Nstate] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                if (if_st_fb == 1) {
                    if (P_nuc[0] > 0) {
                        pop_fb[i * 2 + igrid_cal * Nstate * 2] += creal(correfun_0 * correfun_t[i * Nstate + i]);
                    } else {
                        pop_fb[i * 2 + 1 + igrid_cal * Nstate * 2] += creal(correfun_0 * correfun_t[i * Nstate + i]);
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
    evo_traj_nucP(deltat / 2);
    if (scaleenergy_type == 3) energy_conserve_naf_3(deltat / 2);
    if (strcmp(trim(adjustl(msmodelname)), "LZ") == 0) {
        memcpy(P_nuc, Pinit_LZ, sizeof(double));
    }
    evo_traj_nucR(deltat);
    dV_msmodel(R_nuc, dV);
    V_msmodel(R_nuc, V, t_now);
    if (rep == 1) cal_NACV();
    evo_traj_ele(deltat);
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

