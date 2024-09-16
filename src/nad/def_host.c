#include <complex.h>
#include <stdint.h>
// #include "def.h"
#include "constant.h"
#include "gmath.h"
#include <stdio.h>
#include <math.h>
#include "cJSON.h"
#include "msmodel.h"
#include <stdbool.h>
#include <iofun.h>

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
unsigned long long *mpi_N_nan_sum; // 1D int array: [size1]
double *mpi_real_den;
double *mpi_imag_den;
double *mpi_real_cfeff;
double *mpi_imag_cfeff;




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

char *filepath; // Path of Model input file
char *workpath;

int Ndof1, Ndof2;
int Nstate;
int init_occ, init_occ_adia;

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

bool if_ad_nac;

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

double scale_sqc2;

int type_algorithm;

int type_prop_4cont;

int scaleenergy_type;

int n_step_algo5;

int allow_hop;

int if_traceless_force;



// int forcetype;
// char msmodelname[200];




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

    forcetype = 0;

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

    // iflbd = 0;

    // ifmfsh = 0;

    // thres = 1e-8;

    // ifrw = 0;

    // beta_rw = -1.0;

    ifmodprop = 0;

    ifcorreden = 0;

    // memorylength = 1;

    if_st_eng = 0;

    typeevo_ele = 0;

    // if_switchcv = 0;

    // if_inv_evo = 0;

    type_prop_adia = 0;

    // thres_ms2 = 1e-30;

    // ifmsbranch = 0;

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




void readinp(){

    

    FILE *file = fopen(filepath, "rb");
    if (!file) {
        perror("File opening failed");
        return NULL;
    }
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    char *data = (char *)malloc(length + 1);
    if (data) {
        fread(data, 1, length, file);
        data[length] = '\0';
    }
    fclose(file);
    
    cJSON *json = cJSON_Parse(data);
    if (!json) {
        printf("Error parsing JSON: %s\n", cJSON_GetErrorPtr());
        return;
    }

    cJSON *item = json;
    cJSON *list;
        // cJSON *msmodelname = cJSON_GetObjectItem(item, "msmodelname");
    // if (cJSON_GetObjectItem(item, "msmodelname")!=NULL) {
    //     list=cJSON_GetObjectItem(item, "msmodelname");
    //     strcpy(msmodelname, list->valuestring); // 假设msmodelname_var在def.h中声明为char msmodelname_var[256];
    //     printf(msmodelname,"\n");
    //     item = item->next;
    // }

    if (NULL !=  cJSON_GetObjectItem(item, "msmodelname")){
        list=cJSON_GetObjectItem(item,  "msmodelname");
        strcpy(msmodelname, list->valuestring);
    }

        
        
    while (item) {
        // printf("F1=%d\n",Nstate);

        readinp_msmodel(item, &Ndof1, &Ndof2, &Nstate);
        
        // printf("F2=%d\n",Nstate);

        if (NULL != cJSON_GetObjectItem(item, "beta")) {
            list = cJSON_GetObjectItem(item, "beta");
            if (list->type == cJSON_Number) {
                beta = list->valuedouble; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "ntraj")) {
            list = cJSON_GetObjectItem(item, "ntraj");
            if (list->type == cJSON_Number) {
                Ntraj= list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "dt")) {
            list = cJSON_GetObjectItem(item, "dt");
            if (list->type == cJSON_Number) {
                dt = list->valuedouble;
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "ttot")) {
            list = cJSON_GetObjectItem(item, "ttot");
            if (list->type == cJSON_Number) {
                ttot = list->valuedouble; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "nbreak")) {
            list = cJSON_GetObjectItem(item, "nbreak");
            if (list->type == cJSON_Number) {
                Nbreak = list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "method")) {
            list = cJSON_GetObjectItem(item, "method");
            if (list->type == cJSON_String) {
                strcpy(method, list->valuestring); 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "init_occ")) {
            list = cJSON_GetObjectItem(item, "init_occ");
            if (list->type == cJSON_Number) {
                init_occ = list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "rep")) {
            list = cJSON_GetObjectItem(item, "rep");
            if (list->type == cJSON_Number) {
                rep = list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "outputtype")) {
            list = cJSON_GetObjectItem(item, "outputtype");
            if (list->type == cJSON_Number) {
                outputtype = list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "forcetype")) {
            list = cJSON_GetObjectItem(item, "forcetype");
            if (list->type == cJSON_Number) {
                forcetype = list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "calforcetype")) {
            list = cJSON_GetObjectItem(item, "calforcetype");
            if (list->type == cJSON_Number) {
                calforcetype = list->valueint; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "ifoutputmpi")) {
            list = cJSON_GetObjectItem(item, "ifoutputmpi");
            if (list->type == cJSON_Number) {
                ifoutputmpi = list->valueint; 
            }
        }

        // if (NULL != cJSON_GetObjectItem(item, "nproc_sw")) {
        //     list = cJSON_GetObjectItem(item, "nproc_sw");
        //     if (list->type == cJSON_Number) {
        //         nproc_sw = list->valueint; 
        //     }
        // }

        item = item->next;
    }

    cJSON_Delete(json);
    // free(item);
    // free(list);
    // printf("msmodelname: %s\n",msmodelname);

   
    // printf("beta: %f\n", beta);
    // printf("Ntraj: %d\n", Ntraj);
    // printf("dt: %f\n", dt);
    // printf("ttot: %f\n", ttot);
    // printf("Nbreak: %d\n", Nbreak);
    // printf("method: %s\n", method);
    // printf("init_occ: %d\n", init_occ);
    // printf("rep: %d\n", rep);
    // printf("outputtype: %d\n", outputtype);
    // printf("forcetype: %d\n", forcetype);
    // printf("calforcetype: %d\n", calforcetype);

}



void init_host(){

    fi_time_grid = (double *)malloc(Ngrid * sizeof(double));

    mpi_N_nan_sum = (unsigned long long *)malloc(Ngrid * sizeof(unsigned long long));
    memset(mpi_N_nan_sum, 0, Ngrid * sizeof(unsigned long long));
    
    if (if_allcf == 0) {
        if (outputtype >= 0) {
        // mpi_den = (double complex *)malloc(Nstate * Nstate * Ngrid * sizeof(double complex));
        // memset(mpi_den, 0, Nstate * Nstate * Ngrid * sizeof(double complex));

        mpi_real_den = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        mpi_imag_den = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));  
        memset(mpi_real_den, 0, Nstate * Nstate * Ngrid * sizeof(double));
        memset(mpi_imag_den, 0, Nstate * Nstate * Ngrid * sizeof(double));   
        }

        if (outputtype != 0){
            mpi_population = (double *)malloc(Nstate * Ngrid * sizeof(double));
            memset(mpi_population,0, Nstate * Ngrid * sizeof(double));
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
    // } else if (strcmp(method, "ms2") == 0 || strcmp(method, "MS2") == 0) {
    //     printf("Method: MS2\n");
    //     printf("thres_ms2=%f\n", thres_ms2);
    // } else if (strcmp(method, "mstraj") == 0) {
    //     printf("Method: mstraj\n");
    //     printf("thres_ms2=%f\n", thres_ms2);
    // } else if (strcmp(method, "mstraj2") == 0) {
    //     printf("Method: mstraj2\n");
    //     printf("thres_ms2=%f\n", thres_ms2);
    } else if (strcmp(method, "ms3") == 0 || strcmp(method, "MS3") == 0) {
        printf("Method: MS3\n");
    } else if (strcmp(method, "mash-mf") == 0 || strcmp(method, "MASH-MF") == 0) {
        printf("Method: Mapping Approach to Surface Hopping with mean field trajectories (MASH-MF)\n");
        printf("Related Pulication: J. Chem. Phys. 2023, in press\n");
    } else if (strcmp(method, "mash-mf3") == 0 || strcmp(method, "MASH-MF3") == 0) {
        printf("Method: Mapping Approach to Surface Hopping with mean field trajectories version 3 (MASH-MF3)\n");
        printf("Related Pulication: J. Chem. Phys. 2023, in press\n");
    // } else if (strcmp(method, "cvsh") == 0 || strcmp(method, "CVSH") == 0) {
    //     printf("Method: commutator variable surface hopping\n");
    //     printf("scheme_cvsh=%d\n", scheme_cvsh);
    //     if (scheme_cvsh == 1 || scheme_cvsh == 2 || scheme_cvsh == 3 || scheme_cvsh == 4) printf("ZPE gamma parameter: %f\n", gamma_zpe);
    //     if (scheme_cvsh == 9 || scheme_cvsh == 10 || scheme_cvsh == 11 || scheme_cvsh == 12) printf("alpha=%f, beta=%f\n", alpha_gdtwa, beta_gdtwa);
    //     if (scheme_cvsh == 13 || scheme_cvsh == 14 || scheme_cvsh == 15 || scheme_cvsh == 16) printf("scheme_cvk=%d\n", scheme_cvk);
    } else if (strcmp(method, "nmsse") == 0) {
        printf("Method: NMSSE\n");
    } else if (strcmp(method, "rdmfocus") == 0) {
        printf("Method: random focus\n");
    } else if (strcmp(method, "switch") == 0) {
        printf("Method: switch (unfinished)\n");
    // } else if (strcmp(method, "focus_jump") == 0) {
    //     printf("Method: focus_jump, type_jump=%d\n", type_jump);
    } else if (strcmp(method, "test") == 0) {
        printf("!@TestTestTestTestTestTest===test module===TestTestTestTestTestTest@!\n");
    }

    printf("--------------------------------------------------------------------\n");
        switch (ifcv) {
            case 0:
                // printf("ifcv=0: original version\n");
                break;
            case 1:
                printf("ifcv=1: Commutator Variables will be used for adjustment\n");
                printf("Related Publication: J. Phys. Chem. A 2021, 125(31), 6845-6863\n");
                break;
            case -1:
                printf("ifcv=-1: time-independent ZPE adjustment will be used\n");
                printf("Related Publication: J. Chem. Phys. 2019, 150, 194110\n");
                break;
            case 2:
                printf("ifcv=2: Commutator Variables will be used for adjustment\n");
                printf("But using GDTWA form\n");
                break;
            case -2:
                printf("ifcv=-2: time-independent ZPE adjustment will be used\n");
                printf("But using GDTWA form\n");
                break;
        }
        printf("--------------------------------------------------------------------\n");
        switch (type_evo) {
            case 1:
                printf("type_evo=1: using density form for evolution\n");
                break;
            case 2:
                printf("type_evo=2: x,p and den_e\n");
                break;
            case 3:
                printf("type_evo=3: den_e and den_e4nuc\n");
                break;
            case 4:
                printf("type_evo=4: x,p and G\n");
                break;
            case 5:
                printf("type_evo=5: x_mb, p_mb\n");
                break;
        }

        switch (if_ref) {
            case 1:
                printf("if_ref=1: Reference Twins Trajectory will be used\n");
                break;
        }
        switch (ifBA) {
            case 1:
                printf("ifBA=1: Using approximated Baeck-An nonadiabatic coupling\n");
                break;
        }
        printf("--------------------------------------------------------------------\n");
        printf("Number of trajectories: %d\n", Ntraj);
        if (strcmp(unit_t, "au") == 0) {
            printf("total time: %f au\n", ttot);
            printf("time step: %f au\n", dt);
        } else {
            printf("total time: %f %s = %f au\n", ttot / unittrans_t, unit_t, ttot);
            printf("time step: %f %s = %f au\n", dt / unittrans_t, unit_t, dt);
        }
        printf("Nbreak: %d; number of grids: %d\n", Nbreak, Ngrid);
        printf("beta= %f  (%f K)\n", beta, 1.0 / (kb * beta));
        printf("=====================================================================\n");

        if (calforcetype == 1) {
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            printf("         Warning:   Detected calforcetype = 1  \n");
            printf("In this case, only diagonal nuclear force will be calculated.\n");
            printf("It should only be used for the model whose electronic off-diagonal\n");
            printf("elements are not depend on the nuclear DOFs. (e.g., spin-boson model\n");
            printf("or site-exciton model)\n");
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        }
        if (if_st_nan == 1) {
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            printf("         Warning:   Detected if_st_nan = 1  \n");
            printf("unfinished  \n");
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        }
        if (if_allcf == 1) {
            printf("if_allcf = 1: All time correlation functions will be calculated.\n");
        } else if (if_allcf == 2) {
            printf("if_allcf = %d: All time correlation functions will be calculated,\n", if_allcf);
            printf("But only a weighted effective time correlation function will be saved.\n");
        }
        if (mean_nuc == 1) {
            printf("mean_nuc = 1: Output mean nuclear DOFs.\n");
        }
        if (if_st_eng == 1) {
            printf("if_st_eng = 1: estimate total energy \n");
        }
        if (if_engconsv == 1) {
            printf("if_engconsv = 1: Adjustment for Energy Conservation \n");
        }
        // if (iflbd > 0) {
        //     printf("iflbd = %d: Lindblad equation for evolution\n", iflbd);
        //     printf("rate parameter= %f\n", rate_para_lbd);
        // }
        // if (ifmfsh > 0) {
        //     printf("ifmfsh = %d: mean-field surface hopping\n", ifmfsh);
        //     printf("memory length = %d\n", memorylength);
        // }
        // if (ifrw > 0) {
        //     printf("ifrw = %d: rewighting adjustment \n", ifrw);
        //     printf("beta_rw= %f  (%f K)\n", beta_rw, 1.0 / (kb * beta_rw));
        // }
        if (typeevo_ele > 0) {
            printf("typeevo_ele= %d\n", typeevo_ele);
        }
        if (type_prop_adia > 0) {
            printf("type_prop_adia= %d\n", type_prop_adia);
        }
        // if (ifmsbranch > 0) {
        //     printf("msbranch, ifmsbranch= %d\n", ifmsbranch);
        // }
        if (ifswitchforce > 0) {
            printf("ifswitchforce= %d\n", ifswitchforce);
            printf("ifreflp= %d\n", ifreflp);
            if (ifscalegamma == 1) printf("scale gamma to %f\n", gamma_rescale);
            printf("direction for adjustment P: %d\n", direc_padj);
        }
        if (ifmashforce > 0) {
            printf("ifmashforce= %d\n", ifmashforce);
            printf("ifreflp= %d\n", ifreflp);
            if (ifscalegamma == 1) printf("scale gamma to %f\n", gamma_rescale);
            printf("direction for adjustment P: %d\n", direc_padj);
        }
        if (ifscaleenergy > 0) {
            printf("ifscaleenergy= %d\n", ifscaleenergy);
        }
        if (ifcount == 1) {
            printf("ifcount= %d\n", ifcount);
        }
        if (ifzpecorr > 0) {
            printf("ZPE correction, ifzpecorr= %d\n", ifzpecorr);
        }
        if (iflangevin == 1) {
            printf("Using Langevin dynamics, eta= %f\n", eta_langevin);
        }
        printf("index of algorithm: %d\n", type_algorithm);
        if (type_algorithm == 5) printf("n_step_algo5= %d\n", n_step_algo5);
        if (allow_hop != 0) {
            printf("allow_hop= %d\n", allow_hop);
        }
        if (if_traceless_force != 0) {
            printf("if_traceless_force= %d\n", if_traceless_force);
        }

}






void fileout() {
    int i, totn;
    char outname[256];


    if (if_allcf == 0) {
        printf("output type= %d\n", outputtype);
        if (outputtype == 0) {
            printf("Density matrix will be given in *.den.\n");
        } else if (outputtype > 0) {
            printf("Both density matrix and population data will be given in *.den and *.pop, respectively.\n");
        } else if (outputtype < 0) {
            printf("Only population data will be given in *.pop.\n");
        }
    } else if (if_allcf == 1) {
        printf("if_allcf = 1: All time correlation functions will be given in *.cf\n");
    } else if (if_allcf == 2 || if_allcf == 3) {
        printf("if_allcf = %d: Effective weighted correlation function will be given in *.cfeff\n", if_allcf);
    }

    size_t len = strlen(filepath);
    

    if (mpi_real_den != NULL) {
        strncpy(outname, filepath, len - 5);
        strcpy(outname + len - 5, ".den");
        FILE *den_file = fopen(outname, "w");
        totn = 2 * Nstate * Nstate + 1;
        for (i = 0; i < Ngrid; i++) {
            fprintf(den_file, "%18.8E", fi_time_grid[i] / unittrans_t);
            for (int j = 0; j < Nstate * Nstate; j++) {
                fprintf(den_file, "%18.8E", mpi_real_den[j * Ngrid + i]);
            }
            for (int j = 0; j < Nstate * Nstate; j++) {
                fprintf(den_file, "%18.8E", mpi_imag_den[j * Ngrid + i]);
            }
            fprintf(den_file, "\n");
        }
        fclose(den_file);
    }

    if (mpi_population != NULL) {
        strncpy(outname, filepath, len - 5);
        strcpy(outname + len - 5, ".pop");
        FILE *pop_file = fopen(outname, "w");
        totn = Nstate + 1;
        for (i = 0; i < Ngrid; i++) {
            fprintf(pop_file, "%18.8E", fi_time_grid[i] / unittrans_t);
            for (int j = 0; j < Nstate; j++) {
                fprintf(pop_file, "%18.8E", mpi_population[j*Ngrid+i]);
            }
            fprintf(pop_file, "\n");
        }
        fclose(pop_file);
    }

//     if (if_st_nan == 1) {
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".nfailtraj");
//         printf("Number of failed trajectories will be given in *.nfailtraj\n");
//         FILE *nfailtraj_file = fopen(outname, "w");
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(nfailtraj_file, "%18.8E %d\n", timegrid[i] / unittrans_t, N_nan_sum[i]);
//         }
//         fclose(nfailtraj_file);
//     }

//     if (if_st_fb == 1) {
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".fbpop");
//         printf("Population for forward/backward trajectories will be given in *.fbpop\n");
//         FILE *fbpop_file = fopen(outname, "w");
//         totn = Nstate * 2 + 1;
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(fbpop_file, "%18.8E", timegrid[i] / unittrans_t);
//             for (int j = 0; j < Nstate; j++) {
//                 fprintf(fbpop_file, "%18.8E%18.8E", pop_fb[j * Ngrid * 2 + i * 2], pop_fb[j * Ngrid * 2 + i * 2 + 1]);
//             }
//             fprintf(fbpop_file, "\n");
//         }
//         fclose(fbpop_file);
//     }

//     // if (cfall != NULL) {
//     //     FILE *cf_file = fopen(strcat(filepath, ".cf"), "w");
//     //     totn = 2 * Nstate * Nstate * Nstate * Nstate + 1;
//     //     for (i = 0; i < Ngrid; i++) {
//     //         fprintf(cf_file, "%18.8E", time[i] / unittrans_t);
//     //         for (int j = 0; j < Nstate * Nstate * Nstate * Nstate; j++) {
//     //             fprintf(cf_file, "%18.8E%18.8E", creal(cfall[j * Ngrid + i]), cimag(cfall[j * Ngrid + i]));
//     //         }
//     //         fprintf(cf_file, "\n");
//     //     }
//     //     fclose(cf_file);
//     // }

//     if (cfeff != NULL) {
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5,".cfeff");
//         FILE *cfeff_file = fopen(outname, "w");
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(cfeff_file, "%18.8E%18.8E%18.8E\n", timegrid[i] / unittrans_t, creal(cfeff[i]), cimag(cfeff[i]));
//         }
//         fclose(cfeff_file);
//     }

//     if (P_nuc_mean != NULL) {
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".Pmean");
//         FILE *Pmean_file = fopen(outname, "w");
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".Rmean");
//         FILE *Rmean_file = fopen(outname, "w");
//         totn = Ndof1 * Ndof2 + 1;
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(Pmean_file, "%18.8E", timegrid[i] / unittrans_t);
//             fprintf(Rmean_file, "%18.8E", timegrid[i] / unittrans_t);
//             for (int j = 0; j < Ndof1 * Ndof2; j++) {
//                 fprintf(Pmean_file, "%18.8E", P_nuc_mean[j * Ngrid + i]);
//                 fprintf(Rmean_file, "%18.8E", R_nuc_mean[j * Ngrid + i]);
//             }
//             fprintf(Pmean_file, "\n");
//             fprintf(Rmean_file, "\n");
//         }
//         fclose(Pmean_file);
//         fclose(Rmean_file);

//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".P2mean");
//         FILE *P2mean_file = fopen(outname, "w");
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".R2mean");
//         FILE *R2mean_file = fopen(outname, "w");
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(P2mean_file, "%18.8E", timegrid[i] / unittrans_t);
//             fprintf(R2mean_file, "%18.8E", timegrid[i] / unittrans_t);
//             for (int j = 0; j < Ndof1 * Ndof2; j++) {
//                 fprintf(P2mean_file, "%18.8E", P2_nuc_mean[j * Ngrid + i]);
//                 fprintf(R2mean_file, "%18.8E", R2_nuc_mean[j * Ngrid + i]);
//             }
//             fprintf(P2mean_file, "\n");
//             fprintf(R2mean_file, "\n");
//         }
//         fclose(P2mean_file);
//         fclose(R2mean_file);
//     }

//     if (energy_est != NULL) {
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5, ".energy");
//         FILE *energy_file = fopen(outname, "w");
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(energy_file, "%18.8E %18.8E\n", timegrid[i] / unittrans_t, energy_est[i]);
//         }
//         fclose(energy_file);
//     }

//     if (count_st != NULL) {
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5,  ".count");
//         FILE *count_file = fopen(outname, "w");
//         for (i = 0; i < Ngrid; i++) {
//             fprintf(count_file, "%18.8E", timegrid[i] / unittrans_t);
//             for (int j = 0; j < 5; j++) {
//                 fprintf(count_file, " %10d", count_st[j * Ngrid + i]);
//             }
//             fprintf(count_file, "\n");
//         }
//         fclose(count_file);
//     }

//     if (if_Pdis == 1) {   
//         strncpy(outname, filepath, len - 5);
//         strcpy(outname + len - 5,  ".expisp");
//         FILE *file = fopen(outname, "w");
//         for (int i = 0; i < s_N; i++) {
//             fprintf(file, "%f %f %f\n", s[i], real_expisP[i], imag_expisP[i]);
//         }
//         fclose(file);
//     }
}






void fileout_mpi(int id) {
    int i, totn;
    char outname[256];
    char cid[20];

    snprintf(cid, sizeof(cid), "%d", id);

    // if (if_allcf == 0) {
    //     printf("output type= %d\n", outputtype);
    //     if (outputtype == 0) {
    //         printf("Density matrix will be given in *.den.\n");
    //     } else if (outputtype > 0) {
    //         printf("Both density matrix and population data will be given in *.den and *.pop, respectively.\n");
    //     } else if (outputtype < 0) {
    //         printf("Only population data will be given in *.pop.\n");
    //     }
    // } else if (if_allcf == 1) {
    //     printf("if_allcf = 1: All time correlation functions will be given in *.cf\n");
    // } else if (if_allcf == 2 || if_allcf == 3) {
    //     printf("if_allcf = %d: Effective weighted correlation function will be given in *.cfeff\n", if_allcf);
    // }

    size_t len = strlen(filepath);
    

    if (mpi_real_den != NULL) {
        strncpy(outname, filepath, len - 5);
        // strcpy(outname + len - 5, ".den");
        outname[len - 5] = '\0'; // 确保字符串以null结尾
        strcpy(outname + len - 5, "_mpi");
        strcpy(outname + len - 5 + strlen("_mpi"), cid);
        strcpy(outname + len - 5 + strlen("_mpi") + strlen(cid), ".den");
        FILE *den_file = fopen(outname, "w");
        totn = 2 * Nstate * Nstate + 1;
        for (i = 0; i < Ngrid; i++) {
            fprintf(den_file, "%18.8E", fi_time_grid[i] / unittrans_t);
            for (int j = 0; j < Nstate * Nstate; j++) {
                fprintf(den_file, "%18.8E", mpi_real_den[j * Ngrid + i]/Ntraj*mpi_size);
            }
            for (int j = 0; j < Nstate * Nstate; j++) {
                fprintf(den_file, "%18.8E", mpi_imag_den[j * Ngrid + i]/Ntraj*mpi_size);
            }
            fprintf(den_file, "\n");
        }
        fclose(den_file);
    }

    if (mpi_population != NULL) {
        strncpy(outname, filepath, len - 5);
        // strcpy(outname + len - 5, ".pop");
        outname[len - 5] = '\0'; // 确保字符串以null结尾
        strcpy(outname + len - 5, "_mpi");
        strcpy(outname + len - 5 + strlen("_mpi"), cid);
        strcpy(outname + len - 5 + strlen("_mpi") + strlen(cid), ".pop");
        FILE *pop_file = fopen(outname, "w");
        totn = Nstate + 1;
        for (i = 0; i < Ngrid; i++) {
            fprintf(pop_file, "%18.8E", fi_time_grid[i] / unittrans_t);
            for (int j = 0; j < Nstate; j++) {
                fprintf(pop_file, "%18.8E", mpi_population[j*Ngrid+i]/Ntraj*mpi_size);
            }
            fprintf(pop_file, "\n");
        }
        fclose(pop_file);
    }


}
