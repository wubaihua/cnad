

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

#ifdef x86
    void qm_msmodel(double *R, struct set_host *setm, struct set_slave *sets);
#endif

void initial_vari(struct set_slave *sets,struct set_host *seth) {
    int i;


    // printf("11111111111\n");
    // printf("N=%d %d %d\n",seth->Ndof1,seth->Ndof2,seth->Nstate);

    // 分配内存
    
   
    sets->R_nuc = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->P_nuc = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->mass = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->force = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));


    sets->R_nuc_init = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    // printf("22222\n");

    if (seth->forcetype == 1) {
        sets->force_nuc = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    }

    
    sets->V = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
    sets->dV = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));

    sets->R_nuc_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->P_nuc_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->force_old = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
    sets->V_old = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
    sets->dV_old = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
    
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
    } else if (seth->if_allcf == 1) {
        // sets->cf0 = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
        // memset(sets->cf0, 0, seth->Nstate * seth->Nstate * sizeof(double complex ));
        // // cfall = (double *)malloc(seth->Nstate * seth->Nstate * seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double ));
        // // memset(cfall, 0, seth->Nstate * seth->Nstate * seth->Nstate * seth->Nstate * seth->Ngrid * sizeof(double));
        // seth->type_evo = 1;
    } else if (seth->if_allcf >= 2) {
        sets->cf0 = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex ));
        memset(sets->cf0, 0, seth->Nstate * seth->Nstate * sizeof(double complex ));
        sets->cfeff = (double  complex *)malloc(seth->Ngrid * sizeof(double complex ));
        memset(sets->cfeff, 0, seth->Ngrid * sizeof(double complex ));
        seth->type_evo = 1;
        sets->weight0 = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
        sets->weightt = (double *)malloc(seth->Nstate * seth->Nstate * sizeof(double));
    }
    
    // if (seth->if_st_eng == 1) {
    //     sets->energy_est = (double *)malloc(seth->Ngrid * sizeof(double));
    // }
    //   printf("44444\n");
    
    // if (iflbd > 0) {
    //     rate_lbd = (double *)malloc(seth->Nstate * sizeof(double));
    //     memset(rate_lbd, 0, seth->Nstate * sizeof(double));
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    //     occ_lbd = sets->init_occ;
    

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

    
    //     }

    if (strcmp(seth->method, "FSSH") == 0 || strcmp(seth->method, "fssh") == 0 || 
        strcmp(seth->method, "FS-NAF") == 0 || strcmp(seth->method, "fs-naf") == 0){

        if(seth->if_default == 0) {
            seth->type_hop = 0;
            seth->direc_padj = 1;
            seth->ifreflp = 1;
        }

    } if (strcmp(seth->method, "GFSH") == 0 || strcmp(seth->method, "gfsh") == 0 ){

        if(seth->if_default == 0) {
            seth->type_hop = 2;
            seth->direc_padj = 1;
            seth->ifreflp = 1;
        }
           

    } else if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "MS-MASH") == 0 || strcmp(seth->method, "ms-mash") == 0 ||
               strcmp(seth->method, "msmash") == 0 || strcmp(seth->method, "MSMASH") == 0 ||
               strcmp(seth->method, "MASH-RM") == 0 || strcmp(seth->method, "mash-rm") == 0 ||
               strcmp(seth->method, "MA-NAF-MR") == 0 || strcmp(seth->method, "ma-naf-mr") == 0 ||
               strcmp(seth->method, "MA-NAF-RM") == 0 || strcmp(seth->method, "ma-naf-rm") == 0 ||
               strcmp(seth->method, "msmash2") == 0 || strcmp(seth->method, "MSMASH2") == 0 ||
               strcmp(seth->method, "msmash3") == 0 || strcmp(seth->method, "MSMASH3") == 0) {

        if(seth->if_default == 0) {       
            seth->type_hop = 1;
            seth->direc_padj = 0;
            seth->ifreflp = 0;
        }
    
    // } else if (strcmp(seth->method, "GDTWA") == 0 || strcmp(seth->method, "eGDTWA") == 0) {
    //     if (seth->type_evo != 3) seth->type_evo = 1;
    //     if (seth->if_inv_focus == 1) {
    //         sets->inverse_kernel = (double  complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double  complex ));
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
    }


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
    if (seth->rep == 1 || seth->sampletype == 3 || seth->rep == 2 || seth->rep == 3) {
        sets->U_d2a = (double  complex*)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
        sets->U_ref = (double  complex*)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
        sets->U_d2a_old = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
        sets->U_ref_old = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
        memset(sets->U_ref, 0, seth->Nstate * seth->Nstate * sizeof(double complex));
        sets->E_adia = (double *)malloc(seth->Nstate * sizeof(double));
        sets->dv_adia = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
        sets->nac = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
        sets->E_adia_old = (double *)malloc(seth->Nstate * sizeof(double));
        sets->dv_adia_old = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
        sets->nac_old = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
        // P_kin = (double *)malloc(seth->Ndof1 * seth->Ndof2 * sizeof(double));
        sets->nac_check = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
        sets->nac_check_old = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
        if (seth->type_prop_adia > 0) {
            sets->overlap_adia = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
        }
        memset(sets->nac_check, 0, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
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
    if (seth->if_Pdis == 1) {
        sets->expisp = (double  complex *)malloc(seth->s_N * sizeof(double complex ));
        memset(sets->expisp, 0, seth->s_N * sizeof(double complex ));
    }
   
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

    if (seth->mean_nuc == 1) {
        sets->R_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        sets->P_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        sets->R2_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        sets->P2_nuc_mean = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        memset(sets->R_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        memset(sets->P_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        memset(sets->R2_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        memset(sets->P2_nuc_mean, 0, seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        // if (ifmsbranch > 0) {
        //     sets->R_nuc_oldtraj = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        //     sets->P_nuc_oldtraj = (double *)malloc(seth->Ndof1 * seth->Ndof2 * seth->Ngrid * sizeof(double));
        // }
    }

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
    if (seth->ifswitchforce > 0 || seth->type_hop == 1) {
        if (seth->rep == 0) {
            sets->U_d2a = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
            sets->E_adia = (double *)malloc(seth->Nstate * sizeof(double));
            sets->U_d2a_old = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
            sets->E_adia_old = (double *)malloc(seth->Nstate * sizeof(double));
            sets->U_ref = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
            sets->U_ref_old = (double complex *)malloc(seth->Nstate * seth->Nstate * sizeof(double complex));
            memset(sets->U_ref, 0, seth->Nstate * seth->Nstate * sizeof(double complex));
            if (seth->direc_padj == 0 || seth->direc_padj == 1){
                sets->dv_adia = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
                sets->nac = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
                sets->dv_adia_old = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
                sets->nac_old = (double complex *)malloc(seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));
            }
        
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
    double complex tempcv1[seth->Nstate], tempcv2[seth->Nstate];

    if (seth->if_restart == 1) {
        return;
    }

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
               strcmp(seth->method, "CMM") == 0 || strcmp(seth->method, "cmm") == 0 ||
               strcmp(seth->method, "scmm") == 0 || strcmp(seth->method, "SCMM") == 0 ||
               strcmp(seth->method, "wMM") == 0 || strcmp(seth->method, "wmm") == 0 ) {
        
        // 调用 random_prob 函数
        random_prob(seth->Nstate, action);
        // sets->init_occ=1;
        // for (int i = 0; i < seth->Nstate; i++) {
        //     action[i] = (1.0)/(seth->Nstate);
        //     theta[i] = 1.5;
        //     // printf("action[%d]=%f\n",i,action[i]);
        // }
        // action[sets->init_occ - 1] = 0.5;
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
        if (seth->if_allcf != 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->cf0[i * seth->Nstate + j] = seth->Nstate * (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
                }
            }
        }
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

        // 计算 sets->cf0
        if (seth->if_allcf != 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->cf0[i * seth->Nstate + j] = seth->Nstate * (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
                }
            }
        }
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

        // debug
        // action[0]=1.0279163790587629,action[1]=0.6866200133888645;
        // theta[0]=4.435113217752501,theta[1]=1.0379228857653318;
        // sets->P_nuc[0] = 7.855740656524664;
        // printf("%d %18.8e  %18.8e  %18.8e  %18.8e \n",seth->mpi_rank,action[0],action[1],theta[0],theta[1]);
      

        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }

     
        
        if(seth->if_scale_sqc == 0)seth->gamma_zpe = 1.0 / 3.0 ;

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


        // 计算 sets->cf0
        if (seth->if_allcf != 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    if (i == j) {
                        if (i == sets->init_occ - 1){
                            sets->cf0[i * seth->Nstate + i] = seth->Nstate; 
                        } else {
                            sets->cf0[i * seth->Nstate + i] = 0.0; 
                        }
                    } else {
                        if (i == sets->init_occ - 1 || j == sets->init_occ - 1){
                            sets->cf0[i * seth->Nstate + j] = seth->Nstate * 0.6 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                        } else {
                            sets->cf0[i * seth->Nstate + j] = 0.0;
                        }
                    }
                }
            }
        }

        if (seth->if_scale_sqc == 1){
            sets->scale_sqc2 = 0;
            for (i = 0; i < seth->Nstate; i++){
                sets->scale_sqc2 += (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5;
            }
            sets->scale_sqc2 = (1.0 + seth->Nstate * seth->gamma_zpe) / sets->scale_sqc2;
            for (i = 0; i < seth->Nstate; i++){
                sets->xe[i] *= sqrt(sets->scale_sqc2);
                sets->pe[i] *= sqrt(sets->scale_sqc2);
            }
            if (seth->type_evo >= 1){
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    sets->den_e[i] *= sqrt(sets->scale_sqc2);
                }
            }
        }


        
        

    } else if (strcmp(seth->method, "fssh") == 0 || strcmp(seth->method, "FSSH") == 0 ||
               strcmp(seth->method, "fs-naf") == 0 || strcmp(seth->method, "FS-NAF") == 0 ||
        strcmp(seth->method, "GFSH") == 0 || strcmp(seth->method, "gfsh") == 0 ) {
        for (i = 0; i < seth->Nstate; i++) {
            action[i] = 0;
        }
        action[sets->init_occ-1] = 1;
    
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
            
        }
        
        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        sets->correfun_0 = 1.0;
        
        sets->id_state = sets->init_occ-1;

        if (seth->type_evo >= 1) {
            for (i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->den_e[i] = 0.0 + 0.0 * I;
            }
            sets->den_e[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] = 1.0 + 0.0 * I;
        }

    } else if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "ma-naf-mr") == 0 || strcmp(seth->method, "MA-NAF-MR") == 0 ) {
        
        random_prob(seth->Nstate, action);
        
        if(seth->if_flighttime_tully == 1) {
            sets->id_state = maxloc(action, seth->Nstate);
                while (sets->id_state != sets->init_occ - 1) {
                random_prob(seth->Nstate, action);
                sets->id_state = maxloc(action, seth->Nstate);
            }
        }


        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }
        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        sets->measure_mash = fabs(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0]
                                 -sets->xe[1] * sets->xe[1] - sets->pe[1] * sets->pe[1]);

        memset(sets->rho0_mash,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        if(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0] >= 1){
            sets->rho0_mash[0] = 1;
            sets->id_state = 0;
            if(seth->if_flighttime_tully == 1) sets->rho0_mash[0] = 0.5;
        } else {
            sets->rho0_mash[3] = 1;
            sets->id_state = 1;
            if(seth->if_flighttime_tully == 1) sets->rho0_mash[3] = 0.5;
        }
        sets->rho0_mash[1] = 0.5 * (sets->xe[0] + I * sets->pe[0]) * (sets->xe[1] - I * sets->pe[1]);
        sets->rho0_mash[2] = 0.5 * (sets->xe[0] - I * sets->pe[0]) * (sets->xe[1] + I * sets->pe[1]);

        sets->correfun_0 = 1.0;


    } else if (strcmp(seth->method, "MASH-MF") == 0 || strcmp(seth->method, "mash-mf") == 0 ||
               strcmp(seth->method, "mashmf") == 0 || strcmp(seth->method, "MASHMF") == 0 ||
               strcmp(seth->method, "mr") == 0 || strcmp(seth->method, "MR") == 0) {
        random_prob(seth->Nstate, action);
        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }
        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        sets->measure_mash = fabs(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0]
                                 -sets->xe[1] * sets->xe[1] - sets->pe[1] * sets->pe[1]);

        memset(sets->rho0_mash,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        if(sets->xe[sets->init_occ - 1] * sets->xe[sets->init_occ - 1] + sets->pe[sets->init_occ - 1] * sets->pe[sets->init_occ - 1] >= 1){
            sets->correfun_0 = 1.0;
        } else {
            sets->correfun_0 = 0.0;
        }
        
    
    } else if (strcmp(seth->method, "MS-MASH") == 0 || strcmp(seth->method, "ms-mash") == 0 ||
               strcmp(seth->method, "msmash") == 0 || strcmp(seth->method, "MSMASH") == 0 ||
               strcmp(seth->method, "MASH-RM") == 0 || strcmp(seth->method, "mash-rm") == 0 ||
               strcmp(seth->method, "MS-MASH-MF") == 0 || strcmp(seth->method, "ms-mash-mf") == 0 ||
               strcmp(seth->method, "msmashmf") == 0 || strcmp(seth->method, "MSMASHMF") == 0 ||
               strcmp(seth->method, "RM") == 0 || strcmp(seth->method, "rm") == 0 ||
               strcmp(seth->method, "MA-NAF-RM") == 0 || strcmp(seth->method, "ma-naf-rm") == 0 ) {
        
        sets->id_state=-10;
        while (sets->id_state != sets->init_occ - 1) {
            random_prob(seth->Nstate, action);
            sets->id_state = maxloc(action, seth->Nstate);
        }
        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }

        seth->gamma_zpe = 0.0;

        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        // for (int i = 0; i < seth->Nstate; i++) {
        //     for (int j = 0; j < seth->Nstate; j++) {
        //         sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
        //     }
        // }

        sets->correfun_0 = 1.0;


        if (seth->type_evo >= 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }


        if (seth->if_allcf != 0) {
            alpha_mash = 0;
            for (i = 1; i < seth->Nstate + 1; i++){
                alpha_mash += 1.0 / ((double) i); 
            }
            alpha_mash = (double) (seth->Nstate - 1.0) / (alpha_mash - 1.0);
            beta_mash = (1.0 - alpha_mash) / seth->Nstate;
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    if (i == j) {
                        if( i == sets->init_occ - 1 ){
                           sets->cf0[i * seth->Nstate + j] = seth->Nstate;
                        } else {
                           sets->cf0[i * seth->Nstate + j] = 0;
                        }
                    } else {
                        sets->cf0[i * seth->Nstate + j] = seth->Nstate * (1.0 + seth->Nstate)/((1 - seth->Nstate * beta_mash)) * sets->den_e[i * seth->Nstate + j];
                    }
                }
            }
        }

    } else if (strcmp(seth->method, "msmash2") == 0 || strcmp(seth->method, "MSMASH2") == 0 ) {
        
        sets->id_state=-10;
        while (sets->id_state != sets->init_occ - 1) {
            random_prob(seth->Nstate, action);
            sets->id_state = maxloc(action, seth->Nstate);
        }
        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        }

        seth->gamma_zpe = 0.0;

        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        // for (int i = 0; i < seth->Nstate; i++) {
        //     for (int j = 0; j < seth->Nstate; j++) {
        //         sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
        //     }
        // }

        alpha_mash = 0;
        x2 = 0;
        for (i = 1; i < seth->Nstate + 1; i++){
            alpha_mash += 1.0 / ((double) i); 
            x2 += 1.0 / ((double) i * i); 
        }
        beta_mash = ((seth->Nstate + 1) * alpha_mash - alpha_mash * alpha_mash - x2) / ((seth->Nstate + 1) * seth->Nstate * seth->Nstate * (seth->Nstate - 1)); // \Gamma
        alpha_mash = alpha_mash / (seth->Nstate * seth->Nstate ); // I = (\sum_{n=1}^F 1/n) /F^2

        sets->correfun_0 = 1.0 / (2 * alpha_mash) / seth->Nstate;


        if (seth->type_evo >= 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }


        if (seth->if_allcf != 0) {
            alpha_mash = 0;
            x2 = 0;
            for (i = 1; i < seth->Nstate + 1; i++){
                alpha_mash += 1.0 / ((double) i); 
                x2 += 1.0 / ((double) i * i); 
            }
            beta_mash = ((seth->Nstate + 1) * alpha_mash - alpha_mash * alpha_mash - x2) / ((seth->Nstate + 1) * seth->Nstate * seth->Nstate * (seth->Nstate - 1)); // \Gamma
            alpha_mash = alpha_mash / (seth->Nstate * seth->Nstate ); // I = (\sum_{n=1}^F 1/n) /F^2
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    if (i == j) {
                        if( i == sets->init_occ - 1 ){
                           sets->cf0[i * seth->Nstate + j] = 1.0 / (2 * alpha_mash);
                        } else {
                           sets->cf0[i * seth->Nstate + j] = 0;
                        }
                    } else {
                        if (i == sets->id_state || j == sets->id_state){
                            sets->cf0[i * seth->Nstate + j] = 1.0 / (2 *beta_mash) * sets->den_e[i * seth->Nstate + j];
                        } else {
                            sets->cf0[i * seth->Nstate + j] = 0;
                        }
                    }
                }
            }
            
           

        }


    } else if (strcmp(seth->method, "msmash3") == 0 || strcmp(seth->method, "MSMASH3") == 0 ) {
        
        // sets->id_state=-10;
        // while (sets->id_state != sets->init_occ - 1) {
        //     random_prob(seth->Nstate, action);
        //     sets->id_state = maxloc(action, seth->Nstate);
        // }
        // for (int i = 0; i < seth->Nstate; i++) {
        //     sets->xe[i] = sqrt(2 * action[i]) * cos(theta[i]);
        //     sets->pe[i] = sqrt(2 * action[i]) * sin(theta[i]);
        // }

        // seth->gamma_zpe = 0.0;

        // memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        // // for (int i = 0; i < seth->Nstate; i++) {
        // //     for (int j = 0; j < seth->Nstate; j++) {
        // //         sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
        // //     }
        // // }

        // alpha_mash = 0;
        // x2 = 0;
        // for (i = 1; i < seth->Nstate + 1; i++){
        //     alpha_mash += 1.0 / ((double) i); 
        //     x2 += 1.0 / ((double) i * i); 
        // }
        // beta_mash = ((seth->Nstate + 1) * alpha_mash - alpha_mash * alpha_mash - x2) / ((seth->Nstate + 1) * seth->Nstate * seth->Nstate * (seth->Nstate - 1)); // \Gamma
        // alpha_mash = alpha_mash / (seth->Nstate * seth->Nstate ); // I = (\sum_{n=1}^F 1/n) /F^2

        // sets->correfun_0 = 1.0 / (2 * alpha_mash) / seth->Nstate;


        // if (seth->type_evo >= 1) {
        //     for (int i = 0; i < seth->Nstate; i++) {
        //         for (int j = 0; j < seth->Nstate; j++) {
        //             sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
        //         }
        //     }
        // }


        // if (seth->if_allcf != 0) {
        //     if (seth->rep == 0) {
        //         V_msmodel(sets->R_nuc, sets->V, 0.0, seth);
                
        //         memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
        //         transpose(sets->U_d2a,tempdm1,seth->Nstate);
        //         dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        //         cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                
        //         memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
        //         dc_matmul(tempdm1,gamma_cv_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        //         cd_matmul(tempcm1,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        //     }

        //     for (int i = 0; i < seth->Nstate; i++) {

        //         c_main[i] = creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]);

        //     }

        //     if (seth->ifscalegamma == 1) {
        //         for (int i = 0; i < seth->Nstate; i++) {

        //             c_main[i] = creal(sets->den_e[i * seth->Nstate + i] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + i]);

        //         }
        //     }

        //     sets->id_state = maxloc(c_main, seth->Nstate);



        //     alpha_mash = 0;
        //     x2 = 0;
        //     for (i = 1; i < seth->Nstate + 1; i++){
        //         alpha_mash += 1.0 / ((double) i); 
        //         x2 += 1.0 / ((double) i * i); 
        //     }
        //     beta_mash = ((seth->Nstate + 1) * alpha_mash - alpha_mash * alpha_mash - x2) / ((seth->Nstate + 1) * seth->Nstate * seth->Nstate * (seth->Nstate - 1)); // \Gamma
        //     alpha_mash = alpha_mash / (seth->Nstate * seth->Nstate ); // I = (\sum_{n=1}^F 1/n) /F^2
        //     for (int i = 0; i < seth->Nstate; i++) {
        //         for (int j = 0; j < seth->Nstate; j++) {
        //             if (i == j) {
        //                 if( i == sets->id_state){
        //                    sets->cf0[i * seth->Nstate + j] = 1.0 / (2 * alpha_mash);
        //                 } else {
        //                    sets->cf0[i * seth->Nstate + j] = 0;
        //                 }
        //             } else {
        //                 if (i == sets->id_state || j == sets->id_state){
        //                     sets->cf0[i * seth->Nstate + j] = 1.0 / (2 *beta_mash) * sets->den_e[i * seth->Nstate + j];
        //                 } else {
        //                     sets->cf0[i * seth->Nstate + j] = 0;
        //                 }
        //             }
        //         }
        //     }

        //     transpose(sets->U_d2a,tempdm1,seth->Nstate);
        //     dc_matmul(sets->U_d2a,sets->cf0,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        //     cd_matmul(tempcm1,tempdm1,sets->cf0,seth->Nstate,seth->Nstate,seth->Nstate);


        //     if (seth->rep == 0) {
        //         if (seth->type_evo == 0) {
        //             memcpy(sets->xe,xe_save,seth->Nstate*sizeof(double));
        //             memcpy(sets->pe,pe_save,seth->Nstate*sizeof(double));
        //         } else if (seth->type_evo == 1) {
        //             memcpy(sets->den_e,den_e_save,seth->Nstate * seth->Nstate * sizeof(double complex));
        //         }
        //         memcpy(sets->gamma_cv,gamma_cv_save,seth->Nstate * seth->Nstate * sizeof(double complex));
        //     }   

        // }

    } else if (strcmp(seth->method, "ms-mash-mf2") == 0 || strcmp(seth->method, "MS-MASH-MF2") == 0 ||
               strcmp(seth->method, "mf2") == 0 || strcmp(seth->method, "MF2") == 0 ||
               strcmp(seth->method, "CW1") == 0 || strcmp(seth->method, "cw1") == 0 ) {
        
        
        random_prob(seth->Nstate, action);
        alpha_mash = 0;
        for (i = 1; i < seth->Nstate + 1; i++){
            alpha_mash += 1.0 / ((double) i); 
        }
        
        alpha_mash = (double) (seth->Nstate - 1.0) / (alpha_mash - 1.0);
        beta_mash = (1.0 - alpha_mash) / seth->Nstate;
        seth->gamma_zpe = -1.0 * beta_mash;


        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * (1.0 + seth->Nstate * seth->gamma_zpe) * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * (1.0 + seth->Nstate * seth->gamma_zpe) * action[i]) * sin(theta[i]);
        }

        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        for (int i = 0; i < seth->Nstate; i++) {
            sets->gamma_cv[i * seth->Nstate + i] = seth->gamma_zpe;
        }

        sets->correfun_0 = seth->Nstate * (0.5 * (sets->xe[sets->init_occ-1] * sets->xe[sets->init_occ-1] + sets->pe[sets->init_occ-1] * sets->pe[sets->init_occ-1]) - seth->gamma_zpe);


        if (seth->type_evo >= 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

        if (seth->if_allcf != 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->cf0[i * seth->Nstate + j] = seth->Nstate * (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
                }
            }
        }

    } else if (strcmp(seth->method, "mf3") == 0 || strcmp(seth->method, "MF3") == 0 ||
               strcmp(seth->method, "CW2") == 0 || strcmp(seth->method, "cw2") == 0 ) {
        
        
        random_prob(seth->Nstate, action);
        
        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * (1.0 + seth->Nstate * seth->gamma_zpe) * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * (1.0 + seth->Nstate * seth->gamma_zpe) * action[i]) * sin(theta[i]);
        }

        memset(sets->gamma_cv,0,seth->Nstate*seth->Nstate*sizeof(double complex));
        for (int i = 0; i < seth->Nstate; i++) {
            sets->gamma_cv[i * seth->Nstate + i] = seth->gamma_zpe;
        }

        sets->correfun_0 = seth->Nstate * (0.5 * (sets->xe[sets->init_occ-1] * sets->xe[sets->init_occ-1] + sets->pe[sets->init_occ-1] * sets->pe[sets->init_occ-1]) - seth->gamma_zpe);

        if (seth->type_evo >= 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

        if (seth->if_allcf != 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->cf0[i * seth->Nstate + j] = seth->Nstate * (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
                }
            }
        }

    }  else if (strcmp(seth->method, "NW") == 0 || strcmp(seth->method, "nw") == 0 ) {
        
        sets->id_state=-10;

        while (sets->id_state != sets->init_occ - 1) {
            random_prob(seth->Nstate, action);
            sets->id_state = maxloc(action, seth->Nstate);
        }
        for (int i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = sqrt(2 * (1 + seth->Nstate * seth->gamma_zpe) * action[i]) * cos(theta[i]);
            sets->pe[i] = sqrt(2 * (1 + seth->Nstate * seth->gamma_zpe) * action[i]) * sin(theta[i]);
        }

        
        for (int i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                sets->gamma_cv[i * seth->Nstate + j] = (i == j) ? seth->gamma_zpe : 0.0;
            }
        }

        if (seth->type_evo >= 1) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->den_e[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

        sets->correfun_0 = 1.0;


        if (seth->if_allcf != 0) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    if (i == j) {
                        if (i == sets->init_occ - 1){
                            sets->cf0[i * seth->Nstate + i] = seth->Nstate; 
                        } else {
                            sets->cf0[i * seth->Nstate + i] = 0.0; 
                        }
                    } else {
                        sets->cf0[i * seth->Nstate + j] = seth->Nstate * 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                    }
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
                if (i == j || i == (sets->init_occ - 1) || j == (sets->init_occ -1)) continue;
                sets->gamma_cv[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                if (seth->ifscalegamma == 1) {
                    sets->gamma_cv[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe);
                }
            }
        }
    }



    sets->if_ad_nac = 0;
    if (seth->sampletype == 2){
        
        dV_msmodel(sets->R_nuc, sets->dV,seth);
        V_msmodel(sets->R_nuc, sets->V, 0.0,seth);
        cal_NACV(sets,seth);

        // sets->xe = matmul(transpose(sets->U_d2a), sets->xe)
        transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
        // memcpy(tempdv1, sets->xe, seth->Nstate * sizeof(double));
        // dd_matmul(tempdm1, tempdv1, sets->xe, seth->Nstate, seth->Nstate, 1);
        // memcpy(tempdv1, sets->pe, seth->Nstate * sizeof(double));
        // dd_matmul(tempdm1, tempdv1, sets->pe, seth->Nstate, seth->Nstate, 1);
        for (i = 0; i < seth->Nstate; i++) {
            tempcv1[i] = sets->xe[i] + I * sets->pe[i];
        }
        cc_matmul(tempcm1, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = creal(tempcv2[i]);
            sets->pe[i] = cimag(tempcv2[i]);
        }

        // sets->pe = matmul(transpose(sets->U_d2a), sets->pe)
        cc_matmul(tempcm1,sets->gamma_cv,tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
        cc_matmul(tempcm2, sets->U_d2a, sets->gamma_cv, seth->Nstate, seth->Nstate, seth->Nstate);
        // sets->gamma_cv = matmul(matmul(transpose(sets->U_d2a), sets->gamma_cv), sets->U_d2a)
        
        // 这里需要实现矩阵乘法和转置操作

        if (seth->type_evo == 1 || seth->type_evo == 3) {
            // sets->den_e = matmul(matmul(transpose(sets->U_d2a), sets->den_e), sets->U_d2a)
            cc_matmul(tempcm1,sets->den_e,tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
            cc_matmul(tempcm2, sets->U_d2a, sets->den_e, seth->Nstate, seth->Nstate, seth->Nstate);
        }
        // if (sets->den_e4nuc != NULL) {
        //     // sets->den_e4nuc = matmul(matmul(transpose(sets->U_d2a), sets->den_e4nuc), sets->U_d2a)
        // }
        // if (seth->type_evo == 4) {
        //     // G_xpconfg = matmul(transpose(sets->U_d2a), G_xpconfg)
        // }

        if (strcmp(seth->method, "fssh") == 0 || strcmp(seth->method, "FSSH") == 0 || strcmp(seth->method, "SC-FSSH") == 0 || strcmp(seth->method, "sc-fssh") == 0 ||
            strcmp(seth->method, "CC-FSSH") == 0 || strcmp(seth->method, "cc-fssh") == 0 || strcmp(seth->method, "fsshswitch") == 0 || strcmp(seth->method, "pcsh") == 0 ||
            strcmp(seth->method, "PCSH") == 0 || strcmp(seth->method, "PCSH-NAF") == 0 || strcmp(seth->method, "pcsh-naf") == 0 || strcmp(seth->method, "BCSH") == 0 ||
            strcmp(seth->method, "bcsh") == 0 || strcmp(seth->method, "BCSH-NAF") == 0 || strcmp(seth->method, "bcsh-naf") == 0 ||
            strcmp(seth->method, "fs-naf") == 0 || strcmp(seth->method, "FS-NAF") == 0 ||
            strcmp(seth->method, "GFSH") == 0 || strcmp(seth->method, "gfsh") == 0 ) {

            double x2 = (double)rand() / RAND_MAX;
            double ps1, ps2;
            for (int i = 0; i < seth->Nstate; i++) {
                sets->id_state = i;
                
                if (i == 0) {
                    ps1 = 0;
                    // ps2 = sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + 0] * sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + 0];
                    ps2 = cabs(sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + 0]) * cabs(sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + 0]);
                } else {
                    ps1 = 0;
                    for (int k = 0; k < i - 1; k++) {
                        // ps1 += sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + k] * sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + k];
                        ps1 += cabs(sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + k]) * cabs(sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + k]);
                    }
                    // ps2 = ps1 + sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + i] * sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + i];
                    ps2 = ps1 + cabs(sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + i]) * cabs(sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + i]);
                }
                if (x2 >= ps1 && x2 < ps2) {
                    break;
                }
            }
        
        } else if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "ma-naf-mr") == 0 || strcmp(seth->method, "MA-NAF-MR") == 0) {
            
            memset(sets->rho0_mash,0,seth->Nstate*seth->Nstate*sizeof(double complex));
            if(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0] >= 1){
                sets->rho0_mash[0] = 1;
                sets->id_state = 0;
            } else {
                sets->rho0_mash[3] = 1;
                sets->id_state = 1;
            }
            sets->rho0_mash[1] = 0.5 * (sets->xe[0] + I * sets->pe[0]) * (sets->xe[1] - I * sets->pe[1]);
            sets->rho0_mash[2] = 0.5 * (sets->xe[0] - I * sets->pe[0]) * (sets->xe[1] + I * sets->pe[1]);

            sets->measure_mash = fabs(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0]
                                 -sets->xe[1] * sets->xe[1] - sets->pe[1] * sets->pe[1]);

            memcpy(sets->U0_mash,sets->U_d2a,seth->Nstate * seth->Nstate * sizeof(double complex));
            
        } else if (strcmp(seth->method, "MS-MASH") == 0 || strcmp(seth->method, "ms-mash") == 0 ||
               strcmp(seth->method, "msmash") == 0 || strcmp(seth->method, "MSMASH") == 0 ||
               strcmp(seth->method, "MASH-RM") == 0 || strcmp(seth->method, "mash-rm") == 0 ||
               strcmp(seth->method, "MA-NAF-RM") == 0 || strcmp(seth->method, "ma-naf-rm") == 0 ||
               strcmp(seth->method, "msmash2") == 0 || strcmp(seth->method, "MSMASH2") == 0 ||
               strcmp(seth->method, "msmash3") == 0 || strcmp(seth->method, "MSMASH3") == 0) {
            
            for (int i = 0; i < seth->Nstate; i++) {
                if(seth->type_evo == 0) c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
                if(seth->type_evo == 1) c_main[i] = creal(sets->den_e[i * seth->Nstate + i]);
            }
            sets->id_state = maxloc(c_main, seth->Nstate);
    
        }
    }

   
    if (seth->ifswitchforce > 0 && seth->type_hop == 1 || seth->type_hop == 1) {
        if (seth->rep == 0 || seth->rep == 3) {
            
            if (strcmp(seth->msmodelname, "mole") == 0 ) {
                #ifdef x86
                qm_msmodel(sets->R_nuc, seth, sets); 
                #endif  
            } else {
                V_msmodel(sets->R_nuc, sets->V, 0.0, seth);
            }
            #ifdef sunway
            int slavecore_id = athread_get_id(-1);
            #endif
          
            
            dia_hermitemat(seth->Nstate, sets->V, sets->E_adia, sets->U_d2a);

            // printf("%d\n",slavecore_id);

            if (seth->type_evo == 0) {
                memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
                transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
                for (i = 0; i < seth->Nstate; i++) {
                    tempcv1[i] = sets->xe[i] + I * sets->pe[i];
                }
                cc_matmul(tempcm1, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
                for (i = 0; i < seth->Nstate; i++) {
                    sets->xe[i] = creal(tempcv2[i]);
                    sets->pe[i] = cimag(tempcv2[i]);
                }
            } else if (seth->type_evo == 1) {
                memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
                cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            }
            memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
            cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
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

        if (seth->rep == 0 || seth->rep == 3) {
            if (seth->type_evo == 0) {
                memcpy(sets->xe,xe_save,seth->Nstate*sizeof(double));
                memcpy(sets->pe,pe_save,seth->Nstate*sizeof(double));
            } else if (seth->type_evo == 1) {
                memcpy(sets->den_e,den_e_save,seth->Nstate * seth->Nstate * sizeof(double complex));
            }
            memcpy(sets->gamma_cv,gamma_cv_save,seth->Nstate * seth->Nstate * sizeof(double complex));
        }

    } else if (seth->ifswitchforce > 0 && seth->type_hop != 1) {
        if (seth->rep == 0 || seth->rep == 3) {
            
            if (strcmp(seth->msmodelname, "mole") == 0 ) {
                #ifdef x86
                qm_msmodel(sets->R_nuc, seth, sets);  
                #endif 
            } else {
                V_msmodel(sets->R_nuc, sets->V, 0.0, seth);
            }
            #ifdef sunway
            int slavecore_id = athread_get_id(-1);
            #endif
          
            
            dia_hermitemat(seth->Nstate, sets->V, sets->E_adia, sets->U_d2a);

            // printf("%d\n",slavecore_id);

            if (seth->type_evo == 0) {
                memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
                transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
                for (i = 0; i < seth->Nstate; i++) {
                    tempcv1[i] = sets->xe[i] + I * sets->pe[i];
                }
                cc_matmul(tempcm1, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
                for (i = 0; i < seth->Nstate; i++) {
                    sets->xe[i] = creal(tempcv2[i]);
                    sets->pe[i] = cimag(tempcv2[i]);
                }
            } else if (seth->type_evo == 1) {
                memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
                cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            }
            memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
            cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        }

        if(seth->type_hop == 3 || seth->type_hop == 4) {
            for (int i = 0; i < seth->Nstate; i++) {
                if (seth->type_evo == 0) {
                    c_main[i] = fabs((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) / 2 - creal(sets->gamma_cv[i * seth->Nstate + i]));
                } else if (seth->type_evo == 1) {
                    c_main[i] = fabs(creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]));
                }
            }
            double sum1 = 0.0;
            for (int i = 0; i < seth->Nstate; i++) {
                sum1 += c_main[i];
            }
            for (int i = 0; i < seth->Nstate; i++) {
                c_main[i] = c_main[i] / sum1;
            }

        }


        double r1 = (double)rand() / RAND_MAX;
        double ps1, ps2;
        for (int i = 0; i < seth->Nstate; i++) {
            sets->id_state = i;
            
            if (i == 0) {
                ps1 = 0;
                // ps2 = sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + 0] * sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + 0];
                ps2 = c_main[0];
            } else {
                ps1 = 0;
                for (int k = 0; k < i - 1; k++) {
                    // ps1 += sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + k] * sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + k];
                    ps1 += c_main[k];
                }
                // ps2 = ps1 + sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + i] * sets->U_d2a[(sets->init_occ - 1) * seth->Nstate + i];
                ps2 = ps1 + c_main[i];
            }
            if (r1 >= ps1 && r1 < ps2) {
                break;
            }
        }

    

        if (seth->rep == 0 || seth->rep == 3) {
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
    double tempdm[seth->Nstate*seth->Nstate], tempdm1[seth->Nstate * seth->Nstate],tempdm2[seth->Nstate*seth->Nstate];
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate*seth->Nstate], tempcm3[seth->Nstate*seth->Nstate];
    double c_main[seth->Nstate];
     double xe_save[seth->Nstate], pe_save[seth->Nstate];
    double complex gamma_cv_save[seth->Nstate * seth->Nstate], den_e_save[seth->Nstate * seth->Nstate];
    double complex tempcv1[seth->Nstate], tempcv2[seth->Nstate];
    #ifdef sunway
    int slavecore_id = athread_get_id(-1);
    #endif
    
    
    // tempmat1[seth->Nstate*seth->Nstate],tempmat2[seth->Nstate*seth->Nstate],tempmat3[seth->Nstate*seth->Nstate]



    if (seth->sampletype == 2) {
        
        // memcpy(tempv,sets->xe,seth->Nstate * sizeof(double));
        // dd_matmul(sets->U_d2a,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
        // memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
        // dd_matmul(sets->U_d2a,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
        for (i = 0; i < seth->Nstate; i++) {
            tempcv1[i] = sets->xe[i] + I * sets->pe[i];
        }
        cc_matmul(sets->U_d2a, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = creal(tempcv2[i]);
            sets->pe[i] = cimag(tempcv2[i]);
        }
        
        
        memcpy(tempcm,sets->gamma_cv,seth->Nstate * seth->Nstate*sizeof(double complex));
        
        // transpose(sets->U_d2a,tempdm,seth->Nstate);
        transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
        cc_matmul(sets->gamma_cv,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
        cc_matmul(sets->U_d2a,tempcm2,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        

        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            cc_matmul(sets->den_e,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(sets->U_d2a,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
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

    } else if (strcmp(seth->method, "SCMM") == 0 || strcmp(seth->method, "scmm") == 0 ||
        strcmp(seth->method, "wMM") == 0 || strcmp(seth->method, "wmm") == 0) {
        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
                sets->correfun_t[i] = sets->den_e[i];
            }
            for (int i = 0; i < seth->Nstate; i++) {
                sets->correfun_t[i * seth->Nstate + i] -= seth->gamma_zpe;
            }
        } else {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    sets->correfun_t[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - seth->gamma_zpe* (i == j ? 1 : 0);
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

        if (seth->if_scale_sqc == 1){
            for (i = 0; i < seth->Nstate; i++){
                sets->xe[i] /= sqrt(sets->scale_sqc2);
                sets->pe[i] /= sqrt(sets->scale_sqc2);
            }
            if (seth->type_evo >= 1){
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    sets->den_e[i] /= sqrt(sets->scale_sqc2);
                }
            }
        }

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

        if (seth->if_scale_sqc == 1){
            for (i = 0; i < seth->Nstate; i++){
                sets->xe[i] *= sqrt(sets->scale_sqc2);
                sets->pe[i] *= sqrt(sets->scale_sqc2);
            }
            if (seth->type_evo >= 1){
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    sets->den_e[i] *= sqrt(sets->scale_sqc2);
                }
            }
        }

    } else if (strcmp(seth->method, "fssh") == 0 || strcmp(seth->method, "FSSH") == 0 ||
                strcmp(seth->method, "fs-naf") == 0 || strcmp(seth->method, "FS-NAF") == 0 ||
                strcmp(seth->method, "GFSH") == 0 || strcmp(seth->method, "gfsh") == 0 ) {
        if(seth->sampletype == 2){
            // transpose(sets->U_d2a,tempdm,seth->Nstate);
            // memcpy(tempv,sets->xe,seth->Nstate*sizeof(double));
            // dd_matmul(tempdm,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
            // memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
            // dd_matmul(tempdm,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
            transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
            for (i = 0; i < seth->Nstate; i++) {
                tempcv1[i] = sets->xe[i] + I * sets->pe[i];
            }
            cc_matmul(tempcm1, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
            for (i = 0; i < seth->Nstate; i++) {
                sets->xe[i] = creal(tempcv2[i]);
                sets->pe[i] = cimag(tempcv2[i]);
            }
            if (seth->type_evo == 1 || seth->type_evo == 3) {
                cc_matmul(tempcm1, sets->den_e,  tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
                cc_matmul(tempcm2, sets->U_d2a, sets->den_e, seth->Nstate, seth->Nstate, seth->Nstate);
            }
        }
        for (i = 0; i < seth->Nstate; i++) {
            for (j = 0; j < seth->Nstate; j++) {
                if( i == j && i == sets->id_state){
                    sets->correfun_t[i * seth->Nstate + i] = 1.0;
                }
                else if ( i == j && i != sets->id_state){
                    sets->correfun_t[i * seth->Nstate + i] = 0.0;
                }
                else {
                    sets->correfun_t[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                    if (seth->type_evo == 1 || seth->type_evo == 3) {
                        sets->correfun_t[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j];
                    }
                }
            }
        }


        if(seth->sampletype == 2){
            // transpose(sets->U_d2a,tempdm,seth->Nstate);
            // cd_matmul(sets->correfun_t,tempdm,tempcm,seth->Nstate,seth->Nstate,seth->Nstate);
            // dc_matmul(sets->U_d2a,tempcm,sets->correfun_t,seth->Nstate,seth->Nstate,seth->Nstate);
            transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
            cc_matmul(sets->correfun_t, tempcm1, tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
            cc_matmul(sets->U_d2a, tempcm2, sets->correfun_t, seth->Nstate, seth->Nstate, seth->Nstate);
            // memcpy(tempv,sets->xe,seth->Nstate * sizeof(double));
            // dd_matmul(sets->U_d2a,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
            // memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
            // dd_matmul(sets->U_d2a,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
            for (i = 0; i < seth->Nstate; i++) {
                tempcv1[i] = sets->xe[i] + I * sets->pe[i];
            }
            cc_matmul(sets->U_d2a, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
            for (i = 0; i < seth->Nstate; i++) {
                sets->xe[i] = creal(tempcv2[i]);
                sets->pe[i] = cimag(tempcv2[i]);
            }
            if (seth->type_evo == 1 || seth->type_evo == 3) {
                cc_matmul(sets->den_e,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(sets->U_d2a,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            }
        }
    } else if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "ma-naf-mr") == 0 || strcmp(seth->method, "MA-NAF-MR") == 0) {
        
        if(seth->sampletype == 2){
            // transpose(sets->U_d2a,tempdm,seth->Nstate);
            // memcpy(tempv,sets->xe,seth->Nstate*sizeof(double));
            // dd_matmul(tempdm,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
            // memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
            // dd_matmul(tempdm,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
            transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
            for (i = 0; i < seth->Nstate; i++) {
                tempcv1[i] = sets->xe[i] + I * sets->pe[i];
            }
            cc_matmul(tempcm1, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
            for (i = 0; i < seth->Nstate; i++) {
                sets->xe[i] = creal(tempcv2[i]);
                sets->pe[i] = cimag(tempcv2[i]);
            }
        }

        memset(sets->correfun_t,0,seth->Nstate * seth->Nstate * sizeof(double complex));


        if(seth->sampletype == 2){
            memset(sets->rhot_mash,0,seth->Nstate * seth->Nstate * sizeof(double complex));
            if(sets->xe[0] * sets->xe[0] + sets->pe[0] * sets->pe[0] >= 1){
                sets->rhot_mash[0] = 1;
            } else {
                sets->rhot_mash[3] = 1;
            }
            sets->rhot_mash[1] = 0.5 * (sets->xe[0] + I * sets->pe[0]) * (sets->xe[1] - I * sets->pe[1]);
            sets->rhot_mash[2] = 0.5 * (sets->xe[0] - I * sets->pe[0]) * (sets->xe[1] + I * sets->pe[1]);

            

            for (k = 0; k < seth->Nstate; k++){
                for (l = 0; l < seth->Nstate; l++){
                    if (k == l) {
                        sets->mea_mat_mash[0] = sets->measure_mash;
                        sets->mea_mat_mash[3] = sets->measure_mash;
                        sets->mea_mat_mash[1] = 2;
                        sets->mea_mat_mash[2] = 2;
                    } else {
                        sets->mea_mat_mash[0] = 2;
                        sets->mea_mat_mash[3] = 2;
                        sets->mea_mat_mash[1] = 3;
                        sets->mea_mat_mash[2] = 3;
                    }

                    for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                        tempcm2[i] = sets->rhot_mash[i] * sets->mea_mat_mash[i];
                    }
                    // transpose(sets->U_d2a,tempdm,seth->Nstate);
                    // cd_matmul(tempcm2,tempdm,tempcm,seth->Nstate,seth->Nstate,seth->Nstate);
                    // dc_matmul(sets->U_d2a,tempcm,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                    transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
                    cc_matmul(tempcm2, tempcm1, tempcm3, seth->Nstate, seth->Nstate, seth->Nstate);
                    cc_matmul(sets->U_d2a, tempcm3, tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
                   

                    for (i = 0; i < seth->Nstate; i++){
                        for (j = 0; j < seth->Nstate; j++){
                            sets->correfun_t[i * seth->Nstate + j] += 
                            2 * sets->U0_mash[(sets->init_occ - 1) * seth->Nstate + k] * conj(sets->U0_mash[(sets->init_occ - 1) * seth->Nstate + l])
                            * sets->rho0_mash[k * seth->Nstate + l] * tempcm2[i * seth->Nstate + j];
                        }
                    }
                }
            }

           
        }
        else{
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    if(i == j){
                        if(sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i] >= 1){
                            sets->correfun_t[i * seth->Nstate + i] = 2 * sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
                        } else{
                            sets->correfun_t[i * seth->Nstate + i] = 0.0;
                        }
                    }
                    else{
                        sets->correfun_t[i * seth->Nstate + j] = 2 * 2 * sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                    }
                }
            }

        }



        if(seth->sampletype == 2){
            // memcpy(tempv,sets->xe,seth->Nstate * sizeof(double));
            // dd_matmul(sets->U_d2a,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
            // memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
            // dd_matmul(sets->U_d2a,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
            for (i = 0; i < seth->Nstate; i++) {
                tempcv1[i] = sets->xe[i] + I * sets->pe[i];
            }
            cc_matmul(sets->U_d2a, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
            for (i = 0; i < seth->Nstate; i++) {
                sets->xe[i] = creal(tempcv2[i]);
                sets->pe[i] = cimag(tempcv2[i]);
            }
        }
    

    } else if (strcmp(seth->method, "MASH-MF") == 0 || strcmp(seth->method, "mash-mf") == 0 ||
               strcmp(seth->method, "mashmf") == 0 || strcmp(seth->method, "MASHMF") == 0 ||
               strcmp(seth->method, "mr") == 0 || strcmp(seth->method, "MR") == 0) {
        
        
        for (i = 0; i < seth->Nstate; i++) {
            for (j = 0; j < seth->Nstate; j++) {
                if(i == j){
                    if(sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i] >= 1){
                        sets->correfun_t[i * seth->Nstate + i] = 2 * sets->correfun_0 * sets->measure_mash;
                    } else{
                        sets->correfun_t[i * seth->Nstate + i] = 0.0;
                    }
                }
                else{
                    sets->correfun_t[i * seth->Nstate + j] = 2 * 2 * sets->correfun_0 * 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                }
            }
        }

        


    
    } else if (strcmp(seth->method, "MS-MASH") == 0 || strcmp(seth->method, "ms-mash") == 0 ||
               strcmp(seth->method, "msmash") == 0 || strcmp(seth->method, "MSMASH") == 0 ||
               strcmp(seth->method, "MASH-RM") == 0 || strcmp(seth->method, "mash-rm") == 0 ||
               strcmp(seth->method, "MS-MASH-MF") == 0 || strcmp(seth->method, "ms-mash-mf") == 0 ||
               strcmp(seth->method, "msmashmf") == 0 || strcmp(seth->method, "MSMASHMF") == 0 ||
               strcmp(seth->method, "RM") == 0 || strcmp(seth->method, "rm") == 0 ||
                strcmp(seth->method, "MA-NAF-RM") == 0 || strcmp(seth->method, "ma-naf-rm") == 0) {
        
        alpha_mash = 0;
        for (i = 1; i < seth->Nstate + 1; i++){
            alpha_mash += 1.0 / ((double) i); 
        }
        
        alpha_mash = (double) (seth->Nstate - 1.0) / (alpha_mash - 1.0);
        beta_mash = (1.0 - alpha_mash) / seth->Nstate;


        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->correfun_t[i * seth->Nstate + j] = alpha_mash * sets->den_e[i * seth->Nstate + j] + beta_mash * (i == j ? 1 : 0);
                }
            }
        } else {      
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    sets->correfun_t[i * seth->Nstate + j] = alpha_mash * 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) + beta_mash * (i == j ? 1 : 0);
                }
            }
        }

        // printf("%f %f %f\n",creal(sets->correfun_t[0]),creal(sets->correfun_t[1]),creal(sets->correfun_t[3]));

    } else if (strcmp(seth->method, "msmash2") == 0 || strcmp(seth->method, "MSMASH2") == 0) {
        
        // alpha_mash = 0;
        // for (i = 1; i < seth->Nstate + 1; i++){
        //     alpha_mash += 1.0 / ((double) i); 
        // }
        
        // alpha_mash = (double) (seth->Nstate - 1.0) / (alpha_mash - 1.0);
        // beta_mash = (1.0 - alpha_mash) / seth->Nstate;

        for (int i = 0; i < seth->Nstate; i++) {
            if (seth->type_evo == 0) {
                c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
            } else if (seth->type_evo == 1) {
                c_main[i] = creal(sets->den_e[i * seth->Nstate + i]);
            }
        }
        i_st = maxloc(c_main, seth->Nstate);



        if (seth->type_evo == 1 || seth->type_evo == 3) {
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    if (i == i_st || j == i_st){ 
                        sets->correfun_t[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j];
                        if (i == j) sets->correfun_t[i * seth->Nstate + j] *= 2;
                    } else {
                        sets->correfun_t[i * seth->Nstate + j] = 0.0;
                    }
                }
            }
        } else {      
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    if (i == i_st || j == i_st){ 
                        sets->correfun_t[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                        if (i == j) sets->correfun_t[i * seth->Nstate + j] *= 2;
                    } else {
                        sets->correfun_t[i * seth->Nstate + j] = 0.0;
                    }
                }
            }
        }

        // printf("%f %f %f\n",creal(sets->correfun_t[0]),creal(sets->correfun_t[1]),creal(sets->correfun_t[3]));
    } else if (strcmp(seth->method, "msmash3") == 0 || strcmp(seth->method, "MSMASH3") == 0) {
        
        // alpha_mash = 0;
        // for (i = 1; i < seth->Nstate + 1; i++){
        //     alpha_mash += 1.0 / ((double) i); 
        // }
        
        // alpha_mash = (double) (seth->Nstate - 1.0) / (alpha_mash - 1.0);
        // beta_mash = (1.0 - alpha_mash) / seth->Nstate;

        // if (seth->rep == 0) {
            
            
        //     memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
        //     transpose(sets->U_d2a,tempdm1,seth->Nstate);
        //     dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        //     cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            
        //     memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
        //     dc_matmul(tempdm1,gamma_cv_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        //     cd_matmul(tempcm1,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        // }

        // for (int i = 0; i < seth->Nstate; i++) {
        //     if (seth->type_evo == 0) {
        //         c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
        //     } else if (seth->type_evo == 1) {
        //         c_main[i] = creal(sets->den_e[i * seth->Nstate + i]);
        //     }
        // }
        // i_st = maxloc(c_main, seth->Nstate);



        // if (seth->type_evo == 1 || seth->type_evo == 3) {
        //     for (i = 0; i < seth->Nstate; i++) {
        //         for (j = 0; j < seth->Nstate; j++) {
        //             if (i == i_st || j == i_st){ 
        //                 sets->correfun_t[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j];
        //                 if (i == j) sets->correfun_t[i * seth->Nstate + j] *= 2;
        //             } else {
        //                 sets->correfun_t[i * seth->Nstate + j] = 0.0;
        //             }
        //         }
        //     }
        // } else {      
        //     for (i = 0; i < seth->Nstate; i++) {
        //         for (j = 0; j < seth->Nstate; j++) {
        //             if (i == i_st || j == i_st){ 
        //                 sets->correfun_t[i * seth->Nstate + j] = 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
        //                 if (i == j) sets->correfun_t[i * seth->Nstate + j] *= 2;
        //             } else {
        //                 sets->correfun_t[i * seth->Nstate + j] = 0.0;
        //             }
        //         }
        //     }
        // }


        // if (seth->rep == 0) {
        //     transpose(sets->U_d2a,tempdm1,seth->Nstate);
        //     dc_matmul(sets->U_d2a,sets->correfun_t,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
        //     cd_matmul(tempcm1,tempdm1,sets->correfun_t,seth->Nstate,seth->Nstate,seth->Nstate);

        //     if (seth->type_evo == 0) {
        //         memcpy(sets->xe,xe_save,seth->Nstate*sizeof(double));
        //         memcpy(sets->pe,pe_save,seth->Nstate*sizeof(double));
        //     } else if (seth->type_evo == 1) {
        //         memcpy(sets->den_e,den_e_save,seth->Nstate * seth->Nstate * sizeof(double complex));
        //     }
        //     memcpy(sets->gamma_cv,gamma_cv_save,seth->Nstate * seth->Nstate * sizeof(double complex));
        
        // }


    } else if (strcmp(seth->method, "ms-mash-mf2") == 0 || strcmp(seth->method, "MS-MASH-MF2") == 0 ||
               strcmp(seth->method, "mf2") == 0 || strcmp(seth->method, "MF2") == 0 ||
               strcmp(seth->method, "CW1") == 0 || strcmp(seth->method, "cw1") == 0 ) {
        
        for (int i = 0; i < seth->Nstate; i++) {
            if (seth->type_evo == 0) {
                c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
            } else if (seth->type_evo == 1) {
                c_main[i] = creal(sets->den_e[i * seth->Nstate + i]);
            }
        }
        i_st = maxloc(c_main, seth->Nstate);

        for (i = 0; i < seth->Nstate; i++) {
            for (j = 0; j < seth->Nstate; j++) {
                if (i == j){
                    if (i == i_st) {
                        sets->correfun_t[i * seth->Nstate + i] = 1.0;
                    } else {
                        sets->correfun_t[i * seth->Nstate + i] = 0.0;
                    }
                } else {
                    if (seth->type_evo == 0) {
                        sets->correfun_t[i * seth->Nstate + j] = (1.0 + seth->Nstate) / (2 * pow(1 + seth->Nstate * seth->gamma_zpe, 2)) * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                    } else if (seth->type_evo == 1){
                        sets->correfun_t[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j] * (1.0 + seth->Nstate) / pow(1 + seth->Nstate * seth->gamma_zpe, 2);
                    }
                }
                    
            }
        }

            
    } else if (strcmp(seth->method, "mf3") == 0 || strcmp(seth->method, "MF3") == 0 ||
               strcmp(seth->method, "CW2") == 0 || strcmp(seth->method, "cw2") == 0 ) {
        
        if(seth->type_evo == 0){
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    if (i == j){
                        if (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i] >= 2.0) {
                            sets->correfun_t[i * seth->Nstate + i] = pow((1.0 + seth->Nstate * seth->gamma_zpe)/(seth->Nstate * seth->gamma_zpe),seth->Nstate - 1.0) / seth->Nstate ;
                        } else {
                            sets->correfun_t[i * seth->Nstate + i] = 0.0;
                        }
                    } else {
                        sets->correfun_t[i * seth->Nstate + j] = (1.0 + seth->Nstate) / (2 * pow(1 + seth->Nstate * seth->gamma_zpe, 2)) * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]);
                    }
                        
                }
            }
        } else if (seth->type_evo == 1){
            for (i = 0; i < seth->Nstate; i++) {
                for (j = 0; j < seth->Nstate; j++) {
                    if (i == j){
                        if (creal(sets->den_e[i * seth->Nstate + i]) >= 1.0) {
                            sets->correfun_t[i * seth->Nstate + i] = pow((1.0 + seth->Nstate * seth->gamma_zpe)/(seth->Nstate * seth->gamma_zpe),seth->Nstate - 1.0) / seth->Nstate ;
                        } else {
                            sets->correfun_t[i * seth->Nstate + i] = 0.0;
                        }
                    } else {
                        sets->correfun_t[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j] * (1.0 + seth->Nstate) / pow(1 + seth->Nstate * seth->gamma_zpe, 2);
                    }
                        
                }
            }
        }
        

    } else if (strcmp(seth->method, "NW") == 0 || strcmp(seth->method, "nw") == 0 ) {
        
        for (int i = 0; i < seth->Nstate; i++) {
            if (seth->type_evo == 0) {
                c_main[i] = (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]);
            } else if (seth->type_evo == 1) {
                c_main[i] = creal(sets->den_e[i * seth->Nstate + i]);
            }
        }
        i_st = maxloc(c_main, seth->Nstate);
        alpha_mash = 0;
        for (i = 1; i < seth->Nstate + 1; i++){
            alpha_mash += 1.0 / ((double) i); 
        }
        alpha_mash = (double) (seth->Nstate - 1.0) / (alpha_mash - 1.0);

        for (i = 0; i < seth->Nstate; i++) {
            for (j = 0; j < seth->Nstate; j++) {
                if (i == j) {
                    if (i == i_st) {
                        x2 = 1.0;
                        for (k = 0; k < seth->Nstate; k++){
                            if(k == i) continue;
                            if(seth->type_evo == 0)x2 *= 0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i] - sets->xe[k] * sets->xe[k] - sets->pe[k] * sets->pe[k]) / (1.0 + seth->Nstate * seth->gamma_zpe);
                            if(seth->type_evo == 1)x2 *= creal(sets->den_e[i * seth->Nstate + i] - sets->den_e[k * seth->Nstate + k]) / (1.0 + seth->Nstate * seth->gamma_zpe);
                        }
                        sets->correfun_t[i * seth->Nstate + i] = pow(x2, (double) 3.0/(7*(seth->Nstate - 1)) + (double) 60.0/(7*(seth->Nstate + 13)));
                    } else {
                        sets->correfun_t[i * seth->Nstate + i] = 0.0; 
                    }
                } else {          
                    if(seth->type_evo == 0)sets->correfun_t[i * seth->Nstate + j] = alpha_mash * 0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) / (1.0 + seth->Nstate * seth->gamma_zpe);
                    if(seth->type_evo == 1)sets->correfun_t[i * seth->Nstate + j] = alpha_mash * sets->den_e[i * seth->Nstate + j] / (1.0 + seth->Nstate * seth->gamma_zpe);
                }
          
            }
        }

    }





    if (seth->sampletype == 2) {
        
        // transpose(sets->U_d2a,tempdm,seth->Nstate);
        // memcpy(tempv,sets->xe,seth->Nstate*sizeof(double));
        // dd_matmul(tempdm,tempv,sets->xe,seth->Nstate,seth->Nstate,1);
        // memcpy(tempv,sets->pe,seth->Nstate*sizeof(double));
        // dd_matmul(tempdm,tempv,sets->pe,seth->Nstate,seth->Nstate,1);
        transpose_conjugate(sets->U_d2a, tempcm1, seth->Nstate);
        for (i = 0; i < seth->Nstate; i++) {
            tempcv1[i] = sets->xe[i] + I * sets->pe[i];
        }
        cc_matmul(tempcm1, tempcv1, tempcv2, seth->Nstate, seth->Nstate, 1);
        for (i = 0; i < seth->Nstate; i++) {
            sets->xe[i] = creal(tempcv2[i]);
            sets->pe[i] = cimag(tempcv2[i]);
        }
       
        // memcpy(tempcm2,sets->gamma_cv,seth->Nstate*seth->Nstate*sizeof(double complex));
        // cd_matmul(sets->gamma_cv,sets->U_d2a,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
        // dc_matmul(tempdm,tempcm2,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
        cc_matmul(tempcm1, sets->gamma_cv,  tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
        cc_matmul(tempcm2, sets->U_d2a, sets->gamma_cv, seth->Nstate, seth->Nstate, seth->Nstate);

        
        if (seth->type_evo == 1 || seth->type_evo == 3) {
            // memcpy(tempcm,sets->den_e,seth->Nstate*seth->Nstate*sizeof(double complex));
            // cd_matmul(sets->den_e,sets->U_d2a,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            // dc_matmul(tempdm,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(tempcm1, sets->den_e,  tempcm2, seth->Nstate, seth->Nstate, seth->Nstate);
            cc_matmul(tempcm2, sets->U_d2a, sets->den_e, seth->Nstate, seth->Nstate, seth->Nstate);
        }
        // if (seth->type_evo == 4) {
        //     // G_xpconfg = matmul(sets->U_d2a, G_xpconfg)
        //     matmul(sets->U_d2a, G_xpconfg, seth->Nstate, seth->Nstate, 1);
        // }
    }


}


void cal_propagator(int Nstate, double complex *H, double dt, double complex *U,struct set_slave *sets,struct set_host *seth) {
    int i, j, imax;
    double E[Nstate];
    double complex C[Nstate * Nstate];
    double sineig[Nstate * Nstate], coseig[Nstate * Nstate];
    double real_pro[Nstate * Nstate], img_pro[Nstate * Nstate];
    double complex overlap[seth->Nstate * seth->Nstate];
    double overlap2[seth->Nstate * seth->Nstate], vmax;
    double tempdm1[seth->Nstate*seth->Nstate],tempdm2[seth->Nstate*seth->Nstate],tempdm3[seth->Nstate*seth->Nstate],tempdm4[seth->Nstate*seth->Nstate], tempdv1[seth->Nstate];
    double complex tempcm1[seth->Nstate*seth->Nstate],tempcm2[seth->Nstate*seth->Nstate];
    int id_max[seth->Nstate], idloc, id1, id2;
    double complex eiet[seth->Nstate * seth->Nstate];
    
    // dia_symmat(Nstate, H, E, C);
    dia_hermitemat(Nstate, H, E, C);
    
    memset(sineig,0,Nstate*Nstate*sizeof(double));
    memset(coseig,0,Nstate*Nstate*sizeof(double));
    
    for (i = 0; i < Nstate; i++) {
        sineig[i * Nstate + i] = sin(E[i] * dt);
        coseig[i * Nstate + i] = cos(E[i] * dt);
    }

    for (i = 0; i < seth->Nstate * seth->Nstate; i++){
        eiet[i] = coseig[i] - I * sineig[i];
    }
    
    transpose_conjugate(C,tempcm1,seth->Nstate);
    cc_matmul(eiet,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
    cc_matmul(C,tempcm2,U,seth->Nstate,seth->Nstate,seth->Nstate);

    
    // transpose(C, tempdm1, Nstate);
    // dd_matmul(coseig,tempdm1,tempdm2,Nstate,Nstate,Nstate);
    // dd_matmul(C,tempdm2,real_pro,Nstate,Nstate,Nstate);
    // dd_matmul(sineig,tempdm1,tempdm2,Nstate,Nstate,Nstate);
    // dd_matmul(C,tempdm2,img_pro,Nstate,Nstate,Nstate);
    
    // // printf("z  3333333 \n");
    
    // // 计算U
    // for (i = 0; i < Nstate * Nstate; i++) {
    //     U[i] = real_pro[i] + I * img_pro[i];
    // }
    //  printf("z  4444444444 \n");
    // 检查sets->E_adia是否已分配
    // printf("%d %d\n",seth->ifswitchforce, seth->rep);
    
    if (seth->ifswitchforce > 0 && seth->rep == 0 || seth->type_hop == 1 && seth->rep == 0) {
        // printf("11111\n");
        memcpy(sets->U_d2a,C,Nstate*Nstate*sizeof(double complex));
        memcpy(sets->E_adia,E,Nstate*sizeof(double));
        
        
        if (sets->if_ad_nac) {
            // transpose(sets->U_d2a,tempdm1,seth->Nstate);
            // dd_matmul(tempdm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);
            transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
            cc_matmul(tempcm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);

            memset(overlap2, 0, seth->Nstate * seth->Nstate * sizeof(double));
        
            for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                tempdm1[i] = cabs(overlap[i]);
            }

            for (i = 0; i < seth->Nstate; i++) {
                idloc=maxloc(tempdm1, seth->Nstate * seth->Nstate);
                id1 = idloc / seth->Nstate; 
                id2 = idloc % seth->Nstate; 
                overlap2[idloc] = (creal(overlap[idloc]) >= 0.0) ? 1.0 : -1.0;
                for (j = 0; j < seth->Nstate; j++) {
                    tempdm1[id1 * seth->Nstate + j] = 0;
                    tempdm1[j * seth->Nstate + id2] = 0;
                }
            }

            cd_matmul(sets->U_d2a, overlap2, tempcm1, seth->Nstate, seth->Nstate, seth->Nstate);
            memcpy(sets->U_d2a,tempcm1,seth->Nstate * seth->Nstate * sizeof(double complex));

            for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                tempdm2[i] = fabs(overlap2[i]);
            }
            dd_matmul(sets->E_adia, tempdm2, tempdv1, 1, seth->Nstate, seth->Nstate);
            memcpy(sets->E_adia, tempdv1, seth->Nstate * sizeof(double));
        }
        memcpy(sets->U_d2a_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double complex));
        memcpy(sets->U_ref, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double complex));
        sets->if_ad_nac = 1;
    }
    //  printf("z 55555  \n");
}
        
    
void evo_traj_ele(double deltat,struct set_slave *sets,struct set_host *seth, int para) {
    double x0[seth->Nstate], p0[seth->Nstate];
    double x0_mb[2 * (seth->Nstate - 1)], p0_mb[2 * (seth->Nstate - 1)];
    double complex propagator_unsmash[2 * 2 * (seth->Nstate - 1)];
    int i;
    double complex tempv1[seth->Nstate],tempv2[seth->Nstate];
    double complex tempcm1[seth->Nstate*seth->Nstate],tempcm2[seth->Nstate*seth->Nstate];
   
    switch (seth->type_evo) {
        case 0:
            
            if (seth->rep == 0) cal_propagator(seth->Nstate, sets->V, deltat, sets->propagator,sets,seth);
            
            if (seth->rep == 1) cal_propagator_adia(seth->Nstate, deltat, sets->propagator,sets,seth,para);

            if (seth->rep == 2 || seth->rep == 3) cal_propagator_gen(seth->Nstate, sets->V, deltat, sets->propagator,sets,seth);

            memcpy(x0, sets->xe, seth->Nstate * sizeof(double));
            memcpy(p0, sets->pe, seth->Nstate * sizeof(double));
            
            cd_matmul(sets->propagator,x0,tempv1,seth->Nstate,seth->Nstate,1);
            cd_matmul(sets->propagator,p0,tempv2,seth->Nstate,seth->Nstate,1);
            
            for(i=0;i<seth->Nstate;i++){
                sets->xe[i]=creal(tempv1[i])-cimag(tempv2[i]);
                sets->pe[i]=creal(tempv2[i])+cimag(tempv1[i]);
            }
           
            break;
        
        case 1:
            if (seth->rep == 0) cal_propagator(seth->Nstate, sets->V, deltat, sets->propagator,sets,seth);
            if (seth->rep == 1) cal_propagator_adia(seth->Nstate, deltat, sets->propagator,sets,seth,para);
            if (seth->rep == 2 || seth->rep == 3) cal_propagator_gen(seth->Nstate, sets->V, deltat, sets->propagator,sets,seth);
            
            transpose_conjugate(sets->propagator,tempcm1,seth->Nstate);
            cc_matmul(sets->den_e,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(sets->propagator,tempcm2,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            
            
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
        transpose_conjugate(sets->propagator,tempcm1,seth->Nstate);
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

    if (seth->ifscaleenergy == 3) energy_conserve_naf_3(deltat,sets,seth);
}

void evo_traj_nucR(double deltat,struct set_slave *sets,struct set_host *seth) {
    double x1, x2;

    
            for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
                sets->R_nuc[i] += sets->P_nuc[i] * deltat / sets->mass[i];
            }
        
}


void evo_traj_calProp(int igrid_cal,struct set_slave *sets,struct set_host *seth) {
    int i, j, icfall;
    double x2;
    double complex csum1,csum2;
    double complex tempcm1[seth->Nstate * seth->Nstate];
   
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

    if (seth->if_allcf == 0 && seth->outputtype >= 0) {
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

   

    if (seth->if_allcf == 0 && seth->outputtype != 0) {
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
                    if (strcmp(seth->msmodelname, "retinal") == 0 ) {
                        if (cos(sets->R_nuc[24]) < 0) {
                            sets->pop_fb[i * seth->Ngrid *2  + igrid_cal * 2 + 0] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                        } else {
                            sets->pop_fb[i * seth->Ngrid *2  + igrid_cal * 2 + 1] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                        } 
                    } else {
                        if (sets->P_nuc[0] > 0) {
                            sets->pop_fb[i * seth->Ngrid *2  + igrid_cal * 2 + 0] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                        } else {
                            sets->pop_fb[i * seth->Ngrid *2  + igrid_cal * 2 + 1] += creal(sets->correfun_0 * sets->correfun_t[i * seth->Nstate + i]);
                        }
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

    if (seth->if_allcf >= 2) {
        if (seth->if_allcf == 2) {
            cfweight_msmodel(sets->weight0,sets->weightt, seth->beta, sets->R_nuc_init, sets->P_nuc, 0, seth);

            if (strcmp(seth->method, "sqc") == 0 || strcmp(seth->method, "SQC") == 0) {
                // Q_{nnkl}, k \ne l
                memset(tempcm1,0,seth->Nstate*seth->Nstate*sizeof(double complex));
                for (i = 0; i < seth->Nstate; i++) {
                    tempcm1[i * seth->Nstate + i] = sets->cf0[i * seth->Nstate + i];
                }
                
                csum1 = 0, csum2 = 0;
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    csum1 += sets->weight0[i] * tempcm1[i];
                    csum2 += sets->weightt[i] * sets->correfun_t[i];
                }
                sets->cfeff[igrid_cal] += csum1 * csum2;

                // Q_{nmkl}, n \ne m
                if (seth->if_scale_sqc == 1){
                    for (i = 0; i < seth->Nstate; i++){
                        sets->xe[i] /= sqrt(sets->scale_sqc2);
                        sets->pe[i] /= sqrt(sets->scale_sqc2);
                    }
                    if (seth->type_evo >= 1){
                        for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                            sets->den_e[i] /= sqrt(sets->scale_sqc2);
                        }
                    }
                }

                for (i = 0; i < seth->Nstate ; i++) {
                    for (j = 0; j < seth->Nstate; j++){
                       sets->correfun_t[i*seth->Nstate+j] = sets->den_e[i*seth->Nstate+j];
                    }
                }


                memcpy(tempcm1,sets->cf0,seth->Nstate*seth->Nstate*sizeof(double complex));
                for (i = 0; i < seth->Nstate; i++) {
                    tempcm1[i * seth->Nstate + i] = 0;
                }
                
                csum1 = 0, csum2 = 0;
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    csum1 += sets->weight0[i] * tempcm1[i];
                    csum2 += sets->weightt[i] * sets->correfun_t[i];
                }
                sets->cfeff[igrid_cal] += csum1 * csum2;


                if (seth->if_scale_sqc == 1){
                    for (i = 0; i < seth->Nstate; i++){
                        sets->xe[i] *= sqrt(sets->scale_sqc2);
                        sets->pe[i] *= sqrt(sets->scale_sqc2);
                    }
                    if (seth->type_evo >= 1){
                        for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                            sets->den_e[i] *= sqrt(sets->scale_sqc2);
                        }
                    }
                }

            } else if (strcmp(seth->method, "NW") == 0 || strcmp(seth->method, "nw") == 0 ) {
                // Q_{nnkl}, k \ne l
                memset(tempcm1,0,seth->Nstate*seth->Nstate*sizeof(double complex));
                for (i = 0; i < seth->Nstate; i++) {
                    tempcm1[i * seth->Nstate + i] = sets->cf0[i * seth->Nstate + i];
                }
                
                csum1 = 0, csum2 = 0;
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    csum1 += sets->weight0[i] * tempcm1[i];
                    csum2 += sets->weightt[i] * sets->correfun_t[i];
                }
                sets->cfeff[igrid_cal] += csum1 * csum2;

                // Q_{nmkl}, n \ne m

                for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
                    sets->correfun_t[i] = sets->den_e[i] * (1.0 + seth->Nstate) / pow(1 + seth->Nstate * seth->gamma_zpe, 2);
                }
                for (int i = 0; i < seth->Nstate; i++) {
                    sets->correfun_t[i * seth->Nstate + i] -= (1.0 - seth->gamma_zpe) / (1 + seth->Nstate * seth->gamma_zpe);
                }


                memcpy(tempcm1,sets->cf0,seth->Nstate*seth->Nstate*sizeof(double complex));
                for (i = 0; i < seth->Nstate; i++) {
                    tempcm1[i * seth->Nstate + i] = 0;
                }
                
                csum1 = 0, csum2 = 0;
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    csum1 += sets->weight0[i] * tempcm1[i];
                    csum2 += sets->weightt[i] * sets->correfun_t[i];
                }
                sets->cfeff[igrid_cal] += csum1 * csum2;


            } else {
          
                csum1 = 0, csum2 = 0;
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    csum1 += sets->weight0[i] * sets->cf0[i];
                    csum2 += sets->weightt[i] * sets->correfun_t[i];
                }
                sets->cfeff[igrid_cal] += csum1 * csum2;

            }

        } else if (seth->if_allcf == 3) {
            for (icfall = 1; icfall < seth->allcf_times + 1; icfall++) {
                cfweight_msmodel(sets->weight0, sets->weightt, seth->beta, sets->R_nuc_init, sets->P_nuc, icfall, seth);        
                csum1 = 0, csum2 = 0;
                for (i = 0; i < seth->Nstate * seth->Nstate; i++){
                    csum1 += sets->weight0[i] * sets->cf0[i];
                    csum2 += sets->weightt[i] * sets->correfun_t[i];
                }
                sets->cfeff[igrid_cal] += csum1 * csum2;
            }
        // } else if (seth->if_allcf == 4) {
        //     for (icfall = 0; icfall < seth->allcf_times; icfall++) {
        //         cfweight_msmodel(weight0, weightt, seth->beta, icfall);
        //         sets->cfeff[igrid_cal] += sum(weight0, seth->Nstate) * sum(weightt, seth->Nstate) * sum(sets->correfun_t, seth->Nstate * seth->Nstate);
        //     }
        }
    }

    if (seth->mean_nuc == 1) {
        if (strcmp(seth->method, "sqc") == 0 ||  strcmp(seth->method, "SQC") == 0 
            || strcmp(seth->method, "mf3") == 0 || strcmp(seth->method, "MF3") == 0
            || strcmp(seth->method, "CW2") == 0 || strcmp(seth->method, "cw2") == 0 
            || strcmp(seth->method, "NW") == 0 || strcmp(seth->method, "nw") == 0 ) {
            x2 = 0;
            for (i = 0; i < seth->Nstate; i++) {
                x2 += creal(sets->correfun_t[i * seth->Nstate + i]);
            }
            for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
                sets->R_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * creal(sets->correfun_0) * x2;
                sets->P_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * creal(sets->correfun_0) * x2;
                sets->R2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0) * x2;
                sets->P2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0) * x2;
            }
        } else if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "ma-naf-mr") == 0 || strcmp(seth->method, "MA-NAF-MR") == 0) {
            for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
                // sets->R_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
                // sets->P_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
                // sets->R2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);
                // sets->P2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash * creal(sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)]);

                sets->R_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * (creal(sets->correfun_0 * sets->correfun_t[0]) + creal(sets->correfun_0 * sets->correfun_t[3]));
                sets->P_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * (creal(sets->correfun_0 * sets->correfun_t[0]) + creal(sets->correfun_0 * sets->correfun_t[3]));
                sets->R2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * sets->R_nuc[i] * (creal(sets->correfun_0 * sets->correfun_t[0]) + creal(sets->correfun_0 * sets->correfun_t[3]));
                sets->P2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * sets->P_nuc[i] * (creal(sets->correfun_0 * sets->correfun_t[0]) + creal(sets->correfun_0 * sets->correfun_t[3]));
            
            }
        // } else if (strcmp(seth->method, "unsmash") == 0 || strcmp(seth->method, "UNSMASH") == 0 || strcmp(seth->method, "unSMASH") == 0) {
        //     for (i = 0; i < Ndof; i++) {
        //         sets->R_nuc_mean[i + igrid_cal * Ndof] += sets->R_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
        //         sets->P_nuc_mean[i + igrid_cal * Ndof] += sets->P_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
        //         sets->R2_nuc_mean[i + igrid_cal * Ndof] += sets->R_nuc[i] * sets->R_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
        //         sets->P2_nuc_mean[i + igrid_cal * Ndof] += sets->P_nuc[i] * sets->P_nuc[i] * real(sets->correfun_0) * seth->Nstate * sets->measure_mash * rho0_unsmash[sets->init_occ * seth->Nstate + sets->init_occ];
        //     }
        } else if (strcmp(seth->method, "MASH-MF") == 0 || strcmp(seth->method, "mash-mf") == 0 ||
                   strcmp(seth->method, "mashmf") == 0 || strcmp(seth->method, "MASHMF") == 0 ||
                   strcmp(seth->method, "mr") == 0 || strcmp(seth->method, "MR") == 0) {
            for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
                sets->R_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
                sets->P_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
                sets->R2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
                sets->P2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0) * 2 * sets->measure_mash;
            }
        } else {
            for (i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
                sets->R_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * creal(sets->correfun_0);
                sets->P_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * creal(sets->correfun_0);
                sets->R2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->R_nuc[i] * sets->R_nuc[i] * creal(sets->correfun_0);
                sets->P2_nuc_mean[i * seth->Ngrid + igrid_cal] += sets->P_nuc[i] * sets->P_nuc[i] * creal(sets->correfun_0);
            }
        }
    }

    

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


void energy_conserve_naf_3(double deltat,struct set_slave *sets,struct set_host *seth) {
    double dE_naf, x2 = 0.0;
    for (int i = 0; i < seth->Nstate; i++) {
        for (int j = 0; j < seth->Nstate; j++) {
            if (i == j) continue;
            double sum_nac = 0.0;
            for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
                sum_nac += sets->nac[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k] * sets->P_nuc[k] / sets->mass[k];
            }
            x2 += (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i*seth->Nstate+j])) * (sets->E_adia[j] - sets->E_adia[i]) * sum_nac;
        }
    }
    double P_nuc_sum = 0.0;
    for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
        P_nuc_sum += sets->P_nuc[i] * sets->P_nuc[i] / sets->mass[i];
    }
    double factor = x2 / P_nuc_sum * deltat;
    for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
        sets->P_nuc[i] += sets->P_nuc[i] * factor;
    }
}


// exact propagator of nonadiabatic force proposed by Cheng 
void energy_conserve_naf_exact(double deltat,struct set_slave *sets,struct set_host *seth) {
    double eps = 1e-20;
    double complex Q_dia[seth->Nstate * seth->Nstate];
    double xe_save[seth->Nstate], pe_save[seth->Nstate];
    double complex gamma_cv_save[seth->Nstate * seth->Nstate], den_e_save[seth->Nstate * seth->Nstate];
    double tempdm1[seth->Nstate * seth->Nstate], tempdm2[seth->Nstate * seth->Nstate], tempdv1[seth->Nstate], tempdv2[seth->Nstate];
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate * seth->Nstate];
    double complex tempcv1[seth->Nstate], tempcv2[seth->Nstate];
    double complex commu_d_V[seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2];
    double c1 = 0.0, c2 = 0.0;
    double complex csum1, csum2;


    if (seth->Ndof1 * seth->Ndof2 == 1) { // only one nuclear degree of freedom
        return;
    }
  
    
    double K = 0.0;
    for (int i = 0; i < seth->Ndof1*seth->Ndof2; i++) {
        K += 0.5 * sets->P_nuc[i] * sets->P_nuc[i] / sets->mass[i];
    }
    double vpi[seth->Ndof1*seth->Ndof2], vb[seth->Ndof1*seth->Ndof2]; // vector \pi and vector \tilde B
    memset(vpi, 0, seth->Ndof1*seth->Ndof2 * sizeof(double));
    memset(vb, 0, seth->Ndof1*seth->Ndof2 * sizeof(double));

    if (seth->rep == 0 || seth->rep == 3) {
        if (seth->type_evo == 0) {
            memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
            memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
            // transpose(sets->U_d2a,tempdm1,seth->Nstate);
            // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
            // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
            transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
            for (int i = 0; i < seth->Nstate; i++) {
                tempcv1[i] = xe_save[i] + I * pe_save[i];
            }
            cc_matmul(tempcm1,tempcv1,tempcv2,seth->Nstate,seth->Nstate,1);
            for (int i = 0; i < seth->Nstate; i++) {
                sets->xe[i] = creal(tempcv2[i]);
                sets->pe[i] = cimag(tempcv2[i]);
            }
        } else if (seth->type_evo == 1) {
            memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
            // transpose(sets->U_d2a,tempdm1,seth->Nstate);
            // dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
            // cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
            transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
            cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
        }
        memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
        cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
        cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
                
    }
    for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
        Q_dia[i] = 0;
    }
    for (int i = 0; i < seth->Nstate; i++) {
        for (int j = 0; j < seth->Nstate; j++) {
            if (i == j) continue;
            if (seth->type_evo == 0) {
                if (seth->ifscalegamma == 0) {
                    // Q_dia[i * seth->Nstate + j] = (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + j]);
                    Q_dia[i * seth->Nstate + j] = (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * 0.5 - sets->gamma_cv[i * seth->Nstate + j];
                } else {
                    Q_dia[i * seth->Nstate + j] = (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * 0.5 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j];
                }
            } else if (seth->type_evo == 1) {
                if (seth->ifscalegamma == 0) {
                    Q_dia[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j];
                } else {
                    Q_dia[i * seth->Nstate + j] = sets->den_e[i * seth->Nstate + j] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j];
                }
            }
        }
    }
    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
    // dd_matmul(sets->U_d2a,Q_dia,tempdm2,seth->Nstate,seth->Nstate,seth->Nstate);
    // dd_matmul(tempdm2,tempdm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);
    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
    cc_matmul(sets->U_d2a,Q_dia,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
    cc_matmul(tempcm2,tempcm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);
    if (seth->rep == 0 || seth->rep == 3) {
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

    if(seth->rep == 0){
        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            vpi[k] = sets->P_nuc[k] / sqrt(sets->mass[k]);
            if(seth->calforcetype == 1) {
                for (int i = 0; i < seth->Nstate; i++) {
                    vb[k] += creal(Q_dia[i * seth->Nstate + i] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]);
                }
            } else {
                for (int i = 0; i < seth->Nstate; i++) {
                    for (int j = 0; j < seth->Nstate; j++) {
                        vb[k] += creal(Q_dia[i * seth->Nstate + j] * sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]);
                    }
                }
            }
            vb[k] = vb[k] / sqrt(sets->mass[k]);
        }

    } else if(seth->rep ==1 ){
        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            vpi[k] = sets->P_nuc[k] / sqrt(sets->mass[k]);
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    if (i == j) continue;
                    if (seth->type_evo == 0) {
                        if (seth->ifscalegamma == 0) {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];


                                vb[k] += creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - sets->gamma_cv[i * seth->Nstate + j]) *
                                                (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);

                        } else {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];

                                vb[k] += creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j]) *
                                                (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);

                        }
                    } else if (seth->type_evo == 1) {
                        if (seth->ifscalegamma == 0) {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];


                                vb[k] += creal( (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) *
                                                (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);

                        } else {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];


                                vb[k] += creal(  (sets->den_e[i * seth->Nstate + j] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j]) *
                                (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);

                        }
                    }

                }
            }
            vb[k] = vb[k] / sqrt(sets->mass[k]);
        }
    } else if(seth->rep == 2 ){

        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    csum1 = 0.0;
                    for (int l = 0; l < seth->Nstate; l++) {
                        csum1 += sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + l * seth->Ndof1 * seth->Ndof2 + k] * sets->V[l * seth->Nstate + j]
                                 - sets->V[i * seth->Nstate + l] * sets->nac[l * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                    }
                    commu_d_V[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] = csum1;
                }
            }
        }

        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            vpi[k] = sets->P_nuc[k] / sqrt(sets->mass[k]);
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    if (i == j) continue;
                    if (seth->type_evo == 0) {
                        if (seth->ifscalegamma == 0) {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];


                                vb[k] += creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - sets->gamma_cv[i * seth->Nstate + j]) *
                                                ( commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]));

                        } else {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];

                                vb[k] += creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j]) *
                                                ( commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]));

                        }
                    } else if (seth->type_evo == 1) {
                        if (seth->ifscalegamma == 0) {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];


                                vb[k] += creal( (sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) *
                                                (commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]));

                        } else {

                                // if(seth->rep==0) vb[k] += Q_dia[i * seth->Nstate + j] * sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];


                                vb[k] += creal(  (sets->den_e[i * seth->Nstate + j] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j]) *
                                (commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]));

                        }
                    }

                }
            }
            vb[k] = vb[k] / sqrt(sets->mass[k]);
        }
    }  else if(seth->rep == 3){

        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    csum1 = 0.0;
                    for (int l = 0; l < seth->Nstate; l++) {
                        csum1 += sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + l * seth->Ndof1 * seth->Ndof2 + k] * sets->V[l * seth->Nstate + j]
                                 - sets->V[i * seth->Nstate + l] * sets->nac[l * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                    }
                    commu_d_V[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] = csum1;
                }
            }
        }

        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            vpi[k] = sets->P_nuc[k] / sqrt(sets->mass[k]);
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {

                    vb[k] += creal(Q_dia[i * seth->Nstate + j] *
                            ( commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]));

                }
            }
            vb[k] = vb[k] / sqrt(sets->mass[k]);
        }
    } 
    
    

    double e_vb[seth->Ndof1*seth->Ndof2], norm_vb = 0.0;
    for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
        norm_vb += vb[k] * vb[k];
    }
    norm_vb = sqrt(norm_vb);
   
    if (norm_vb < eps) { // if \vec B is too small, using the taylor expansion
        c1 = 0.0;
        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            c1 += vb[k] * sets->P_nuc[k] / sqrt(sets->mass[k]);
        }
        c1 = c1 * deltat / (2 * K) + 1.0;
        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            sets->P_nuc[k] = c1 * sets->P_nuc[k] - deltat * sqrt(sets->mass[k]) * vb[k];
        }
        return;
    }
    memcpy(e_vb, vb, seth->Ndof1*seth->Ndof2 * sizeof(double));
    for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
        e_vb[k] = e_vb[k] / norm_vb;
    }
    
    

    double vpi_B[seth->Ndof1*seth->Ndof2], vpi_ver[seth->Ndof1*seth->Ndof2]; // vector \pi_B and vector \pi_ver

    double vdot = 0.0; // \pi_{B,0}
    for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
        vdot += vpi[k] * vb[k]/norm_vb;
    }
    
    for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
        vpi_B[k] = e_vb[k] * vdot;
        vpi_ver[k] = vpi[k] - vpi_B[k];
    }

   
    double sintheta = 1.0 - (vdot * vdot) / (2 * K);
    if(fabs(sintheta) < eps) { // \vec B // M^{-1/2}P 
        
        return;
    }


    
    // c1 = sqrt(2 * K) * (vdot - sqrt(2 * K) * tanh(norm_vb * deltat /  sqrt(2 * K))) / (sqrt(2 * K) - vdot * tanh(norm_vb * deltat / sqrt(2 * K)));
    // c2 = sqrt(2 * K) / (sqrt(2 * K) * cosh(norm_vb * deltat / sqrt(2 * K)) - vdot * sinh(norm_vb * deltat / sqrt(2 * K)));
    c1 = sqrt(2 * K) * ((vdot - sqrt(2 * K)) + (vdot + sqrt(2 * K)) * exp(-2 * norm_vb * deltat/sqrt(2 * K))) / ((sqrt(2 * K) - vdot) + (vdot + sqrt(2 * K)) * exp(-2 * norm_vb * deltat/sqrt(2 * K)));
    c2 = sqrt(2 * K) * 2 * exp(-1 * norm_vb * deltat/sqrt(2 * K)) / ((sqrt(2 * K) - vdot) + (vdot + sqrt(2 * K)) * exp(-2 * norm_vb * deltat/sqrt(2 * K)));

    for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
        sets->P_nuc[k] = c1 * sqrt(sets->mass[k]) * e_vb[k] + c2 * sqrt(sets->mass[k]) * vpi_ver[k];
    }
   

}



// P-E-R-P
void evo_traj_algorithm1(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(deltat,sets,seth,2);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}



// P-R-P-E
void evo_traj_algorithm2(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    evo_traj_ele(deltat,sets,seth,1);
    cal_force(sets,seth,2);
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucR(deltat,sets,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(0.0,sets,seth,3);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}



// E-P-R-P


void evo_traj_algorithm3(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    cal_force(sets,seth,2);
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucR(deltat,sets,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(0.0,sets,seth,3);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    evo_traj_ele(deltat,sets,seth,1);
    
}


// P-R-E-P

void evo_traj_algorithm4(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_ele(deltat,sets,seth,1);
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(0.0,sets,seth,3);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}

// P-R-E-R-P

void evo_traj_algorithm5(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucR(deltat / 2,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(deltat,sets,seth,1);
    evo_traj_nucR(deltat / 2,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(0.0,sets,seth,3);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}

// E-P-R-P-E

void evo_traj_algorithm6(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
   
    evo_traj_ele(deltat / 2,sets,seth,1);
    cal_force(sets,seth,2);
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(0.0,sets,seth,3);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    evo_traj_ele(deltat / 2,sets,seth,1);
}


// P-E-R-E-P

void evo_traj_algorithm7(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    

    
    // cal_force(sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_ele(deltat / 2,sets,seth,1);
    evo_traj_nucR(deltat,sets,seth);
    if (strcmp(seth->msmodelname, "mole") == 0 ) {
        #ifdef x86
        qm_msmodel(sets->R_nuc, seth, sets); 
        #endif
        corre_trajprop(sets, seth);
    } else {
        dV_msmodel(sets->R_nuc, sets->dV,seth);
        V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
        if (seth->rep == 1) cal_NACV(sets,seth);
        if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    }
    evo_traj_ele(deltat / 2,sets,seth,2);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}

// P-E-R-P-E
void evo_traj_algorithm8(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    

    evo_traj_ele(deltat / 2,sets,seth,1);
    cal_force(sets,seth,2);
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(deltat / 2,sets,seth,2);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    
}

// E-P-R-E-P
void evo_traj_algorithm9(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    cal_force(sets,seth,2);
    evo_traj_nucP(deltat / 2,sets,seth);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_ele(deltat / 2,sets,seth,1);   
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);
    if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    evo_traj_ele(0.0,sets,seth,3);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / 2,sets,seth);
    evo_traj_nucP(deltat / 2,sets,seth);
    evo_traj_ele(deltat / 2,sets,seth,1);
    
    
}


// (P-E)^N-R-(E-P)^N
void evo_traj_algorithm10(double deltat,struct set_slave *sets,struct set_host *seth) {
    #ifdef sunway
    int slavecore_id=athread_get_id(-1);
    #endif
    double tempv[seth->Nstate];
    double tempdm[seth->Nstate*seth->Nstate];
    
    for (int i = 0; i < seth->n_step_algom; i++) {
        cal_force(sets,seth,2);
        evo_traj_nucP(deltat / (2 * seth->n_step_algom),sets,seth);
        if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / (2 * seth->n_step_algom),sets,seth);
        evo_traj_ele(deltat / (2 * seth->n_step_algom),sets,seth,1);
    }
   
    evo_traj_nucR(deltat,sets,seth);
    dV_msmodel(sets->R_nuc, sets->dV,seth);
    V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
    if (seth->rep == 1) cal_NACV(sets,seth);

    evo_traj_ele(deltat / (2 * seth->n_step_algom),sets,seth,2);
    cal_force(sets,seth,1);
    if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / (2 * seth->n_step_algom),sets,seth);
    evo_traj_nucP(deltat / (2 * seth->n_step_algom),sets,seth);
    if (seth->n_step_algom > 1) {
        for (int i = 0; i < seth->n_step_algom - 1; i++) {
            evo_traj_ele(deltat / (2 * seth->n_step_algom),sets,seth,1);
            cal_force(sets,seth,2);
            if(seth->ifscaleenergy == 7) energy_conserve_naf_exact(deltat / (2 * seth->n_step_algom),sets,seth);
            evo_traj_nucP(deltat / (2 * seth->n_step_algom),sets,seth);
        }
    }
    

   
    
    
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
    memcpy(sets->V_old, sets->V, seth->Nstate * seth->Nstate * sizeof(double complex));
    memcpy(sets->dV_old, sets->dV, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
    // // printf("xxxxxxxxxxx\n");
    
    
    if (seth->rep == 1 || seth->rep == 2 || seth->rep == 3) {
        if (sets->E_adia != NULL) memcpy(sets->E_adia_old, sets->E_adia, seth->Nstate * sizeof(double));
        if (sets->dv_adia != NULL) memcpy(sets->dv_adia_old, sets->dv_adia, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
        if (sets->nac != NULL) memcpy(sets->nac_old, sets->nac, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
        if (sets->U_d2a != NULL) memcpy(sets->U_d2a_old, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double complex));
        if (sets->U_ref != NULL) memcpy(sets->U_ref_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double complex));
        if (sets->nac_check != NULL) memcpy(sets->nac_check_old, sets->nac_check, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
    }

    if (seth->ifswitchforce > 0 && seth->rep == 0 || seth->type_hop == 1 && seth->rep == 0){
        if (sets->E_adia != NULL) memcpy(sets->E_adia_old, sets->E_adia, seth->Nstate * sizeof(double));
        if (sets->U_d2a != NULL) memcpy(sets->U_d2a_old, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double complex));
        if (sets->U_ref != NULL) memcpy(sets->U_ref_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double complex));
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
    memcpy(sets->V, sets->V_old, seth->Nstate * seth->Nstate * sizeof(double complex));
    memcpy(sets->dV, sets->dV_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
    
    if (seth->rep == 1 || seth->rep == 2 || seth->rep == 3) {
        memcpy(sets->E_adia, sets->E_adia_old, seth->Nstate * sizeof(double));
        memcpy(sets->dv_adia, sets->dv_adia_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
        memcpy(sets->nac, sets->nac_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
        memcpy(sets->U_d2a, sets->U_d2a_old, seth->Nstate * seth->Nstate * sizeof(double complex));
        memcpy(sets->U_ref, sets->U_ref_old, seth->Nstate * seth->Nstate * sizeof(double complex));
        memcpy(sets->nac_check, sets->nac_check_old, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 *sizeof(double complex));
    }


    if (seth->ifswitchforce > 0 && seth->rep == 0 || seth->type_hop == 1 && seth->rep == 0){
        memcpy(sets->E_adia, sets->E_adia_old, seth->Nstate * sizeof(double));
        memcpy(sets->U_d2a, sets->U_d2a_old, seth->Nstate * seth->Nstate * sizeof(double complex));
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
    double R00,P00;
    // double E_diag[seth->Nstate * seth->Nstate];
    // double P_s[seth->Ndof1 * seth->Ndof2], P_s_all[seth->Nstate * seth->Ndof1 * seth->Ndof2], P_s_main[seth->Ndof1 * seth->Ndof2], R_s[seth->Ndof1 * seth->Ndof2], f_s[seth->Ndof1 * seth->Ndof2], sets->pex, pt, proj[seth->Nstate];
    double Etot, Ekin, Epot, dE;
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
    int nstep_small, istep_small;
    double t_now_small;
    // bool alive;
    int slavecore_id;
    #ifdef sunway
    slavecore_id=athread_get_id(-1);
    #endif
    FILE *traj_Rnuc = NULL;
    FILE *traj_Pnuc = NULL;
    FILE *traj_ele= NULL;
    FILE *traj_occ= NULL;
    FILE *traj_kinetic= NULL;
    FILE *traj_potential= NULL;
    FILE *traj_nac= NULL;
    FILE *traj_dv = NULL;
    FILE *traj_Ud2a = NULL;
    FILE *traj_den = NULL;
    
    // if_bak = false;
    // itime_save = 0;

    // count_sets->pertraj = 0;
    sets->if_recal_qm = 1;

    

    sets->t_now = 0;
    nstep = (int)(seth->ttot / seth->dt) + 1;

    // V_msmodel(sets->R_nuc, sets->V, 0.0,seth);
    // dV_msmodel(sets->R_nuc, sets->dV,seth);
    // if (seth->rep == 1) cal_NACV(sets,seth);
    if (strcmp(seth->msmodelname, "mole") == 0) {
        #ifdef x86
        if (seth->if_restart == 0) qm_msmodel(sets->R_nuc, seth, sets); 
        #endif   
    } else {
        dV_msmodel(sets->R_nuc, sets->dV,seth);
        V_msmodel(sets->R_nuc, sets->V, sets->t_now,seth);
        if (seth->rep == 1) cal_NACV(sets,seth);
        if (seth->rep == 2 || seth->rep == 3) nac_msmodel(sets->R_nuc, sets->nac, seth);
    }
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

    if (seth->ifscaleenergy > 0 && seth->ifscaleenergy < 7) {
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
            case 6:
                seth->scaleenergy_type = 1;
                break;
            case 3:
                seth->scaleenergy_type = 3;
                break;
        }
    }

    
    if(seth->if_restart == 0) evo_traj_savetraj(sets,seth);
    if(seth->if_restart == 0) cal_force(sets,seth,1);

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


    memcpy(sets->R_nuc_init,sets->R_nuc,seth->Ndof1 * seth->Ndof2 * sizeof(double));
   
    R00 = sets->R_nuc[0];
    P00 = sets->P_nuc[0];

    if (seth->if_printtraj == 1 && seth->if_restart == 0){
        traj_Rnuc = fopen("traj_Rnuc.xyz","w");
        traj_Pnuc = fopen("traj_Pnuc.xyz","w");
        traj_ele = fopen("traj_ele.dat","w");
        traj_occ = fopen("traj_occ.dat","w");
        traj_kinetic = fopen("traj_kinetic.dat","w");
        traj_potential = fopen("traj_potential.dat","w");
        traj_nac = fopen("traj_nac.dat","w");
        traj_dv = fopen("traj_dv.dat","w");
        if (seth->rep == 2 || seth->rep == 3) {
            traj_Ud2a = fopen("traj_Ud2a.dat","w");
        }
        traj_Ud2a = fopen("traj_Ud2a.dat","w");
        traj_den = fopen("traj_den.dat","w");
    } else if (seth->if_printtraj == 1 && seth->if_restart == 1) {
        traj_Rnuc = fopen("traj_Rnuc.xyz","a+");
        traj_Pnuc = fopen("traj_Pnuc.xyz","a+");
        traj_ele = fopen("traj_ele.dat","a+");
        traj_occ = fopen("traj_occ.dat","a+");
        traj_kinetic = fopen("traj_kinetic.dat","a+");
        traj_potential = fopen("traj_potential.dat","a+");
        traj_nac = fopen("traj_nac.dat","a+");
        traj_dv = fopen("traj_dv.dat","a+");
        if (seth->rep == 2 || seth->rep == 3) {
            traj_Ud2a = fopen("traj_Ud2a.dat","a+");
        }
        traj_Ud2a = fopen("traj_Ud2a.dat","a+");
        traj_den = fopen("traj_den.dat","a+");
        read_restart(&itime, &i_re, &igrid, sets, seth);

        printf("Restarting from time: %18.8E, itime: %d, i_re: %d, igrid: %d\n", sets->t_now, itime, i_re, igrid);
        printf("correfun_0: %18.8E %18.8E\n", creal(sets->correfun_0), cimag(sets->correfun_0));
        printf("idocc=%d\n", sets->id_state);

        sets->if_recal_qm = 0;

        // printf("R=\n");
        // for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
        //     printf("%18.8E ", sets->R_nuc[i]);
        // }
        // printf("\n");
        // printf("P=\n");
        // for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
        //     printf("%18.8E ", sets->P_nuc[i]);
        // }
        // printf("\n");

        // printf("xe=\n");
        // for (i = 0; i < seth->Nstate; i++) {
        //     printf("%18.8E ", sets->xe[i]);
        // }
        // printf("\n");
        // printf("pe=\n");
        // for (i = 0; i < seth->Nstate; i++) {
        //     printf("%18.8E ", sets->pe[i]);
        // }
        // printf("\n");
        // printf("gamma=\n");
        // for (i = 0; i < seth->Nstate; i++) {
        //     for (j = 0; j < seth->Nstate; j++) {
        //         printf("%18.8E ", creal(sets->gamma_cv[i * seth->Nstate + j]));
        //     }
        //     printf("\n");
        // }
        // for (i = 0; i < seth->Nstate; i++) {
        //     for (j = 0; j < seth->Nstate; j++) {
        //         printf("%18.8E ", cimag(sets->gamma_cv[i * seth->Nstate + j]));
        //     }
        //     printf("\n");
        // }
        // printf("V=\n");
        // for (i = 0; i < seth->Nstate; i++) {
        //     for (j = 0; j < seth->Nstate; j++) {
        //         printf("%18.8E ", creal(sets->V[i * seth->Nstate + j]));
        //     }
        //     printf("\n");
        // }
        // for (i = 0; i < seth->Nstate; i++) {
        //     for (j = 0; j < seth->Nstate; j++) {
        //         printf("%18.8E ", cimag(sets->V[i * seth->Nstate + j]));
        //     }
        //     printf("\n");
        // }
        // printf("dV=\n");
        // for (int i = 0; i < seth->Natom_mole * 3; i++){
        //     for (int j = 0; j < seth->Nstate; j++){
        //         if (seth->rep == 1){
        //             printf("%18.8E", creal(sets->dv_adia[j * seth->Nstate *  seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + i]));
        //         } else if (seth->rep == 2 || seth->rep == 3){
        //             printf("%18.8E", creal(sets->dV[j * seth->Nstate *  seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + i]));
        //         }
                
        //     }
        //     printf("\n");
        // }
        // printf("nac=\n");
        // for (int i = 0; i < seth->Natom_mole * 3; i++){
        //     for (int j = 0; j < seth->Nstate * seth->Nstate; j++){
        //            printf("%18.8E", creal(sets->nac[j * seth->Natom_mole * 3 + i]));
        //     }
        //     printf( "\n");
        // }
        

        // printf( "force:\n");
        // for (int i = 0; i < seth->Natom_mole * 3; i++){
        //     printf( "%18.8E",sets->force[i]);
        // }
        // printf( "\n");


        // if (seth->rep == 2 || seth->rep == 3) {
        //     printf( "Eadia:\n");
        //     for (int i = 0; i < seth->Nstate; i++){
        //         printf( "%18.8E", sets->E_adia[i]);
        //     }
        //     printf( "\n");
        //     printf("Ud2a:\n");
        //     for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
        //         printf("%18.8E", creal(sets->U_d2a[i]));
        //     }
        //     for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
        //         printf("%18.8E", cimag(sets->U_d2a[i]));
        //     }
        //     printf("\n");
        //     // fprintf(restart, "Uref:\n");
        //     // for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
        //     //     fprintf(restart, "%18.8E", creal(sets->U_ref[i]));
        //     // }
        //     // for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
        //     //     fprintf(restart, "%18.8E", cimag(sets->U_ref[i]));
        //     // }
        //     // fprintf(restart,"\n");
        // }



        // exit(-1);

        evo_traj_savetraj(sets,seth);
        cal_force(sets,seth,1);

    }
 
    while (itime <= nstep) {
        if (i_re >= seth->Nbreak && igrid < seth->Ngrid) {
            // printf("%18.8E %18.8E %18.8E\n",sets->t_now,sets->R_nuc[0],sets->P_nuc[0]);
            
            evo_traj_calProp(igrid,sets,seth);
             
            // if(slavecore_id == 10) printf("%18.8E %18.8E %18.8E\n", sets->t_now, sets->R_nuc[0], sets->P_nuc[0]);
             
            sets->timegrid[igrid] = sets->t_now;
            igrid++;
            i_re = 0;
            
            if (seth->if_printtraj == 1){
                
                print_traj(traj_Rnuc, traj_Pnuc, traj_ele, traj_occ, 
                            traj_kinetic, traj_potential, traj_nac, traj_den, traj_dv, traj_Ud2a, 
                            sets, seth);
                
                print_restart(itime, i_re, igrid, sets, seth);
                
                
            }

        }
        
        // printf("%18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %d %18.8E %18.8E\n",sets->t_now, sets->R_nuc[0],sets->P_nuc[0],
        //                         sets->xe[0],sets->xe[1],sets->pe[0],sets->pe[1],sets->id_state+1,
        //                         0.5*(sets->xe[0]*sets->xe[0]+sets->pe[0]*sets->pe[0])/sets->scale_sqc2,
        //                         0.5*(sets->xe[1]*sets->xe[1]+sets->pe[1]*sets->pe[1])/sets->scale_sqc2);
        if(seth->if_flighttime_tully == 1){
            if((sets->R_nuc[0] < -1.0 * seth->Xb_tully && sets->P_nuc[0] < 0) || 
               (sets->R_nuc[0] > seth->Xb_tully && sets->P_nuc[0] > 0) ) {    
                break;
            }
        }
     

        evo_traj_savetraj(sets,seth);

        dt_evo = seth->dt;

        // if(seth->if_flighttime_tully == 1){
        //     if(fabs(sets->R_nuc[0]+sets->P_nuc[0]/sets->mass[0])>2){
        //         // dt_evo=10;
        //     }
        // }
        

        switch (seth->type_algorithm) {
            case 1:
            
                evo_traj_algorithm1(dt_evo,sets,seth);
                break;
            case 2:
                evo_traj_algorithm2(dt_evo,sets,seth);
                break;
            case 3:
                evo_traj_algorithm3(dt_evo,sets,seth);
                break;
            case 4:
                evo_traj_algorithm4(dt_evo,sets,seth);
                break;
            case 5:
                evo_traj_algorithm5(dt_evo,sets,seth);
                break;
            case 6:
                evo_traj_algorithm6(dt_evo,sets,seth);
                break;
            case 7:
                evo_traj_algorithm7(dt_evo,sets,seth);
                break;
            case 8:
                evo_traj_algorithm8(dt_evo,sets,seth);
                break;
            case 9:
                evo_traj_algorithm9(dt_evo,sets,seth);
                break;
            case 10:
                evo_traj_algorithm10(dt_evo,sets,seth);
                break;
        }

        //debug
        // if(slavecore_id == 0) printf("itime=%d,sets->force=%18.8E\n",itime,sets->force[299]);

  
       
        if (seth->ifscaleenergy > 0) {
            if (seth->scaleenergy_type == 1) energy_conserve_naf_1(sets->E_conserve, &deltaE, sets, seth);

            switch (seth->ifscaleenergy) {
                case 4:
                    if (seth->scaleenergy_type == 1 && deltaE < 0) {
                        nstep_small = 2;
                        dt_evo /= 2;
                        t_now_small = 0;

                    REDO:
                        //debug
                        // printf("before %18.8e %18.8e %18.8e %18.8e %18.8e %18.8e %18.8e\n",sets->R_nuc[0],sets->P_nuc[0],dt_evo, deltaE, sets->E_conserve, sets->E_adia[sets->id_state],sets->E_conserve-sets->E_adia[sets->id_state]);

                        evo_traj_back(sets,seth);
                        for (istep_small = 1; istep_small <= nstep_small; ++istep_small) {
                            switch (seth->type_algorithm) {
                                case 1:
                                    evo_traj_algorithm1(dt_evo,sets,seth);
                                    break;
                                case 2:
                                    evo_traj_algorithm2(dt_evo,sets,seth);
                                    break;
                                case 3:
                                    evo_traj_algorithm3(dt_evo,sets,seth);
                                    break;
                                case 4:
                                    evo_traj_algorithm4(dt_evo,sets,seth);
                                    break;
                                case 5:
                                    evo_traj_algorithm5(dt_evo,sets,seth);
                                    break;
                                case 6:
                                    evo_traj_algorithm6(dt_evo,sets,seth);
                                    break;
                                case 7:
                                    evo_traj_algorithm7(dt_evo,sets,seth);
                                    break;
                                case 8:
                                    evo_traj_algorithm8(dt_evo,sets,seth);
                                    break;
                                case 9:
                                    evo_traj_algorithm9(dt_evo,sets,seth);
                                    break;
                                case 10:
                                    evo_traj_algorithm10(dt_evo,sets,seth);
                                    break;
                            }

                            if (seth->scaleenergy_type == 1) {
                                energy_conserve_naf_1(sets->E_conserve, &deltaE, sets, seth);
                            }
                            // printf("after %18.8e %18.8e %18.8e %18.8e %18.8e %18.8e %18.8e\n",sets->R_nuc[0],sets->P_nuc[0],dt_evo, deltaE, sets->E_conserve, sets->E_adia[sets->id_state],sets->E_conserve-sets->E_adia[sets->id_state]);


                            if (seth->scaleenergy_type == 1 && deltaE < 0) {
                                if (dt_evo > seth->dt / 1024) {
                                    nstep_small *= 2;
                                    dt_evo /= 2;
                                    goto REDO;
                                }
                            }
                            
                        }
                        // dt_evo = seth->dt;
                    }
                   
                    if (seth->scaleenergy_type == 3) {
                        seth->scaleenergy_type = 1;
                    }
                    break;

                case 5:

                    if (seth->scaleenergy_type == 1 && deltaE < 0) {
                        evo_traj_back(sets,seth);
                        for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++){
                            sets->P_nuc[i] = -1.0 * sets->P_nuc[i];
                        }
                        switch (seth->type_algorithm) {
                            case 1:
                                evo_traj_algorithm1(dt_evo,sets,seth);
                                break;
                            case 2:
                                evo_traj_algorithm2(dt_evo,sets,seth);
                                break;
                            case 3:
                                evo_traj_algorithm3(dt_evo,sets,seth);
                                break;
                            case 4:
                                evo_traj_algorithm4(dt_evo,sets,seth);
                                break;
                            case 5:
                                evo_traj_algorithm5(dt_evo,sets,seth);
                                break;
                            case 6:
                                evo_traj_algorithm6(dt_evo,sets,seth);
                                break;
                            case 7:
                                evo_traj_algorithm7(dt_evo,sets,seth);
                                break;
                            case 8:
                                evo_traj_algorithm8(dt_evo,sets,seth);
                                break;
                            case 9:
                                evo_traj_algorithm9(dt_evo,sets,seth);
                                break;
                            case 10:
                                evo_traj_algorithm10(dt_evo,sets,seth);
                                break;
                        }

                        if (seth->scaleenergy_type == 1) {
                            energy_conserve_naf_1(sets->E_conserve, &deltaE, sets, seth);
                        }
                    }
                
                case 6:

                    if (seth->scaleenergy_type == 1 && deltaE < 0) {
                        evo_traj_back(sets,seth);
                        for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++){
                            sets->P_nuc[i] = 0.0;
                        }
                        switch (seth->type_algorithm) {
                            case 1:
                                evo_traj_algorithm1(dt_evo,sets,seth);
                                break;
                            case 2:
                                evo_traj_algorithm2(dt_evo,sets,seth);
                                break;
                            case 3:
                                evo_traj_algorithm3(dt_evo,sets,seth);
                                break;
                            case 4:
                                evo_traj_algorithm4(dt_evo,sets,seth);
                                break;
                            case 5:
                                evo_traj_algorithm5(dt_evo,sets,seth);
                                break;
                            case 6:
                                evo_traj_algorithm6(dt_evo,sets,seth);
                                break;
                            case 7:
                                evo_traj_algorithm7(dt_evo,sets,seth);
                                break;
                            case 8:
                                evo_traj_algorithm8(dt_evo,sets,seth);
                                break;
                            case 9:
                                evo_traj_algorithm9(dt_evo,sets,seth);
                                break;
                            case 10:
                                evo_traj_algorithm10(dt_evo,sets,seth);
                                break;
                        }

                        if (seth->scaleenergy_type == 1) {
                            energy_conserve_naf_1(sets->E_conserve, &deltaE, sets, seth);
                        }
                    }
            }
        }
        
        // Ekin = 0.0;
        // for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
        //     Ekin += sets->P_nuc[i] * sets->P_nuc[i] / sets->mass[i];
        // }
        // Ekin *= 0.5;
        // Epot = sets->E_adia[sets->id_state];
        // Epot = sets->V[sets->id_state * seth->Nstate + sets->id_state];
        
        // Epot = 0.0;
        // for (i = 0; i < seth->Nstate; i++) {
        //     for (j = 0; j < seth->Nstate; j++) {
        //         Epot += creal(0.5 * (sets->xe[i] + I * sets->xe[i]) * (sets->xe[j] - I * sets->xe[j]) * sets->V[i * seth->Nstate + j]);
        //     }
        // }
        // Etot = Ekin + Epot;
        // printf("%18.8E %18.8E \n", sets->t_now, Etot);


        // if(seth->if_flighttime_tully == 1){
        //     sets->t_now += dt_evo;
        // }else{
            sets->t_now = (itime + 1) * seth->dt;
        // }
        i_re++;
        itime++;
        // if (ifzsets->pecorr > 0) zsets->pecorr_msmodel(sets->P_nuc, sets->R_nuc, ifzsets->pecorr);

        if (strcmp(seth->msmodelname, "morse3") == 0 || strcmp(seth->msmodelname, "Morse3") == 0) {
            if (seth->ifhardwall == 1) {
                if (sets->P_nuc[0] < 0 && sets->R_nuc[0] < 0) sets->P_nuc[0] = -sets->P_nuc[0];
            }
        }

        sets->if_recal_qm = 1;
        
    }


    if(seth->if_flighttime_tully == 1){
        // if((sets->R_nuc[0] < -1.0 * seth->Xb_tully && sets->P_nuc[0] < 0) || 
        //    (sets->R_nuc[0] > seth->Xb_tully && sets->P_nuc[0] > 0) ) {
        #ifdef sunway
            // for (int i = 0; i < 64; i++){
            //     if(i == slavecore_id){
            //         printf("%18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %d\n",
            //         sets->t_now,
            //         R00,P00,
            //         sets->R_nuc[0],sets->P_nuc[0],
            //         creal(sets->correfun_0 * sets->correfun_t[0]),
            //         creal(sets->correfun_0 * sets->correfun_t[3]),
            //         sets->id_state + 1);
            //         } 
                
            // }
            // athread_ssync_array();
            // printf(
            // "%18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %d\n",
            // // "%18.8E %18.8E %18.8E %d\n",
            //     sets->t_now,
            //     R00,P00,
            //     sets->R_nuc[0],sets->P_nuc[0],
            //     creal(sets->correfun_0 * sets->correfun_t[0]),
            //     creal(sets->correfun_0 * sets->correfun_t[3]),
            //     sets->id_state + 1);
            // for (int i = 0; i < 64; i++){
            //     if(i == slavecore_id){
                    seth->save_flighttime[0*seth->Ntraj/seth->mpi_size+itraj-1]=sets->t_now;
                    seth->save_flighttime[1*seth->Ntraj/seth->mpi_size+itraj-1]=R00;
                    seth->save_flighttime[2*seth->Ntraj/seth->mpi_size+itraj-1]=P00;
                    seth->save_flighttime[3*seth->Ntraj/seth->mpi_size+itraj-1]=sets->R_nuc[0];
                    seth->save_flighttime[4*seth->Ntraj/seth->mpi_size+itraj-1]=sets->P_nuc[0];
                    seth->save_flighttime[5*seth->Ntraj/seth->mpi_size+itraj-1]=creal(sets->correfun_0 * sets->correfun_t[0]);
                    seth->save_flighttime[6*seth->Ntraj/seth->mpi_size+itraj-1]=creal(sets->correfun_0 * sets->correfun_t[3]);
                    seth->save_flighttime[7*seth->Ntraj/seth->mpi_size+itraj-1]=(float)(sets->id_state + 1);
            //     }
            //     athread_ssync_array(); 
            // }
        // }
        #elif defined(x86)

            // printf(
            // "%18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %18.8E %d\n",
            // // "%18.8E %18.8E %18.8E %d\n",
            //     sets->t_now,
            //     R00,P00,
            //     sets->R_nuc[0],sets->P_nuc[0],
            //     creal(sets->correfun_0 * sets->correfun_t[0]),
            //     creal(sets->correfun_0 * sets->correfun_t[3]),
            //     sets->id_state + 1);
            seth->save_flighttime[0*seth->Ntraj/seth->mpi_size+itraj-1]=sets->t_now;
            seth->save_flighttime[1*seth->Ntraj/seth->mpi_size+itraj-1]=R00;
            seth->save_flighttime[2*seth->Ntraj/seth->mpi_size+itraj-1]=P00;
            seth->save_flighttime[3*seth->Ntraj/seth->mpi_size+itraj-1]=sets->R_nuc[0];
            seth->save_flighttime[4*seth->Ntraj/seth->mpi_size+itraj-1]=sets->P_nuc[0];
            seth->save_flighttime[5*seth->Ntraj/seth->mpi_size+itraj-1]=creal(sets->correfun_0 * sets->correfun_t[0]);
            seth->save_flighttime[6*seth->Ntraj/seth->mpi_size+itraj-1]=creal(sets->correfun_0 * sets->correfun_t[3]);
            seth->save_flighttime[7*seth->Ntraj/seth->mpi_size+itraj-1]=(float)(sets->id_state + 1);

        #endif
    }
     
    

    if (seth->if_Pdis == 1) {
        if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "ma-naf-mr") == 0 || strcmp(seth->method, "MA-NAF-MR") == 0) {
    
            for(i = 0; i < seth->s_N; i++){
                // sets->expisp[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * 2 * sets->rho0_mash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
                sets->expisp[i] += (cos(sets->P_nuc[0] * seth->s_grid[i]) + I * sin(sets->P_nuc[0] * seth->s_grid[i])) * (creal(sets->correfun_0 * sets->correfun_t[0]) + creal(sets->correfun_0 * sets->correfun_t[3]));
            }
                // break;
        // } else if (strcmp(seth->method, "unsmash") == 0 ||
        //           strcmp(seth->method, "UNSMASH") == 0 ||
        //           strcmp(seth->method, "unSMASH") == 0 ){
        //     for(i = 0; i < seth->s_N; i++){
        //         // sets->expisp[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * seth->Nstate * rho0_unsmash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
        //         sets->expisp[i] += sets->correfun_0 * (cos(sets->P_nuc[0] * sets->s[i]) + I * sin(sets->P_nuc[0] * sets->s[i])) * seth->Nstate * sets->rho0_unsmash[(sets->init_occ-1) * seth->Nstate + (sets->init_occ-1)] * sets->measure_mash;
        //     }
        //         // break;
        } else if (strcmp(seth->method, "MASH-MF") == 0 || strcmp(seth->method, "mash-mf") == 0 ||
                   strcmp(seth->method, "mashmf") == 0 || strcmp(seth->method, "MASHMF") == 0 ||
                   strcmp(seth->method, "mr") == 0 || strcmp(seth->method, "MR") == 0) {
            for(i = 0; i < seth->s_N; i++){
                // sets->expisp[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * 2 * sets->measure_mash;
                sets->expisp[i] += sets->correfun_0 * (cos(sets->P_nuc[0] * seth->s_grid[i]) + I* sin(sets->P_nuc[0] * seth->s_grid[i])) * 2 * sets->measure_mash;
            }
                // break;
        } else if (strcmp(seth->method, "sqc") == 0 ||  strcmp(seth->method, "SQC") == 0 
            || strcmp(seth->method, "mf3") == 0 || strcmp(seth->method, "MF3") == 0
            || strcmp(seth->method, "CW2") == 0 || strcmp(seth->method, "cw2") == 0 
            || strcmp(seth->method, "NW") == 0 || strcmp(seth->method, "nw") == 0 ) {
            x2 = 0;
            for (i = 0; i < seth->Nstate; i++) {
                x2 += creal(sets->correfun_t[i * seth->Nstate + i]);
            }
            for(i = 0; i < seth->s_N; i++){
                // sets->expisp[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]) * x2;
                sets->expisp[i] += sets->correfun_0 * ( cos( sets->P_nuc[0] * seth->s_grid[i]) + I * sin( sets->P_nuc[0] * seth->s_grid[i]) ) * x2;
            }
                // break;

        } else {
                for(i = 0; i < seth->s_N; i++){
                    // sets->expisp[i] += sets->correfun_0 * cexp(I * sets->P_nuc[0] * s[i]);
                    sets->expisp[i] += sets->correfun_0 * ( cos( sets->P_nuc[0] * seth->s_grid[i]) + I * sin( sets->P_nuc[0] * seth->s_grid[i]) );
                }
                // break;
        }
    }
    



    if (seth->if_printtraj == 1){
        fclose(traj_Rnuc);
        fclose(traj_Pnuc); 
        fclose(traj_ele);
        fclose(traj_occ);
        fclose(traj_kinetic);
        fclose(traj_potential);
        fclose(traj_nac);
        fclose(traj_dv);
        if (seth->rep == 2 || seth->rep == 3) {
            fclose(traj_Ud2a);
        }
        fclose(traj_den);
    }


}


void cal_force(struct set_slave *sets,struct set_host *seth,int para) {
    int iref, i, j;
    double frdm, x2;

    // sets->force = 0.0;
    
    memset(sets->force,0,seth->Ndof1*seth->Ndof2*sizeof(double));

    if (strcmp(seth->method, "FSSH") == 0 || strcmp(seth->method, "fssh") == 0 ||
        strcmp(seth->method, "GFSH") == 0 || strcmp(seth->method, "gfsh") == 0 ) {
        cal_force_sh(sets,seth,para);
        // sets->force = -sets->dv_adia[sets->id_state][sets->id_state];
    } else if (strcmp(seth->method, "MASH") == 0 || strcmp(seth->method, "mash") == 0 ||
               strcmp(seth->method, "mash-mr") == 0 || strcmp(seth->method, "MASH-MR") == 0 ||
               strcmp(seth->method, "MS-MASH") == 0 || strcmp(seth->method, "ms-mash") == 0 ||
               strcmp(seth->method, "msmash") == 0 || strcmp(seth->method, "MSMASH") == 0 ||
               strcmp(seth->method, "MASH-RM") == 0 || strcmp(seth->method, "mash-rm") == 0 ||
               strcmp(seth->method, "msmash2") == 0 || strcmp(seth->method, "MSMASH2") == 0 ||
               strcmp(seth->method, "msmash3") == 0 || strcmp(seth->method, "MSMASH3") == 0 ) {

        cal_force_sh(sets,seth,para);

    } else {
        
        // if (if_ref == 1) {
        //     cal_force_mf_ref();
        // } else if (if_1st > 0) {
        //     cal_force_mf_1st();
        // } else if (ifmsbranch > 0) {
        //     cal_force_msbranch();
        if (seth->ifswitchforce == 1) {
            if(seth->ifscaleenergy == 7 || seth->ifscaleenergy == 8){
                cal_force_sh(sets,seth,para);
            } else {
                cal_force_switch(sets,seth,para);
            }
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
        } else if (seth->if_eld == 1) {
            cal_force_eld(sets,seth);
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
    double complex commu_d_V[seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2];
    double complex csum1, csum2;
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
                            sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - sets->gamma_cv[i * seth->Nstate + i]));
                        }
                    }
                } else {
                    for (i = 0; i < seth->Nstate; i++) {
                        for (j = 0; j < seth->Nstate; j++) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * ((sets->xe[j] + I * sets->pe[j]) * (sets->xe[i] - I * sets->pe[i]) * 0.5 - sets->gamma_cv[j * seth->Nstate + i]));
                            }
                        }
                    }
                }
            } else if (seth->rep == 1) {
                for (i = 0; i < seth->Nstate; i++) {
                    for (j = 0; j < seth->Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal((0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - creal(sets->gamma_cv[i * seth->Nstate + i])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]);
                                // printf("sets->force(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
                            }
                        } else {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);
                                // sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                            }
                        }
                    }
                }
            } else if (seth->rep == 2) {

                for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
                    for (int i = 0; i < seth->Nstate; i++) {
                        for (int j = 0; j < seth->Nstate; j++) {
                            csum1 = 0.0;
                            for (int l = 0; l < seth->Nstate; l++) {
                                csum1 += sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + l * seth->Ndof1 * seth->Ndof2 + k] * sets->V[l * seth->Nstate + j]
                                         - sets->V[i * seth->Nstate + l] * sets->nac[l * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                            }
                            commu_d_V[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] = csum1;
                        }
                    }
                }
        
                for (i = 0; i < seth->Nstate; i++) {
                    for (j = 0; j < seth->Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] * ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - sets->gamma_cv[i * seth->Nstate + i]));
                                // printf("sets->force(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
                            }
                        } else {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) 
                                                  * (commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]
                                                    + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] ));
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
                            sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * (sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]));
                        }
                    }
                } else {
                    for (i = 0; i < seth->Nstate; i++) {
                        for (j = 0; j < seth->Nstate; j++) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * (sets->den_e[j * seth->Nstate + i] - sets->gamma_cv[j * seth->Nstate + i]));
                            }
                        }
                    }
                }
            } else if (seth->rep == 1) {
                for (i = 0; i < seth->Nstate; i++) {
                    for (j = 0; j < seth->Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal((sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]);
                            }
                        } else {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal((sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) * (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);
                            }
                        }
                    }
                }
            } else if (seth->rep == 2) {

                for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
                    for (int i = 0; i < seth->Nstate; i++) {
                        for (int j = 0; j < seth->Nstate; j++) {
                            csum1 = 0.0;
                            for (int l = 0; l < seth->Nstate; l++) {
                                csum1 += sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + l * seth->Ndof1 * seth->Ndof2 + k] * sets->V[l * seth->Nstate + j]
                                         - sets->V[i * seth->Nstate + l] * sets->nac[l * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                            }
                            commu_d_V[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] = csum1;
                        }
                    }
                }

                for (i = 0; i < seth->Nstate; i++) {
                    for (j = 0; j < seth->Nstate; j++) {
                        if (i == j) {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] * (sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]));
                                // printf("sets->force(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
                            }
                        } else {
                            for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                                sets->force[k] -= creal((sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) 
                                                  * (commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]
                                                    + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] ));
                                // sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
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

void cal_force_switch(struct set_slave *sets, struct set_host *seth, int para) {
    double c_main[seth->Nstate], sumc_main[seth->Nstate];
    double xe_save[seth->Nstate], pe_save[seth->Nstate];
    double complex gamma_cv_save[seth->Nstate * seth->Nstate], den_e_save[seth->Nstate * seth->Nstate];
    double tempdm1[seth->Nstate * seth->Nstate], tempdm2[seth->Nstate * seth->Nstate], tempdv1[seth->Nstate], tempdv2[seth->Nstate];
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate * seth->Nstate];
    double complex tempcv1[seth->Nstate], tempcv2[seth->Nstate];
    double deltavector[seth->Ndof1 * seth->Ndof2],P_para[seth->Ndof1 * seth->Ndof2],P_ver[seth->Ndof1 * seth->Ndof2];
    double sum;
    double complex csum;
    double complex Q_dia[seth->Nstate * seth->Nstate];
    double prob_hop[seth->Nstate], r_hop;
    double deltaE_mash;
    int id_switch;


    switch (seth->type_hop){
        case 0: // FS-NAF
            for (int i = 0; i < seth->Nstate; i++){
                csum = 0;
                for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    csum += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * sets->P_nuc[j] / sets->mass[j];
                }
                // prob_hop[i] = 2.0 * seth->dt * (sets->xe[sets->id_state] * sets->xe[i] + sets->pe[sets->id_state] * sets->pe[i]) * sum / (sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]);
                prob_hop[i] = 2.0 * seth->dt * creal((sets->xe[sets->id_state] - I * sets->pe[sets->id_state]) * (sets->xe[i] + I * sets->pe[i]) * sum) / (sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]);
                if (prob_hop[i] < 0) prob_hop[i] = 0;
                if (prob_hop[i] > 1) prob_hop[i] = 1;
            }
            
            sum = 0;
            for (int i = 0; i < seth->Nstate; i++) {
                sum += prob_hop[i];
            }
            r_hop = ((double) rand() / RAND_MAX);
            id_switch = sets->id_state;
            for (int i = 0; i < seth->Nstate; i++) {
                r_hop -= prob_hop[i];
                if(r_hop < 0){
                    // double a=0, b=0, c=0;
                    // for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    //     a += 0.5 * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j]
                    //                  / sets->mass[j];
                    //     b += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->P_nuc[j] / sets->mass[j];
                    // }
                    // c = sets->E_adia[i] - sets->E[sets->id_state];
                    id_switch = i;
                    break;
                }
            }

            break;
        
        case 1: // NAF
            if (seth->rep == 0) {
                if (seth->type_evo == 0) {
                    memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                    memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                    // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    for (int i = 0; i < seth->Nstate; i++) {
                        tempcv1[i] = xe_save[i] + I * pe_save[i];
                    }
                    cc_matmul(tempcm1,tempcv1,tempcv2,seth->Nstate,seth->Nstate,1);
                    for (int i = 0; i < seth->Nstate; i++) {
                        sets->xe[i] = creal(tempcv2[i]);
                        sets->pe[i] = cimag(tempcv2[i]);
                    }
                } else if (seth->type_evo == 1) {
                    memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
                    // cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                    cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                }
                memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
                cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
          
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

            id_switch = maxloc(c_main, seth->Nstate);
    }
    
    if (para == 1){
        if (sets->id_state != id_switch) {
            switch (seth->direc_padj) {
                case 0:
                    if (seth->rep == 0) cal_NACV(sets,seth);
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
                    if (seth->rep == 0) cal_NACV(sets,seth);
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        deltavector[i] = creal(sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + id_switch * seth->Ndof1 * seth->Ndof2 + i]);
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

    
    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
    // dd_matmul(sets->U_d2a,Q_dia,tempdm2,seth->Nstate,seth->Nstate,seth->Nstate);
    // dd_matmul(tempdm2,tempdm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);
    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
    cc_matmul(sets->U_d2a,Q_dia,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
    cc_matmul(tempcm2,tempcm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);

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
                            sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]));
                        } else {
                            sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]));
                        }
                    } else if (seth->type_evo == 1) {
                        if (seth->ifscalegamma == 0) {
                            sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        (creal(sets->den_e[i * seth->Nstate + i]) - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]));
                        } else {
                            sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *
                                        (creal(sets->den_e[i * seth->Nstate + i]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - creal(sets->gamma_cv[i * seth->Nstate + i]) + Q_dia[i * seth->Nstate + i]));
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
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            ((sets->xe[j] + I * sets->pe[j]) * (sets->xe[i] - I * sets->pe[i]) * 0.5 - sets->gamma_cv[j * seth->Nstate + i] + Q_dia[j * seth->Nstate + i]));
                            } else {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            ((sets->xe[j] + I * sets->pe[j]) * (sets->xe[i] - I * sets->pe[i]) * 0.5 * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j] + Q_dia[i * seth->Nstate + j]));
                            }
                        } else if (seth->type_evo == 1) {
                            if (seth->ifscalegamma == 0) {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            (sets->den_e[j * seth->Nstate + i] - sets->gamma_cv[j * seth->Nstate + i] + Q_dia[j * seth->Nstate + i]));
                            } else {
                                sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *
                                            (sets->den_e[j * seth->Nstate + i] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[j * seth->Nstate + i] + Q_dia[j * seth->Nstate + i]));
                            }
                        }
                    }
                }
            }
        }
    } else if (seth->rep == 1) {
        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
            sets->force[i] = -creal(sets->dv_adia[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + i]);
        }
        for (int i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                if (i == j) continue;
                if (seth->type_evo == 0) {
                    if (seth->ifscalegamma == 0) {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) - sets->gamma_cv[i * seth->Nstate + j]) *
                                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);
                        }
                    } else {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= creal((0.5 * (sets->xe[i] + I * sets->pe[i]) * (sets->xe[j] - I * sets->pe[j]) * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j]) *
                                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);
                        }
                    }
                } else if (seth->type_evo == 1) {
                    if (seth->ifscalegamma == 0) {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= creal((sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) *
                                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);
                        }
                    } else {
                        for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                            sets->force[k] -= creal((sets->den_e[i * seth->Nstate + j] * (1 + seth->Nstate * seth->gamma_rescale) / (1 + seth->Nstate * seth->gamma_zpe) - sets->gamma_cv[i * seth->Nstate + j]) *
                            (sets->E_adia[j] - sets->E_adia[i]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]);
                        }
                    }
                }
            }
        }
    }

   


}




void cal_force_sh(struct set_slave *sets, struct set_host *seth, int para) {
    double c_main[seth->Nstate], sumc_main[seth->Nstate];
    double xe_save[seth->Nstate], pe_save[seth->Nstate];
    double complex gamma_cv_save[seth->Nstate * seth->Nstate], den_e_save[seth->Nstate * seth->Nstate];
    double tempdm1[seth->Nstate * seth->Nstate], tempdm2[seth->Nstate * seth->Nstate], tempdv1[seth->Nstate], tempdv2[seth->Nstate];
    double complex tempcm1[seth->Nstate * seth->Nstate], tempcm2[seth->Nstate * seth->Nstate];
    double complex tempcv1[seth->Nstate], tempcv2[seth->Nstate];
    double deltavector[seth->Ndof1 * seth->Ndof2],P_para[seth->Ndof1 * seth->Ndof2],P_ver[seth->Ndof1 * seth->Ndof2];
    double sum, x1, x2;
    double complex csum;
    double complex Q_dia[seth->Nstate * seth->Nstate];
    double deltaE_mash;
    double prob_hop[seth->Nstate], r_hop;
    int id_switch;
    double complex commu_d_V[seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2];
    double complex csum1, csum2;
    
    
    switch (seth->type_hop){
        case 0: // FSSH
            for (int i = 0; i < seth->Nstate; i++){
                csum = 0;
                for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    csum += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * sets->P_nuc[j] / sets->mass[j];
                }
                // prob_hop[i] = 2.0 * seth->dt * (sets->xe[sets->id_state] * sets->xe[i] + sets->pe[sets->id_state] * sets->pe[i]) * sum / (sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]);
                if (seth->type_evo == 0) prob_hop[i] = 2.0 * seth->dt * creal((sets->xe[sets->id_state] - I * sets->pe[sets->id_state]) * (sets->xe[i] + I * sets->pe[i]) * csum) / (sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]);
                if (seth->type_evo == 1) prob_hop[i] = 2.0 * seth->dt * creal(sets->den_e[i * seth->Nstate + sets->id_state] * csum) / creal(sets->den_e[sets->id_state * seth->Nstate + sets->id_state]);
                if (prob_hop[i] < 0) prob_hop[i] = 0;
                if (prob_hop[i] > 1) prob_hop[i] = 1;
            }
    
            sum = 0;
            for (int i = 0; i < seth->Nstate; i++) {
                sum += prob_hop[i];
            }
            r_hop = ((double) rand() / RAND_MAX);
            id_switch = sets->id_state;
            for (int i = 0; i < seth->Nstate; i++) {
                r_hop -= prob_hop[i];
                if(r_hop < 0){
                    // double a=0, b=0, c=0;
                    // for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    //     a += 0.5 * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j]
                    //                  / sets->mass[j];
                    //     b += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->P_nuc[j] / sets->mass[j];
                    // }
                    // c = sets->E_adia[i] - sets->E[sets->id_state];
                    id_switch = i;
                    break;
                }
            }

            break;
        
        case 1: // MASH & NAF
            if (seth->rep == 0 || seth->rep == 3) {
                if (seth->type_evo == 0) {
                    memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                    memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                    // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    for (int i = 0; i < seth->Nstate; i++) {
                        tempcv1[i] = xe_save[i] + I * pe_save[i];
                    }
                    cc_matmul(tempcm1,tempcv1,tempcv2,seth->Nstate,seth->Nstate,1);
                    for (int i = 0; i < seth->Nstate; i++) {
                        sets->xe[i] = creal(tempcv2[i]);
                        sets->pe[i] = cimag(tempcv2[i]);
                    }
                } else if (seth->type_evo == 1) {
                    memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
                    // cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                    cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                }
                memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
                cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
          
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

            id_switch = maxloc(c_main, seth->Nstate);
            break;

        case 2: // GFSH
            if (seth->rep == 0) {
                if (seth->type_evo == 0) {
                    memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                    memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                    // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    for (int i = 0; i < seth->Nstate; i++) {
                        tempcv1[i] = xe_save[i] + I * pe_save[i];
                    }
                    cc_matmul(tempcm1,tempcv1,tempcv2,seth->Nstate,seth->Nstate,1);
                    for (int i = 0; i < seth->Nstate; i++) {
                        sets->xe[i] = creal(tempcv2[i]);
                        sets->pe[i] = cimag(tempcv2[i]);
                    }
                } else if (seth->type_evo == 1) {
                    memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
                    // cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                    cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                }
                memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
                cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
          
            }

            if (seth->type_evo == 0){
                x1 = 0.5 * ( sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state])
                    - 0.5 * (sets->xe_old[sets->id_state] * sets->xe_old[sets->id_state] + sets->pe_old[sets->id_state] * sets->pe_old[sets->id_state]);
            } else if (seth->type_evo == 1){
                x1 = creal(sets->den_e[sets->id_state * seth->Nstate + sets->id_state]) - creal(sets->den_e_old[sets->id_state * seth->Nstate + sets->id_state]);
            }
            if (x1 >= 0) {
                memset(prob_hop, 0, seth->Nstate * sizeof(double));
            } else {
                memset(prob_hop, 0, seth->Nstate * sizeof(double));
                sum = 0.0;
                for (int i = 0; i < seth->Nstate; i++) {
                    if (seth->type_evo == 0) {
                        if (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i] 
                            < sets->xe_old[i] * sets->xe_old[i] + sets->pe_old[i] * sets->pe_old[i]) {
                            sum += 0.5 * ( sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i])
                                - 0.5 * (sets->xe_old[i] * sets->xe_old[i] + sets->pe_old[i] * sets->pe_old[i]);
                        }
                    } else if (seth->type_evo == 1) {
                        if (creal(sets->den_e[i * seth->Nstate + i]) < creal(sets->den_e_old[i * seth->Nstate + i])) { 
                            sum += creal(sets->den_e[i * seth->Nstate + i]) - creal(sets->den_e_old[i * seth->Nstate + i]);
                        }
                    }
                }

                

                for (int i = 0; i < seth->Nstate; i++){
                    if (i == sets->id_state) continue;
                    if (seth->type_evo == 0){
                        x2 = 0.5 * ( sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i])
                                - 0.5 * (sets->xe_old[i] * sets->xe_old[i] + sets->pe_old[i] * sets->pe_old[i]);
                    } else if (seth->type_evo == 1){
                        x2 = creal(sets->den_e[i * seth->Nstate + i]) - creal(sets->den_e_old[i * seth->Nstate + i]);
                    }
                    if (x2 > 0) {
                        prob_hop[i] = x2 * x1 / sum;
                        if (seth->type_evo == 0)  prob_hop[i] /= 0.5 * ( sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]);
                        if (seth->type_evo == 1)  prob_hop[i] /= creal(sets->den_e[sets->id_state * seth->Nstate + sets->id_state]);
                    }
                    if (prob_hop[i] < 0) prob_hop[i] = 0;
                    if (prob_hop[i] > 1) prob_hop[i] = 1;
                }
            }


            sum = 0;
            for (int i = 0; i < seth->Nstate; i++) {
                sum += prob_hop[i];
            }
            r_hop = ((double) rand() / RAND_MAX);
            id_switch = sets->id_state;
            for (int i = 0; i < seth->Nstate; i++) {
                r_hop -= prob_hop[i];
                if(r_hop < 0){
                    // double a=0, b=0, c=0;
                    // for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    //     a += 0.5 * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j]
                    //                  / sets->mass[j];
                    //     b += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->P_nuc[j] / sets->mass[j];
                    // }
                    // c = sets->E_adia[i] - sets->E[sets->id_state];
                    id_switch = i;
                    break;
                }
            }

            break;
        
        case 3: // NaF(S)
        
            if (seth->rep == 0) {
                if (seth->type_evo == 0) {
                    memcpy(xe_save,sets->xe,seth->Nstate*sizeof(double));
                    memcpy(pe_save,sets->pe,seth->Nstate*sizeof(double));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dd_matmul(tempdm1,xe_save,sets->xe,seth->Nstate,seth->Nstate,1);
                    // dd_matmul(tempdm1,pe_save,sets->pe,seth->Nstate,seth->Nstate,1);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    for (int i = 0; i < seth->Nstate; i++) {
                        tempcv1[i] = xe_save[i] + I * pe_save[i];
                    }
                    cc_matmul(tempcm1,tempcv1,tempcv2,seth->Nstate,seth->Nstate,1);
                    for (int i = 0; i < seth->Nstate; i++) {
                        sets->xe[i] = creal(tempcv2[i]);
                        sets->pe[i] = cimag(tempcv2[i]);
                    }
                } else if (seth->type_evo == 1) {
                    memcpy(den_e_save,sets->den_e,seth->Nstate * seth->Nstate * sizeof(double complex));
                    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
                    // dc_matmul(tempdm1,den_e_save,tempcm1,seth->Nstate,seth->Nstate,seth->Nstate);
                    // cd_matmul(tempcm1,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
                    cc_matmul(tempcm1,den_e_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                    cc_matmul(tempcm2,sets->U_d2a,sets->den_e,seth->Nstate,seth->Nstate,seth->Nstate);
                }
                memcpy(gamma_cv_save,sets->gamma_cv ,seth->Nstate * seth->Nstate * sizeof(double complex));
                cc_matmul(tempcm1,gamma_cv_save,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
                cc_matmul(tempcm2,sets->U_d2a,sets->gamma_cv,seth->Nstate,seth->Nstate,seth->Nstate);
          
            }

    
            for (int i = 0; i < seth->Nstate; i++){
                
                if (seth->type_evo == 0)  prob_hop[i] = cabs(0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - creal(sets->gamma_cv[i * seth->Nstate + i]));
                if (seth->type_evo == 1)  prob_hop[i] = cabs(creal(sets->den_e[i * seth->Nstate + i]) - creal(sets->gamma_cv[i * seth->Nstate + i]));
                
            }
            


            sum = 0;
            for (int i = 0; i < seth->Nstate; i++) {
                sum += prob_hop[i];
            }
            for (int i = 0; i < seth->Nstate; i++) {
                prob_hop[i] /= sum;
            }
            
            r_hop = ((double) rand() / RAND_MAX);
            id_switch = sets->id_state;
            for (int i = 0; i < seth->Nstate; i++) {
                r_hop -= prob_hop[i];
                if(r_hop < 0){
                    // double a=0, b=0, c=0;
                    // for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    //     a += 0.5 * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j]
                    //                  / sets->mass[j];
                    //     b += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->P_nuc[j] / sets->mass[j];
                    // }
                    // c = sets->E_adia[i] - sets->E[sets->id_state];
                    id_switch = i;
                    break;
                }
            }

            break;


        case 4: // NaF(FS)
            for (int i = 0; i < seth->Nstate; i++){
                csum = 0;
                for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    csum += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * sets->P_nuc[j] / sets->mass[j];
                }
                // prob_hop[i] = 2.0 * seth->dt * (sets->xe[sets->id_state] * sets->xe[i] + sets->pe[sets->id_state] * sets->pe[i]) * sum / (sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]);
                if (seth->type_evo == 0) prob_hop[i] = 2.0 * seth->dt * creal((0.5 * (sets->xe[sets->id_state] - I * sets->pe[sets->id_state]) * (sets->xe[i] + I * sets->pe[i]) - sets->gamma_cv[i * seth->Nstate + sets->id_state]) * csum) 
                                                        / (0.5 * (sets->xe[sets->id_state] * sets->xe[sets->id_state] + sets->pe[sets->id_state] * sets->pe[sets->id_state]) - creal(sets->gamma_cv[sets->id_state * seth->Nstate + sets->id_state]));
                if (seth->type_evo == 1) prob_hop[i] = 2.0 * seth->dt * creal((sets->den_e[i * seth->Nstate + sets->id_state] - sets->gamma_cv[i * seth->Nstate + sets->id_state]) * csum) 
                                                        / creal(sets->den_e[sets->id_state * seth->Nstate + sets->id_state] - sets->gamma_cv[sets->id_state * seth->Nstate + sets->id_state]);
                if (prob_hop[i] < 0) prob_hop[i] = 0;
                if (prob_hop[i] > 1) prob_hop[i] = 1;
            }
    
            sum = 0;
            for (int i = 0; i < seth->Nstate; i++) {
                sum += prob_hop[i];
            }
            r_hop = ((double) rand() / RAND_MAX);
            id_switch = sets->id_state;
            for (int i = 0; i < seth->Nstate; i++) {
                r_hop -= prob_hop[i];
                if(r_hop < 0){
                    // double a=0, b=0, c=0;
                    // for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++){
                    //     a += 0.5 * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j]
                    //                  / sets->mass[j];
                    //     b += sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] 
                    //          * sets->P_nuc[j] / sets->mass[j];
                    // }
                    // c = sets->E_adia[i] - sets->E[sets->id_state];
                    id_switch = i;
                    break;
                }
            }

            break;


    }
    
    if (para == 1){
        if (sets->id_state != id_switch) {
            switch (seth->direc_padj) {
                case 0:
                    if (seth->rep == 0) cal_NACV(sets,seth);
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        sets->P_nuc[i] /= sqrt(sets->mass[i]);
                        deltavector[i] = 0;
                    }
                    for (int i = 0; i < seth->Nstate; i++) {
                        for (int j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
                            if (seth->type_evo == 0) deltavector[j] += 1.0 / sqrt(sets->mass[j]) 
                            * creal(0.5 * (sets->xe[i] - I * sets->pe[i]) * (sets->xe[sets->id_state] + I * sets->pe[sets->id_state]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + j] 
                            - 0.5 * (sets->xe[i] - I * sets->pe[i]) * (sets->xe[id_switch] + I * sets->pe[id_switch]) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + id_switch * seth->Ndof1 * seth->Ndof2 + j]);

                            if (seth->type_evo == 1) deltavector[j] += 1.0 / sqrt(sets->mass[j]) 
                            * creal(sets->den_e[i * seth->Nstate + sets->id_state] * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + j] 
                            - sets->den_e[i * seth->Nstate + id_switch] * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + id_switch * seth->Ndof1 * seth->Ndof2 + j]);

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
                    if (seth->rep == 0) cal_NACV(sets,seth);
                    for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
                        deltavector[i] = creal(sets->nac[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + id_switch * seth->Ndof1 * seth->Ndof2 + i]);
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
                    if (seth->rep == 2) {
                        deltaE_mash += (sets->V[sets->id_state * seth->Nstate + sets->id_state] - sets->V[id_switch * seth->Nstate + id_switch]);
                    } else {
                        deltaE_mash += (sets->E_adia[sets->id_state] - sets->E_adia[id_switch]);
                    }
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
    }

    for (int i = 0; i < seth->Nstate * seth->Nstate; i++) {
        Q_dia[i] = 0;
    }
    Q_dia[sets->id_state * seth->Nstate + sets->id_state] += 1;
    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
    // dd_matmul(sets->U_d2a,Q_dia,tempdm2,seth->Nstate,seth->Nstate,seth->Nstate);
    // dd_matmul(tempdm2,tempdm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);
    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
    cc_matmul(sets->U_d2a,Q_dia,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
    cc_matmul(tempcm2,tempcm1,Q_dia,seth->Nstate,seth->Nstate,seth->Nstate);

    if (seth->rep == 0 || seth->rep == 3) {
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
                    sets->force[j] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] *(Q_dia[i * seth->Nstate + i])); 
                }
            }
        } else {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++){
                        sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * (Q_dia[j * seth->Nstate + i]));
                    }
                }
            }
        }
    } else if (seth->rep == 1) {
        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
            sets->force[i] = - creal(sets->dv_adia[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + i]);
        }
        
    } else if (seth->rep == 2) {
        for (int i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
            sets->force[i] = - creal(sets->dV[sets->id_state * seth->Nstate * seth->Ndof1 * seth->Ndof2 + sets->id_state * seth->Ndof1 * seth->Ndof2 + i]);
        }
    } else if (seth->rep == 3) {

        for (int k = 0; k < seth->Ndof1*seth->Ndof2; k++) {
            for (int i = 0; i < seth->Nstate; i++) {
                for (int j = 0; j < seth->Nstate; j++) {
                    csum1 = 0.0;
                    for (int l = 0; l < seth->Nstate; l++) {
                        csum1 += sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + l * seth->Ndof1 * seth->Ndof2 + k] * sets->V[l * seth->Nstate + j]
                                 - sets->V[i * seth->Nstate + l] * sets->nac[l * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                    }
                    commu_d_V[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] = csum1;
                }
            }
        }

        for (int i = 0; i < seth->Nstate; i++) {
            for (int j = 0; j < seth->Nstate; j++) {
                if (i == j) {
                    for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                        sets->force[k] -= creal(sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] * (Q_dia[j * seth->Nstate + i]));
                        // printf("sets->force(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
                    }
                } else {
                    for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
                        sets->force[k] -= creal((Q_dia[i * seth->Nstate + j])
                                          * (commu_d_V[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k]
                                            + sets->dV[j * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] ));
                        // sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
                    }
                }
            }
        }
        
    }

   


}




void cal_force_eld(struct set_slave *sets,struct set_host *seth) {
    // int i, j, k;
    // double force_trace[seth->Ndof1 * seth->Ndof2];
    // double normforce;

    // // if (seth->if_traceless_force == 1) {
    // //     // for (i = 0; i < seth->Ndof1 * seth->Ndof2; i++) {
    // //     //     sets->force_trace[i] = 0;
    // //     // }
    // //     memset(force_trace,0,seth->Ndof1*seth->Ndof2*sizeof(double));
    // //     for (i = 0; i < seth->Nstate; i++) {
    // //         for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    // //             force_trace[j] += sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] / seth->Nstate;
    // //         }
    // //     }
    // //     for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    // //         for (i = 0; i < seth->Nstate; i++) {
    // //             sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] -= force_trace[j];
    // //         }
    // //         sets->force[j] = -force_trace[j];
    // //     }
    // //     // for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    // //     //     sets->force[j] = -sets->force_trace[j];
    // //     // }
    // // }

    // normforce = 0;
    // for (i = 0; i < seth->Nstate; i++){
    //     if(seth->type_evo == 0) normforce += seth->beta * cexp(-seth->beta * sets->E_adia[i]) * (0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - creal(sets->gamma_cv[i * seth->Nstate + i]));
    //     if(seth->type_evo == 1) normforce += seth->beta * cexp(-seth->beta * sets->E_adia[i]) * creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]);
    // }

     
   

    // switch (seth->type_evo) {
    //     case 0:
    //     case 2:
    //         if (seth->rep == 0) {
    //             // if (seth->calforcetype == 1) {
    //             //     for (i = 0; i < seth->Nstate; i++) {
    //             //         for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    //             //             sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * ((sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + i]));
    //             //         }
    //             //     }
    //             // } else {
    //             //     for (i = 0; i < seth->Nstate; i++) {
    //             //         for (j = 0; j < seth->Nstate; j++) {
    //             //             for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
    //             //                 sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * ((sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) * 0.5 - creal(sets->gamma_cv[i * seth->Nstate + j]));
    //             //             }
    //             //         }
    //             //     }
    //             // }
    //         } else if (seth->rep == 1) {
    //             for (i = 0; i < seth->Nstate; i++) {
    //                 for (j = 0; j < seth->Nstate; j++) {
    //                     if (i == j) {
    //                         for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
    //                             sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[i] + sets->pe[i] * sets->pe[i]) - creal(sets->gamma_cv[i * seth->Nstate + i])) 
    //                                               * seth->beta * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] * cexp(-1 * seth->beta * sets->E_adia[i]) / normforce;
    //                             // printf("sets->force(%d)=%18.8E,%d,%d,%d)=%18.8E\n",i,j,k,l,sets->dv_adia[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k*seth->Ndof2+l]);
    //                         }
    //                     } else {
    //                         for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
    //                             sets->force[k] += (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) 
    //                                               * (cexp( -seth->beta * sets->E_adia[j]) - cexp(-seth->beta * sets->E_adia[i])) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] / normforce;
    //                             // sets->force[k] -= (0.5 * (sets->xe[i] * sets->xe[j] + sets->pe[i] * sets->pe[j]) - creal(sets->gamma_cv[i * seth->Nstate + j])) * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         break;
    //     case 1:
    //         if (seth->rep == 0) {
    //             // if (seth->calforcetype == 1) {
    //             //     for (i = 0; i < seth->Nstate; i++) {
    //             //         for (j = 0; j < seth->Ndof1 * seth->Ndof2; j++) {
    //             //             sets->force[j] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + j] * creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]);
    //             //         }
    //             //     }
    //             // } else {
    //             //     for (i = 0; i < seth->Nstate; i++) {
    //             //         for (j = 0; j < seth->Nstate; j++) {
    //             //             for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
    //             //                 sets->force[k] -= sets->dV[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] * creal(sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]);
    //             //             }
    //             //         }
    //             //     }
    //             // }
    //         } else if (seth->rep == 1) {
    //             for (i = 0; i < seth->Nstate; i++) {
    //                 for (j = 0; j < seth->Nstate; j++) {
    //                     if (i == j) {
    //                         for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
    //                             sets->force[k] -= creal(sets->den_e[i * seth->Nstate + i] - sets->gamma_cv[i * seth->Nstate + i]) 
    //                                             * seth->beta * sets->dv_adia[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof1 * seth->Ndof2 + k] * cexp(-1 * seth->beta * sets->E_adia[i]) / normforce;
    //                         }
    //                     } else {
    //                         for (k = 0; k < seth->Ndof1 * seth->Ndof2; k++) {
    //                             sets->force[k] += creal(sets->den_e[i * seth->Nstate + j] - sets->gamma_cv[i * seth->Nstate + j]) 
    //                                             * (cexp( -seth->beta * sets->E_adia[j]) - cexp(-seth->beta * sets->E_adia[i])) * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] / normforce;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         break;
    // }
}







void cal_NACV(struct set_slave *sets,struct set_host *seth){
    int i, j, imax;
    double dotp, norm_nac1, norm_nac2, cosphase;
    double complex overlap[seth->Nstate * seth->Nstate];
    double overlap2[seth->Nstate * seth->Nstate], vmax;
    double gij[seth->Ndof1 * seth->Ndof2], alpha_BA;
    int id_max[seth->Nstate], idloc, id1, id2;
    double complex eiet[seth->Nstate * seth->Nstate];
    double complex tempcm1[seth->Nstate*seth->Nstate],tempcm2[seth->Nstate*seth->Nstate], tempcm3[seth->Nstate*seth->Nstate], tempcm4[seth->Nstate*seth->Nstate];
    double tempdm1[seth->Nstate*seth->Nstate],tempdm2[seth->Nstate*seth->Nstate],tempdm3[seth->Nstate*seth->Nstate],tempdm4[seth->Nstate*seth->Nstate], tempdv1[seth->Nstate];

    // dia_symmat(seth->Nstate, sets->V, sets->E_adia, sets->U_d2a);
    dia_hermitemat(seth->Nstate, sets->V, sets->E_adia, sets->U_d2a);

    


    if (sets->if_ad_nac) {
        // transpose(sets->U_d2a,tempdm1,seth->Nstate);
        // dd_matmul(tempdm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);
        transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
        cc_matmul(tempcm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);

        memset(overlap2, 0, seth->Nstate * seth->Nstate * sizeof(double));
       
        for (i = 0; i < seth->Nstate * seth->Nstate; i++){
            tempdm1[i] = cabs(overlap[i]);
        }

        for (i = 0; i < seth->Nstate; i++) {
            idloc=maxloc(tempdm1, seth->Nstate * seth->Nstate);
            id1 = idloc / seth->Nstate; 
            id2 = idloc % seth->Nstate; 
            overlap2[idloc] = (creal(overlap[idloc]) >= 0.0) ? 1.0 : -1.0;
            for (j = 0; j < seth->Nstate; j++) {
                tempdm1[id1 * seth->Nstate + j] = 0;
                tempdm1[j * seth->Nstate + id2] = 0;
            }
        }

        cd_matmul(sets->U_d2a, overlap2, tempcm1, seth->Nstate, seth->Nstate, seth->Nstate);
        memcpy(sets->U_d2a,tempcm1,seth->Nstate * seth->Nstate * sizeof(double complex));

        for (i = 0; i < seth->Nstate * seth->Nstate; i++){
            tempdm2[i] = fabs(overlap2[i]);
        }
        dd_matmul(sets->E_adia, tempdm2, tempdv1, 1, seth->Nstate, seth->Nstate);
        memcpy(sets->E_adia, tempdv1, seth->Nstate * sizeof(double));
    }

    if (seth->type_prop_adia > 0 && seth->rep == 1) {
        // transpose(sets->U_d2a,tempdm1,seth->Nstate);
        transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
        cc_matmul(tempcm1, sets->U_ref, sets->overlap_adia, seth->Nstate, seth->Nstate, seth->Nstate);
    }

    memcpy(sets->U_d2a_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double complex));
    memcpy(sets->U_ref, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double complex));
    sets->if_ad_nac = 1;

    // transpose(sets->U_d2a,tempdm1,seth->Nstate);
    transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
    for (i = 0; i < seth->Ndof1; i++) {
        for (j = 0; j < seth->Ndof2; j++) {
            for (int k = 0; k < seth->Nstate * seth->Nstate; k++){
                tempcm2[k] = sets->dV[k * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof2 + j];
            }
            cc_matmul(tempcm1,tempcm2,tempcm3,seth->Nstate,seth->Nstate,seth->Nstate);
            cc_matmul(tempcm3,sets->U_d2a,tempcm4,seth->Nstate,seth->Nstate,seth->Nstate);
            for (int k = 0; k < seth->Nstate * seth->Nstate; k++){
                sets->dv_adia[k * seth->Ndof1 * seth->Ndof2 + i * seth->Ndof2 + j] = tempcm4[k];
            }
        }
    }

    memset(sets->nac, 0, seth->Nstate * seth->Nstate * seth->Ndof1 * seth->Ndof2 * sizeof(double complex));

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


void cal_propagator_adia(int Nstate, double dt, double complex *U, struct set_slave *sets, struct set_host *seth, int para){
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
    double tempdm1[seth->Nstate * seth->Nstate];


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
    transpose_conjugate(C,tempcm1,seth->Nstate);
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

        // switch (seth->type_algorithm) {
        //     case 1:
        //     case 3:
        //         // matmul(eiet, sets->overlap_adia, U);
        //         cd_matmul(eiet,sets->overlap_adia,U,seth->Nstate,seth->Nstate,seth->Nstate);
        //         break;
        //     // case 2:
        //     //     matmul(sets->U_d2a, eiet, U);
        //     //     break;
        //     // case 3:
        //     //     matmul(eiet, sets->overlap_adia, U);
        //     //     break;
        //     // case 4:
        //     //     if (tysets->pe_prop_4cont == 1) {
        //     //         matmul(sets->U_d2a, eiet, U);
        //     //     } else if (tysets->pe_prop_4cont == 2) {
        //     //         matmul(eiet, transpose(sets->U_d2a), U);
        //     //     }
        //     //     break;
        //     // case 5:
        //     //     if(para == 1){
        //     //         // matmul(sets->U_d2a, eiet, U);
        //     //         dc_matmul(sets->overlap_adia,eiet,U,seth->Nstate,seth->Nstate,seth->Nstate);
        //     //     }else if(para == 2){
        //     //         // matmul(eiet, transpose(sets->U_d2a), U);
        //     //         cd_matmul(eiet,sets->U_d2a,U,seth->Nstate,seth->Nstate,seth->Nstate);
        //     //     }
        //     case 6:
        //     case 7:
        //     case 8:
        //     case 9:
                if(para == 1){
                    
                    memcpy(U,eiet,seth->Nstate * seth->Nstate * sizeof(double complex));
                }else if(para == 2){
                    
                    cc_matmul(eiet,sets->overlap_adia,U,seth->Nstate,seth->Nstate,seth->Nstate);
                }else if(para == 3){
                    for(i=0;i<seth->Nstate*seth->Nstate;i++){
                        U[i] = sets->overlap_adia[i];
                    }
                }
        // }
    }

    // int slavecore_id = athread_get_id(-1);
    // if(slavecore_id == 0) printf("U=%18.8E %18.8E %18.8E %18.8E\n",creal(U[0]),creal(U[1]),creal(U[2]),creal(U[3]));
}


void cal_propagator_gen(int Nstate, double complex *H, double dt, double complex *U, struct set_slave *sets, struct set_host *seth){
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
    double tempdm1[seth->Nstate*seth->Nstate],tempdm2[seth->Nstate*seth->Nstate],tempdm3[seth->Nstate*seth->Nstate],tempdm4[seth->Nstate*seth->Nstate], tempdv1[seth->Nstate];
    int id_max[seth->Nstate], idloc, id1, id2;
    double complex overlap[seth->Nstate * seth->Nstate];
    double overlap2[seth->Nstate * seth->Nstate], vmax;

    memset(H_eff,0,seth->Nstate * seth->Nstate * sizeof(double complex));

    for (i = 0; i < seth->Nstate; i++) {
        for (j = 0; j < seth->Nstate; j++) {
            sum=0;
            for (int k=0;k<seth->Ndof1*seth->Ndof2;k++){
                sum += sets->P_nuc[k] * sets->nac[i*seth->Nstate*seth->Ndof1*seth->Ndof2+j*seth->Ndof1*seth->Ndof2+k]/sets->mass[k];
            }
            H_eff[i * seth->Nstate + j] = H[i * seth->Nstate + j] - I * sum;
        }
    }

    dia_hermitemat(seth->Nstate, H_eff, E, C);

  
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
    transpose_conjugate(C,tempcm1,seth->Nstate);
    cc_matmul(eiet,tempcm1,tempcm2,seth->Nstate,seth->Nstate,seth->Nstate);
    cc_matmul(C,tempcm2,U,seth->Nstate,seth->Nstate,seth->Nstate);



    // if (seth->rep == 3) {
    //     // printf("11111\n");
    //     dia_hermitemat(seth->Nstate, H, sets->E_adia, sets->U_d2a);
        
    //     if (sets->if_ad_nac) {
    //         // transpose(sets->U_d2a,tempdm1,seth->Nstate);
    //         // dd_matmul(tempdm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);
    //         transpose_conjugate(sets->U_d2a,tempcm1,seth->Nstate);
    //         cc_matmul(tempcm1,sets->U_ref,overlap,seth->Nstate,seth->Nstate,seth->Nstate);

    //         memset(overlap2, 0, seth->Nstate * seth->Nstate * sizeof(double));
        
    //         for (i = 0; i < seth->Nstate * seth->Nstate; i++){
    //             tempdm1[i] = cabs(overlap[i]);
    //         }

    //         for (i = 0; i < seth->Nstate; i++) {
    //             idloc=maxloc(tempdm1, seth->Nstate * seth->Nstate);
    //             id1 = idloc / seth->Nstate; 
    //             id2 = idloc % seth->Nstate; 
    //             overlap2[idloc] = (creal(overlap[idloc]) >= 0.0) ? 1.0 : -1.0;
    //             for (j = 0; j < seth->Nstate; j++) {
    //                 tempdm1[id1 * seth->Nstate + j] = 0;
    //                 tempdm1[j * seth->Nstate + id2] = 0;
    //             }
    //         }

    //         cd_matmul(sets->U_d2a, overlap2, tempcm1, seth->Nstate, seth->Nstate, seth->Nstate);
    //         memcpy(sets->U_d2a,tempcm1,seth->Nstate * seth->Nstate * sizeof(double complex));

    //         for (i = 0; i < seth->Nstate * seth->Nstate; i++){
    //             tempdm2[i] = fabs(overlap2[i]);
    //         }
    //         dd_matmul(sets->E_adia, tempdm2, tempdv1, 1, seth->Nstate, seth->Nstate);
    //         memcpy(sets->E_adia, tempdv1, seth->Nstate * sizeof(double));
    //     }
    //     memcpy(sets->U_d2a_old, sets->U_ref, seth->Nstate * seth->Nstate * sizeof(double complex));
    //     memcpy(sets->U_ref, sets->U_d2a, seth->Nstate * seth->Nstate * sizeof(double complex));
    //     sets->if_ad_nac = 1;
    // }
}

// void cal_propagator_adia_unsmash(){}
// void cal_propagator_dia_unsmash(){}



void corre_trajprop( struct set_slave *sets, struct set_host *seth){
    double norm_nac, norm_nac_old, dot_nac;

    for (int i = 0; i < seth->Nstate; i++){
        for (int j = 0; j < seth->Nstate; j++){
            if (i == j) continue;
            norm_nac = 0;
            norm_nac_old = 0;
            dot_nac = 0;
            for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++){
                norm_nac += pow(sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k],2);
                norm_nac_old += pow(sets->nac_old[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k], 2);
                dot_nac += sets->nac_old[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k]
                          * sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k];
            }

            if (dot_nac < 0.0) {
                for (int k = 0; k < seth->Ndof1 * seth->Ndof2; k++){
                    sets->nac[i * seth->Nstate * seth->Ndof1 * seth->Ndof2 + j * seth->Ndof1 * seth->Ndof2 + k] *= -1; 
                }
            }

        }
    }
}


void print_traj(FILE *traj_Rnuc, FILE *traj_Pnuc,FILE *traj_ele, FILE *traj_occ, 
                FILE *traj_kinetic, FILE *traj_potential, FILE *traj_nac, FILE *traj_den, FILE *traj_dv, FILE *traj_Ud2a,
                struct set_slave *sets, struct set_host *seth){
        
        fprintf(traj_Rnuc, "%d\n",seth->Natom_mole);
        fprintf(traj_Rnuc, "%s %18.8E %s\n","nuclear coordinate (ang), t=",sets->t_now / seth->unittrans_t, seth->unit_t);
        for (int i = 0; i < seth->Natom_mole; i++){
            fprintf(traj_Rnuc, "%s  %18.8E  %18.8E  %18.8E\n",seth->atomlist_mole[i], sets->R_nuc[i * 3 + 0] * au_2_angstrom, sets->R_nuc[i * 3 + 1] * au_2_angstrom, sets->R_nuc[i * 3 + 2] * au_2_angstrom);
        } 
        fflush(traj_Rnuc);


        fprintf(traj_Pnuc, "%d\n",seth->Natom_mole);
        fprintf(traj_Pnuc, "%s %18.8E %s\n","nuclear momentum (au), t=",sets->t_now / seth->unittrans_t, seth->unit_t);
        for (int i = 0; i < seth->Natom_mole; i++){
            fprintf(traj_Pnuc, "%s  %18.8E  %18.8E  %18.8E\n",seth->atomlist_mole[i], sets->P_nuc[i * 3 + 0], sets->P_nuc[i * 3 + 1], sets->P_nuc[i * 3 + 2]);
        } 
        fflush(traj_Pnuc);

        fprintf(traj_ele, "%18.8E",sets->t_now / seth->unittrans_t);
        if(seth->type_evo == 0){
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(traj_ele, "%18.8E",sets->xe[i]);
            }
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(traj_ele, "%18.8E",sets->pe[i]);
            }
        } else if (seth->type_evo == 1){
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(traj_ele, "%18.8E",creal(sets->den_e[i]));
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(traj_ele, "%18.8E",cimag(sets->den_e[i]));
            }
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            fprintf(traj_ele, "%18.8E",creal(sets->gamma_cv[i]));
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            fprintf(traj_ele, "%18.8E",cimag(sets->gamma_cv[i]));
        }
        fprintf(traj_ele,"\n");
        fflush(traj_ele);

        fprintf(traj_occ, "%18.8E", sets->t_now / seth->unittrans_t);
        fprintf(traj_occ, "  %d", sets->id_state + 1);
        if (seth->ifswitchforce > 0) {
            if (seth->rep == 2){
                fprintf(traj_occ, "  %18.8E",creal(sets->V[sets->id_state * seth->Nstate + sets->id_state]));
            } else if (seth->rep == 1 || seth->rep == 3){
                fprintf(traj_occ, "  %18.8E",sets->E_adia[sets->id_state]);
            }
        }
        
        fprintf(traj_occ,"\n");
        fflush(traj_occ);


        double Ekin = 0;
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            Ekin += 0.5 * sets->P_nuc[i] *  sets->P_nuc[i] / sets->mass[i];
        }
        fprintf(traj_kinetic, "%18.8E  %18.8E\n", sets->t_now / seth->unittrans_t, Ekin);
        fflush(traj_kinetic);

        fprintf(traj_potential, "%18.8E",sets->t_now / seth->unittrans_t);
        if (seth->rep == 1){
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(traj_potential, "%18.8E", sets->E_adia[i]);
            }
        } else {
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(traj_potential, "%18.8E", creal(sets->V[i]));
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(traj_potential, "%18.8E", cimag(sets->V[i]));
            }
        }
        fprintf(traj_potential,"\n");
        fflush(traj_potential);


        
        fprintf(traj_nac, "%s %18.8E %s\n","NAC(au), t=",sets->t_now / seth->unittrans_t, seth->unit_t);
         for (int i = 0; i < seth->Natom_mole * 3; i++){
            for (int j = 0; j < seth->Nstate * seth->Nstate; j++){
                fprintf(traj_nac, "%18.8E", creal(sets->nac[j * seth->Natom_mole * 3 + i]));
            }
            fprintf(traj_nac, "\n");
        }
        fflush(traj_nac);


        fprintf(traj_dv, "%s %18.8E %s\n","dV(au), t=",sets->t_now / seth->unittrans_t, seth->unit_t);
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            for (int j = 0; j < seth->Nstate; j++){
                if (seth->rep == 1){
                    fprintf(traj_dv, "%18.8E", creal(sets->dv_adia[j * seth->Nstate *  seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + i]));
                } else if (seth->rep == 2 || seth->rep == 3){
                    fprintf(traj_dv, "%18.8E", creal(sets->dV[j * seth->Nstate *  seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + i]));
                }
                
            }
            fprintf(traj_dv, "\n");
        }
        fflush(traj_dv);


        if (seth->rep == 2 || seth->rep == 3) {
            fprintf(traj_Ud2a, "%s %18.8E %s\n","Ud2a(au), t=",sets->t_now / seth->unittrans_t, seth->unit_t);
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(traj_Ud2a, "%18.8E", creal(sets->U_d2a[i]));
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(traj_Ud2a, "%18.8E", cimag(sets->U_d2a[i]));
            }
            fprintf(traj_Ud2a,"\n");
            fflush(traj_Ud2a);
        }



        fprintf(traj_den, "%18.8E",sets->t_now / seth->unittrans_t);
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            fprintf(traj_den, "%18.8E", creal(sets->correfun_0 * sets->correfun_t[i]));
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            fprintf(traj_den, "%18.8E", cimag(sets->correfun_0 * sets->correfun_t[i]));
        }
        fprintf(traj_den,"\n");
        fflush(traj_den);



        
}




void print_restart(int itime, int i_re, int igrid, struct set_slave *sets, struct set_host *seth){

        
        FILE *restart = fopen("restart.dat", "w");
        fprintf(restart, "%18.8E\n",sets->t_now);
        fprintf(restart, "%d\n", itime);
        fprintf(restart, "%d\n", i_re);
        fprintf(restart, "%d\n", igrid);
        fprintf(restart, "%18.8E %18.8E\n", creal(sets->correfun_0), cimag(sets->correfun_0));
        fprintf(restart, "%d\n",  sets->id_state);
        fprintf(restart, "%18.8E\n",sets->scale_sqc2);
                
        
        fprintf(restart, "R_nuc:\n");
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            fprintf(restart, "%18.8E",sets->R_nuc[i]);
        } 
        fprintf(restart, "\n");

        fprintf(restart, "P_nuc:\n");
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            fprintf(restart, "%18.8E",sets->P_nuc[i]);
        }
        fprintf(restart, "\n");

        fprintf(restart, "eleDOFs:\n");
        if(seth->type_evo == 0){
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(restart, "%18.8E",sets->xe[i]);
            }
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(restart, "%18.8E",sets->pe[i]);
            }
        } else if (seth->type_evo == 1){
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(restart, "%18.8E",creal(sets->den_e[i]));
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(restart, "%18.8E",cimag(sets->den_e[i]));
            }
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            fprintf(restart, "%18.8E",creal(sets->gamma_cv[i]));
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            fprintf(restart, "%18.8E",cimag(sets->gamma_cv[i]));
        }
        fprintf(restart, "\n");


        fprintf(restart, "potential:\n");
        if (seth->rep == 1){
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(restart, "%18.8E", sets->E_adia[i]);
            }
        } else {
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(restart, "%18.8E", creal(sets->V[i]));
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(restart, "%18.8E", cimag(sets->V[i]));
            }
        }
        fprintf(restart,"\n");
        
        fprintf(restart, "dV:\n");
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            for (int j = 0; j < seth->Nstate; j++){
                if (seth->rep == 1){
                    fprintf(restart, "%18.8E", creal(sets->dv_adia[j * seth->Nstate *  seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + i]));
                } else if (seth->rep == 2 || seth->rep == 3){
                    fprintf(restart, "%18.8E", creal(sets->dV[j * seth->Nstate *  seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + i]));
                }
                
            }
            fprintf(restart, "\n");
        }

        fprintf(restart, "nac:\n");
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            for (int j = 0; j < seth->Nstate * seth->Nstate; j++){
                fprintf(restart, "%18.8E", creal(sets->nac[j * seth->Natom_mole * 3 + i]));
            }
            fprintf(restart, "\n");
        }


        fprintf(restart, "force:\n");
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            fprintf(restart, "%18.8E",sets->force[i]);
        }
        fprintf(restart, "\n");
        

        if (seth->rep == 2 || seth->rep == 3) {
            fprintf(restart, "Eadia:\n");
            for (int i = 0; i < seth->Nstate; i++){
                fprintf(restart, "%18.8E", sets->E_adia[i]);
            }
            fprintf(restart, "\n");
            fprintf(restart, "Ud2a:\n");
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(restart, "%18.8E", creal(sets->U_d2a[i]));
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                fprintf(restart, "%18.8E", cimag(sets->U_d2a[i]));
            }
            fprintf(restart,"\n");
            // fprintf(restart, "Uref:\n");
            // for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            //     fprintf(restart, "%18.8E", creal(sets->U_ref[i]));
            // }
            // for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            //     fprintf(restart, "%18.8E", cimag(sets->U_ref[i]));
            // }
            // fprintf(restart,"\n");
        }





        fflush(restart);
        fclose(restart); 



}




void read_restart(int *itime, int *i_re, int *igrid, struct set_slave *sets, struct set_host *seth){
        double dn1, dn2;
        int i, j;
        char line[512];
        
        FILE *restart = fopen("restart.dat", "r");
       
        if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->t_now = dn1;
        if (fscanf(restart, "%d", &i)   != 1) {}; *itime = i;
        if (fscanf(restart, "%d", &i)   != 1) {}; *i_re = i;
        if (fscanf(restart, "%d", &i)   != 1) {}; *igrid = i;
        if (fscanf(restart, "%lf", &dn1)!= 1) {}; if (fscanf(restart, "%lf", &dn2) != 1) {}; sets->correfun_0 = dn1 + I * dn2;
        if (fscanf(restart, "%d", &i)   != 1) {}; sets->id_state = i;
        if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->scale_sqc2 = dn1;
                
        
        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->R_nuc[i] = dn1;
        } 
        

        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->P_nuc[i] = dn1;
        }
        

        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        if(seth->type_evo == 0){
            for (int i = 0; i < seth->Nstate; i++){
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->xe[i] = dn1;
            }
            for (int i = 0; i < seth->Nstate; i++){
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->pe[i] = dn1;
            }
        } else if (seth->type_evo == 1){
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->den_e[i] = dn1;
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->den_e[i] += I * dn1;
            }
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->gamma_cv[i] = dn1;
        }
        for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
            if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->gamma_cv[i] += I * dn1;
        }
        


        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        if (seth->rep == 1){
            for (int i = 0; i < seth->Nstate; i++){
                // fprintf(restart, "%18.8E", sets->E_adia[i]);
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->E_adia[i] = dn1;
            }
        } else {
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                // fprintf(restart, "%18.8E", creal(sets->V[i]));
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->V[i] = dn1;
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                // fprintf(restart, "%18.8E", cimag(sets->V[i]));
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->V[i] += I * dn1;
            }
        }
        // fprintf(restart,"\n");
        
        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        double dPES_read[seth->Nstate * seth->Natom_mole * 3];
        for (int j = 0; j < seth->Natom_mole * 3; j++) {
            for (int i = 0; i < seth->Nstate; i++) {
                if (fscanf(restart, "%lf", &dPES_read[i * seth->Natom_mole * 3 + j]) != 1) {
                }
            }
        }
        
        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        double nacread[seth->Nstate * seth->Nstate * seth->Natom_mole * 3];
        for (int i = 0; i < seth->Natom_mole * 3; i++) {
            for (int j = 0; j < seth->Nstate * seth->Nstate; j++) {
                if (fscanf(restart, "%lf", &nacread[j * seth->Natom_mole * 3 + i]) != 1) {
                }
            }
        }

        if (seth->rep == 1) {
            memset(sets->dv_adia, 0, seth->Nstate * seth->Nstate * seth->Natom_mole * 3 *  sizeof(double complex));
            memset(sets->nac, 0, seth->Nstate * seth->Nstate * seth->Natom_mole * 3 *  sizeof(double complex));
            for (int i = 0; i < seth->Nstate; i++){
                for (int j = 0; j < seth->Nstate; j++){
                    if (i == j) {
                        
                        for (int k = 0; k < seth->Natom_mole * 3; k++){
                            sets->dv_adia[i * seth->Nstate * seth->Natom_mole * 3 + i * seth->Natom_mole * 3 + k] = dPES_read[i * seth->Natom_mole * 3 + k];
                        }
                    } else {
                        for (int k = 0; k < seth->Natom_mole * 3; k++){
                            sets->nac[i * seth->Nstate * seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + k] = nacread[i * seth->Nstate * seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + k];
                        }
                    }
                }
            }
        } else if (seth->rep == 2 || seth->rep == 3) {
            memset(sets->dV, 0, seth->Nstate * seth->Nstate * seth->Natom_mole * 3 *  sizeof(double complex));
            memset(sets->nac, 0, seth->Nstate * seth->Nstate * seth->Natom_mole * 3 *  sizeof(double complex));
            
            for (int i = 0; i < seth->Nstate; i++){
                for (int j = 0; j < seth->Nstate; j++){
                    if (i == j) {
                        for (int k = 0; k < seth->Natom_mole * 3; k++){
                            sets->dV[i * seth->Nstate * seth->Natom_mole * 3 + i * seth->Natom_mole * 3 + k] = dPES_read[i * seth->Natom_mole * 3 + k];
                        }
                    } else {
                        for (int k = 0; k < seth->Natom_mole * 3; k++){
                            sets->nac[i * seth->Nstate * seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + k] = nacread[i * seth->Nstate * seth->Natom_mole * 3 + j * seth->Natom_mole * 3 + k];
                        }
                    }
                }
            }
        }






        for (int i = 0; i < 2; i++){
            if (fgets(line, sizeof(line), restart) != NULL) {};
        }
        for (int i = 0; i < seth->Natom_mole * 3; i++){
            // fprintf(restart, "%18.8E",sets->force[i]);
            if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->force[i] = dn1;
        }
        // fprintf(restart, "\n");
        

        if (seth->rep == 2 || seth->rep == 3) {
            for (int i = 0; i < 2; i++){
                if (fgets(line, sizeof(line), restart) != NULL) {};
            }
            for (int i = 0; i < seth->Nstate; i++){
                
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->E_adia[i] = dn1;
            }
           
            for (int i = 0; i < 2; i++){
                if (fgets(line, sizeof(line), restart) != NULL) {};
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->U_d2a[i] = dn1; sets->U_ref[i] = dn1;
            }
            for (int i = 0; i < seth->Nstate * seth->Nstate; i++){
                
                if (fscanf(restart, "%lf", &dn1)!= 1) {}; sets->U_d2a[i] += I * dn1; sets->U_ref[i] += I * dn1;
            }
           
        }





        fclose(restart); 



}

void free_vari(struct set_slave *sets, struct set_host *seth) {
    free(sets->R_nuc);
    free(sets->P_nuc);
    free(sets->mass);
    free(sets->force);


    free(sets->R_nuc_init);

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
    } else if (seth->if_allcf >= 2) {
        free(sets->cf0);
        free(sets->cfeff);
        free(sets->weight0);
        free(sets->weightt);
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


    if (seth->mean_nuc == 1) {
        free(sets->R_nuc_mean);
        free(sets->P_nuc_mean);
        free(sets->R2_nuc_mean);
        free(sets->P2_nuc_mean);
    }

    if(seth->if_Pdis == 1){
       free(sets->expisp);
    }
    
    if (seth->ifswitchforce > 0 || seth->type_hop == 1) {
        if (seth->rep == 0) {
            free(sets->U_d2a);
            free(sets->U_ref);
            free(sets->U_d2a_old);
            free(sets->U_ref_old);
            free(sets->E_adia);
            free(sets->E_adia_old);

            if (seth->direc_padj == 0 || seth->direc_padj == 1) {
                free(sets->dv_adia);
                free(sets->nac);
                free(sets->dv_adia_old);
                free(sets->nac_old);
            }
        }
        
    }


}
