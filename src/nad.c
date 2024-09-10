

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cJSON.h"

// 假设所有使用的模块已经转换为C库或函数
#include "def.h"
// #include "def_sse.h"
#include "constant.h"
#include "gmath.h"
#include "msmodel.h"

int main(int argc, char *argv[]) {
    int mpi_size, mpi_rank, mpi_ierr;
    char date_now[20], time_now[20];
    double t1=0, t2=0, tt1=0, tt2=0;
    time_t now;

    MPI_Init(&argc, &mpi_ierr);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // mpi_rank=0,mpi_size;
    

    if (mpi_rank == 0) {
        printf("====================================================================\n");
        printf("*         cNAD: Non-Adiabatic Dynamics simulator (C-version)       *\n");
        printf("*                        Module of Liouville                       *\n");
        printf("*                                                                  *\n");
        printf("*               Author: Baihua Wu (wubaihua@pku.edu.cn)            *\n");
        printf("*                           version: 1.0                           *\n");
        printf("*                           Sept. 5, 2024                          *\n");
        printf("====================================================================\n");

        now = time(NULL);
        struct tm *t = localtime(&now);
        strftime(date_now, sizeof(date_now), "%Y%m%d", t);
        strftime(time_now, sizeof(time_now), "%H%M%S", t);
        printf(" NAD start at %c%c%c%c-%c%c-%c%c %c%c:%c%c:%c%c\n",
               date_now[0], date_now[1], date_now[2], date_now[3],
               date_now[4], date_now[5], date_now[6], date_now[7],
               time_now[0], time_now[1], time_now[2], time_now[3], time_now[4], time_now[5]);
    }

    t1 = MPI_Wtime();

    tt1 = MPI_Wtime();
    initial_para();
    tt2 = MPI_Wtime();
    if (mpi_rank == 0) 
    printf("Initialization Finish, using time: %f\n", tt2 - tt1);

    // char filepath[256];
    if (argc > 1) {
        filepath=argv[1];
    } else {
        if (mpi_rank == 0) {
            printf("File path not provided.\n");
        }
        MPI_Finalize();
        return 1;
    }
    if (mpi_rank == 0) {
        printf("Reading files: %s\n", filepath);
    }

    readinp();

    tt2 = MPI_Wtime();
    if (mpi_rank == 0) printf("Reading input file Finish, using time: %f\n", tt2 - tt1);



   
    tt1 = MPI_Wtime();
    initial_vari();
    if (mpi_rank == 0) print_info();
   
    init_msmodel(mass);
    
    // if (if_allcf == 2) {
    //     cfweight_msmodel(weight0, weightt, beta);
    // }

   
    tt2 = MPI_Wtime();
    if (mpi_rank == 0) {
        printf("Initializing model Finish, using time: %f\n", tt2 - tt1);
        printf("=====================================================================\n");
        printf("\n");
        printf("<<<<<<<<<<<<<<<<<<<<<<<<---Simulating--->>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
        printf("\n");
    }

    tt1 = MPI_Wtime();
    // init_seed(mpi_rank);
     
    for (int itraj = 1; itraj <= Ntraj; itraj++) {
        if (itraj % mpi_size == mpi_rank) {
            
            sample_msmodel(P_nuc, R_nuc, beta);

            
            
            // if (if_ref == 1) {
            //     for (iref = 1; iref <= Nref; iref++) {
            //         sample_msmodel(&P_nuc_ref[iref * Ndof1 * Ndof2], &R_nuc_ref[iref * Ndof1 * Ndof2], beta);
            //     }
            // }
            sample_ele();

           
           
            evo_traj_new(itraj);
            
        }
    }
    tt2 = MPI_Wtime();
    printf("MPI process %d using time: %f\n", mpi_rank, tt2 - tt1);
    fflush(stdout);  // 强制刷新输出缓冲区


    // // 以下代码将依次转化为C语言，省略了大部分详细实现
    // printf("Process %d before barrier\n", mpi_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("Process %d after barrier\n", mpi_rank);

    if (mpi_rank == 0) printf("Simulation Finish\n");
    fflush(stdout);  // 强制刷新输出缓冲区
    if (mpi_rank == 0) printf("=====================================================================\n");
    fflush(stdout);  // 强制刷新输出缓冲区
    if (mpi_rank == 0) printf("Outputting data ..\n");
    fflush(stdout);  // 强制刷新输出缓冲区
    tt1 = MPI_Wtime();

    // if (ifoutputmpi == 1) {
    //     if (mpi_rank == 0) printf("ifoutputmpi=1: output data from each mpi process\n");
    //     fileout_mpi(mpi_rank);
    // }

    // // 继续MPI reduce和数据输出的代码转换
    mpi_N_nan_sum = (unsigned long long *)malloc(Ngrid * sizeof(unsigned long long));
    memcpy(mpi_N_nan_sum, N_nan_sum, Ngrid * sizeof(unsigned long long));
    MPI_Reduce(N_nan_sum, mpi_N_nan_sum, Ngrid, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) printf("Number of failed trajectories: %d\n", mpi_N_nan_sum[Ngrid-1]);
    MPI_Barrier(MPI_COMM_WORLD);

    if (den != NULL) {
        mpi_den = (double complex *)malloc(Nstate * Nstate * Ngrid * sizeof(double complex));
        real_rho = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        imag_rho = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        mpi_real_den = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        mpi_imag_den = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));

        memcpy(mpi_den, den, Nstate * Nstate * Ngrid * sizeof(double complex));

        // Extract real and imaginary parts
        for (int i = 0; i < Nstate * Nstate * Ngrid; i++) {
            mpi_real_den[i] = creal(mpi_den[i]);
            mpi_imag_den[i] = cimag(mpi_den[i]);
        }

        MPI_Reduce(mpi_real_den, real_rho, Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mpi_imag_den, imag_rho, Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        free(mpi_den);
        free(mpi_real_den);
        free(mpi_imag_den);
    }

    if (population != NULL) {
        mpi_population = (double *)malloc(Nstate * Ngrid * sizeof(double));
        memcpy(mpi_population, population, Nstate * Ngrid * sizeof(double));
        MPI_Reduce(mpi_population, population, Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (if_st_fb == 1) {
            mpi_pop_fb = (double *)malloc(Nstate * Ngrid * 2 * sizeof(double));
            memcpy(mpi_pop_fb, pop_fb, Nstate * Ngrid * 2 * sizeof(double));
            MPI_Reduce(mpi_pop_fb, pop_fb, Nstate * Ngrid * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            free(mpi_pop_fb);
        }
        free(mpi_population);
    }

    // if (cfall != NULL) {
    //     double *mpi_cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    //     double *real_cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    //     double *imag_cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    //     memcpy(mpi_cfall, cfall, Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    //     MPI_Reduce(mpi_cfall, real_cfall, Nstate * Nstate * Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(mpi_cfall, imag_cfall, Nstate * Nstate * Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     free(mpi_cfall);
    //     free(real_cfall);
    //     free(imag_cfall);
    // }

    if (P_nuc_mean != NULL) {
        mpi_P_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        mpi_R_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        mpi_P2_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        mpi_R2_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memcpy(mpi_P_nuc_mean, P_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memcpy(mpi_R_nuc_mean, R_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memcpy(mpi_P2_nuc_mean, P2_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        memcpy(mpi_R2_nuc_mean, R2_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
        MPI_Reduce(mpi_P_nuc_mean, P_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mpi_R_nuc_mean, R_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mpi_P2_nuc_mean, P2_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mpi_R2_nuc_mean, R2_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        free(mpi_P_nuc_mean);
        free(mpi_R_nuc_mean);
        free(mpi_P2_nuc_mean);
        free(mpi_R2_nuc_mean);
    }

    if (cfeff != NULL) {
        mpi_cfeff = (double  complex *)malloc(Ngrid * sizeof(double complex));
        mpi_real_cfeff = (double *)malloc(Ngrid * sizeof(double));
        mpi_imag_cfeff = (double *)malloc(Ngrid * sizeof(double));
        real_cfeff = (double *)malloc(Ngrid * sizeof(double));
        imag_cfeff = (double *)malloc(Ngrid * sizeof(double));
        memcpy(mpi_cfeff, cfeff, Ngrid * sizeof(double complex));
        for (int i = 0; i <  Ngrid; i++) {
            mpi_real_cfeff[i] = creal(mpi_cfeff[i]);
            mpi_imag_cfeff[i] = cimag(mpi_cfeff[i]);
        }

        MPI_Reduce(mpi_real_cfeff, real_cfeff, Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mpi_imag_cfeff, imag_cfeff, Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        free(mpi_cfeff);
        free(mpi_real_cfeff);
        free(mpi_imag_cfeff);
    }

    if (energy_est != NULL) {
        mpi_energy_est = (double *)malloc(Ngrid * sizeof(double));
        memcpy(mpi_energy_est, energy_est, Ngrid * sizeof(double));
        MPI_Reduce(mpi_energy_est, energy_est, Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        free(mpi_energy_est);
    }

    // if (count_st != NULL) {
    //     int *mpi_count_st = (int *)malloc(5 * Ngrid * sizeof(int));
    //     memcpy(mpi_count_st, count_st, 5 * Ngrid * sizeof(int));
    //     MPI_Reduce(mpi_count_st, count_st, 5 * Ngrid, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //     free(mpi_count_st);
    // }

    if (if_Pdis == 1) {
        mpi_expisP = (double complex *)malloc(s_N * sizeof(double complex));
        real_expisP = (double *)malloc(s_N * sizeof(double));
        imag_expisP = (double *)malloc(s_N * sizeof(double));
        mpi_real_expisP = (double *)malloc(s_N * sizeof(double));
        mpi_imag_expisP = (double *)malloc(s_N * sizeof(double));
        memcpy(mpi_expisP, expisP, s_N * sizeof(double  complex));
        for(int i = 0; i < s_N; i++){
            mpi_real_expisP[i] = creal(mpi_expisP[i]);
            mpi_imag_expisP[i] = cimag(mpi_expisP[i]);
        }
        MPI_Reduce(mpi_real_expisP, real_expisP, s_N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(mpi_imag_expisP, imag_expisP, s_N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            for (int i = 0; i < s_N; i++) {
                real_expisP[i] /= Ntraj;
                imag_expisP[i] /= Ntraj;
            }
        }
        free(mpi_expisP);
        free(mpi_real_expisP);
        free(mpi_imag_expisP);
    }

    if (mpi_rank == 0) {
        if (den != NULL) {
            for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
                den[i] = real_rho[i] + I * imag_rho[i];
            }
            if (if_st_nan == 1) {
                for (int igrid = 0; igrid < Ngrid; igrid++) {
                    for (int i = 0; i < Nstate * Nstate; i++) {
                        den[igrid * Nstate * Nstate + i] /= (Ntraj - N_nan_sum[igrid]);
                    }
                }
            } else {
                for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
                    den[i] /= Ntraj;
                }
            }
            
            // if (ifcorreden == 1) {
            //     correct_den();
            // }
        }

        if (population != NULL) {
            if (if_st_nan == 1) {
                for (int igrid = 0; igrid < Ngrid; igrid++) {
                    for (int i = 0; i < Nstate; i++) {
                        population[igrid * Nstate + i] /= (Ntraj - N_nan_sum[igrid]);
                    }
                }
            } else {
                for (int i = 0; i < Ngrid * Nstate; i++) {
                    population[i] /= Ntraj;
                }
            }
            if (if_st_fb == 1) {
                for (int i = 0; i < 2 * Ngrid * Nstate; i++) {
                    pop_fb[i] /= Ntraj;
                }
            }
            
            
        }

        // if (cfall != NULL) {
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         cfall[i] = real_cfall[i] + I * imag_cfall[i];
        //     }
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         cfall[i] /= Ntraj;
        //     }
        // }

        if (cfeff != NULL) {
            for (int i = 0; i < Ngrid; i++) {
                cfeff[i] = (real_cfeff[i] + I * imag_cfeff[i])/Ntraj;
            }
        }

        if (P_nuc_mean != NULL) {
            for (int i = 0; i < Ngrid * Ndof1 * Ndof2; i++) {
                P_nuc_mean[i] /= Ntraj;
                R_nuc_mean[i] /= Ntraj;
                P2_nuc_mean[i] /= Ntraj;
                R2_nuc_mean[i] /= Ntraj;
            }
        }

        if (energy_est != NULL) {
            for (int i = 0; i < Ngrid; i++) {
                energy_est[i] /= Ntraj;
            }
        }

        fileout();
    }

    tt2 = MPI_Wtime();
    if (mpi_rank == 0) {
        printf("Outputting data Finish, using time: %f\n", tt2 - tt1);
        printf("=====================================================================\n");
    }

    

    t2 = MPI_Wtime();
    if (mpi_rank == 0) {
        printf("Total Running time: %f\n", t2 - t1);
        now = time(NULL);
        time_t t = localtime(&now);
        strftime(date_now, sizeof(date_now), "%Y%m%d", t);
        strftime(time_now, sizeof(time_now), "%H%M%S", t);
        printf(" NAD end at %c%c%c%c-%c%c-%c%c %c%c:%c%c:%c%c\n",
               date_now[0], date_now[1], date_now[2], date_now[3],
               date_now[4], date_now[5], date_now[6], date_now[7],
               time_now[0], time_now[1], time_now[2], time_now[3], time_now[4], time_now[5]);
    }

    MPI_Finalize();
    return 0;
}


