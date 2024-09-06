

// #include <mpi.h>
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

    // MPI_Init(&argc, &mpi_ierr);
    // MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    mpi_rank=0,mpi_size;
    

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

    // t1 = MPI_Wtime();

    // tt1 = MPI_Wtime();
    initial_para();
    // tt2 = MPI_Wtime();
    if (mpi_rank == 0) 
    // printf("Initialization Finish, using time: %f\n", tt2 - tt1);

    // char filepath[256];
    if (argc > 1) {
        filepath=argv[1];
    } else {
        if (mpi_rank == 0) {
            printf("File path not provided.\n");
        }
        // MPI_Finalize();
        return 1;
    }
    printf("Reading files: %s\n", filepath);

    readinp();

    // tt2 = MPI_Wtime();
    // if (mpi_rank == 0) printf("Reading input file Finish, using time: %f\n", tt2 - tt1);



   
    // tt1 = MPI_Wtime();
    initial_vari();
    if (mpi_rank == 0) print_info();
   
    init_msmodel(mass);
    
    // if (if_allcf == 2) {
    //     cfweight_msmodel(weight0, weightt, beta);
    // }

    
    // tt2 = MPI_Wtime();
    if (mpi_rank == 0) {
        printf("Initializing model Finish, using time: %f\n", tt2 - tt1);
        printf("=====================================================================\n");
        printf("\n");
        printf("<<<<<<<<<<<<<<<<<<<<<<<<---Simulating--->>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
        printf("\n");
    }

    // tt1 = MPI_Wtime();
    // init_seed(mpi_rank);
    
    for (int itraj = 1; itraj <= Ntraj; itraj++) {
        // if (itraj % mpi_size == mpi_rank) {
            
            sample_msmodel(P_nuc, R_nuc, beta);

            
            
            // if (if_ref == 1) {
            //     for (iref = 1; iref <= Nref; iref++) {
            //         sample_msmodel(&P_nuc_ref[iref * Ndof1 * Ndof2], &R_nuc_ref[iref * Ndof1 * Ndof2], beta);
            //     }
            // }
            sample_ele();

           
           
            evo_traj_new(itraj);
            
        // }
    }
    // tt2 = MPI_Wtime();
    // printf("MPI process %d using time: %f\n", mpi_rank, tt2 - tt1);
    


    // // 以下代码将依次转化为C语言，省略了大部分详细实现

    // MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) printf("Simulation Finish\n");

    if (mpi_rank == 0) printf("Outputting data ..\n");
    // tt1 = MPI_Wtime();

    // if (ifoutputmpi == 1) {
    //     if (mpi_rank == 0) printf("ifoutputmpi=1: output data from each mpi process\n");
    //     fileout_mpi(mpi_rank);
    // }

    // // 继续MPI reduce和数据输出的代码转换

    if (mpi_rank == 0) {
        // if (den != NULL) {
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         den[i] = real_rho[i] + I * imag_rho[i];
        //     }
        //     if (if_st_nan == 1) {
        //         for (int igrid = 0; igrid < Ngrid; igrid++) {
        //             for (int i = 0; i < Nstate * Nstate; i++) {
        //                 den[igrid * Nstate * Nstate + i] /= (Ntraj - N_nan_sum[igrid]);
        //             }
        //         }
        //     } else {
        //         for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //             den[i] /= Ntraj;
        //         }
        //     }
        //     if (strcmp(method, "sqc") == 0 || strcmp(method, "SQC") == 0) {
        //         // Add specific method handling here
        //     }
        //     if (ifcorreden == 1) {
        //         correct_den();
        //     }
        // }

        if (population != NULL) {
            // if (if_st_nan == 1) {
            //     for (int igrid = 0; igrid < Ngrid; igrid++) {
            //         for (int i = 0; i < Nstate; i++) {
            //             population[igrid * Nstate + i] /= (Ntraj - N_nan_sum[igrid]);
            //         }
            //     }
            // } else {
                for (int i = 0; i < Ngrid * Nstate; i++) {
                    population[i] /= Ntraj;
                }
            // }
            // if (if_st_fb == 1) {
            //     for (int i = 0; i < Nstate; i++) {
            //         pop_fb[i] /= Ntraj;
            //     }
            // }
            
            
        }

        // if (cfall != NULL) {
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         cfall[i] = real_cfall[i] + I * imag_cfall[i];
        //     }
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         cfall[i] /= Ntraj;
        //     }
        // }

        // if (cfeff != NULL) {
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         cfeff[i] = real_cfeff[i] + I * imag_cfeff[i];
        //     }
        //     for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
        //         cfeff[i] /= Ntraj;
        //     }
        // }

        // if (P_nuc_mean != NULL) {
        //     for (int i = 0; i < Ngrid * Nstate; i++) {
        //         P_nuc_mean[i] /= Ntraj;
        //         R_nuc_mean[i] /= Ntraj;
        //         P2_nuc_mean[i] /= Ntraj;
        //         R2_nuc_mean[i] /= Ntraj;
        //     }
        // }

        // if (energy_est != NULL) {
        //     for (int i = 0; i < Ngrid; i++) {
        //         energy_est[i] /= Ntraj;
        //     }
        // }

        fileout();
    }

    // double tt2 = MPI_Wtime();
    if (mpi_rank == 0) {
        printf("Outputting data Finish, using time: %f\n", tt2 - tt1);
        printf("=====================================================================\n");
    }

    // double t2 = MPI_Wtime();

    // t2 = MPI_Wtime();
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

    // MPI_Finalize();
    return 0;
}


