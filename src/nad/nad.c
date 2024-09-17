

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cJSON.h"
#include <memory.h>
#include <complex.h>
// 假设所有使用的模块已经转换为C库或函数
// #include "def.h"
#include "def_host.h"
// #include "def_sse.h"
#include "constant.h"
#include "gmath.h"
#include "msmodel.h"
#include "msmodelio.h"
#include "run_slave.h"
#include <athread.h>


void unsignlong_reduce(unsigned long long *sendbuf, unsigned long long *recvbuf, int count, MPI_Op op, int root, MPI_Comm comm);
void double_reduce(double *sendbuf, double *recvbuf, int count, MPI_Op op, int root, MPI_Comm comm);

int main(int argc, char *argv[]) {
    // int mpi_size, mpi_rank, mpi_ierr;
    char date_now[20], time_now[20];
    double t1=0, t2=0, tt1=0, tt2=0;
    time_t now;
    double dsum;
    unsigned long long lsum;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // mpi_rank=0,mpi_size;
    
    athread_init();

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


    // initial_vari();
    // init_msmodel(mass);
    athread_spawn(init_slave,0);
    athread_join();

    init_host();

    if (mpi_rank == 0) print_info();
    
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

    int run_size = Ntraj/mpi_size;
    athread_spawn(dynamics_slave,run_size);
    athread_join();

     
    // for (int itraj = 1; itraj <= Ntraj; itraj++) {
    //     if (itraj % mpi_size == mpi_rank) {
            
    //         sample_msmodel(P_nuc, R_nuc, beta);

            
            
    //         // if (if_ref == 1) {
    //         //     for (iref = 1; iref <= Nref; iref++) {
    //         //         sample_msmodel(&P_nuc_ref[iref * Ndof1 * Ndof2], &R_nuc_ref[iref * Ndof1 * Ndof2], beta);
    //         //     }
    //         // }
    //         sample_ele();

           
           
    //         evo_traj_new(itraj);
            
    //     }
    // }
    tt2 = MPI_Wtime();
    printf("MPI process %d using time: %f\n", mpi_rank, tt2 - tt1);
    // fflush(stdout);  // 强制刷新输出缓冲区


    // // 以下代码将依次转化为C语言，省略了大部分详细实现
    // printf("Process %d before barrier\n", mpi_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("Process %d after barrier\n", mpi_rank);

    if (mpi_rank == 0) printf("Simulation Finish\n");
    // fflush(stdout);  // 强制刷新输出缓冲区
    if (mpi_rank == 0) printf("=====================================================================\n");
    // fflush(stdout);  // 强制刷新输出缓冲区
    if (mpi_rank == 0) printf("Outputting data ..\n");
    // fflush(stdout);  // 强制刷新输出缓冲区
    tt1 = MPI_Wtime();

    for (int i = 0; i < 64; i++){
        athread_spawn(data_transport,i);
        athread_join();
    }

    
    // athread_spawn(free_slave,0);
    // athread_join();

    for (int i = 0; i < Ngrid; i++){
        fi_time_grid[i] = i*dt*Nbreak;
    }
   
    // printf("1111\n");
    // athread_spawn(data_transport,1);
    // athread_join();
   
    // printf("2222\n");
    // printf("pop=%18.8E\n",mpi_population[0 * Ngrid  + Ngrid -1 ]); // debug

    // printf("den=%18.8E\n",mpi_real_den[0 * Ngrid*Nstate + 0* Ngrid + Ngrid -1 ]); // debug


    if (ifoutputmpi == 1) {
        if (mpi_rank == 0) printf("ifoutputmpi=1: output data from each mpi process\n");
        fileout_mpi(mpi_rank);
    }
   
    // 继续MPI reduce和数据输出的代码转换
    // printf("1111\n");
    // fi_N_nan_sum = (unsigned long long *)malloc(Ngrid * sizeof(unsigned long long));
    // printf("2222\n");
    // MPI_Reduce(mpi_N_nan_sum, mpi_N_nan_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    for (int i = 0; i < Ngrid; i++){
        MPI_Reduce(&mpi_N_nan_sum[i], &mpi_N_nan_sum[i], 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // MPI_Status status;
    // int ierr;
    // ierr = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    // printf("%d\n",ierr);
    // if (ierr != MPI_SUCCESS) {
    //     printf("Error in MPI_Probe\n");
    // }

    // MPI_Request request;
    // MPI_Ireduce(mpi_N_nan_sum, mpi_N_nan_sum, Ngrid, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD,&request);
    // printf("debugr44444%d\n",mpi_rank);
    // MPI_Wait(&request, MPI_STATUS_IGNORE);
    // // MPI_Barrier(MPI_COMM_WORLD);
    // // fflush(stdout);  // 强制刷新输出缓冲区
    // // int bufsize = 300000; // 缓冲区大小
    // // void *buffer = malloc(bufsize);
    // // MPI_Buffer_attach(buffer, bufsize);
    // // unsigned long long *result_N_nan_sum = NULL;
    // if (mpi_rank == 0) {
    //     result_N_nan_sum = (unsigned long long *)malloc(Ngrid * sizeof(unsigned long long));
    // }
    // unsignlong_reduce(mpi_N_nan_sum, result_N_nan_sum, Ngrid, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (mpi_rank == 0) memcpy(mpi_N_nan_sum, result_N_nan_sum,Ngrid * sizeof(unsigned long long));
    // free(result_N_nan_sum);
    // MPI_Buffer_detach(&buffer, &bufsize);
    // free(buffer);

    

    if (mpi_population != NULL) {
        // fi_population = (double *)malloc(Nstate * Ngrid * sizeof(double));
        // MPI_Reduce(mpi_population, mpi_population, Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        // for (int i = 0; i < Nstate * Ngrid; i++){
        //     MPI_Reduce(&mpi_population[i], &dsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //     if(mpi_rank == 0) mpi_population[i] = dsum;
        // }
        // printf("111111111111111111\n");
      
        for (int i = 0; i < Nstate * Ngrid; i++) {
            // printf("aaa: %d, %d, %f, %f \n", mpi_rank, i, mpi_population[i], dsum);
            MPI_Reduce(&mpi_population[i], &mpi_population[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }


        // double *result_population = NULL;
        // if (mpi_rank == 0) {
        //     result_population = (double *)malloc(Nstate * Ngrid * sizeof(double));
        // }
        // double_reduce(mpi_population, result_population, Nstate * Ngrid, MPI_SUM, 0, MPI_COMM_WORLD);
        // if (mpi_rank == 0) memcpy(mpi_population, result_population,Nstate * Ngrid * sizeof(double));
        // free(result_population);

        if (if_st_fb == 1) {
            // fi_pop_fb = (double *)malloc(Nstate * Ngrid * 2 * sizeof(double));
            
            // MPI_Reduce(mpi_pop_fb, mpi_pop_fb, Nstate * Ngrid * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            for (int i = 0; i < Nstate * Ngrid * 2; i++){
                MPI_Reduce(&mpi_pop_fb[i], &mpi_pop_fb[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            // free(mpi_pop_fb);
        }
        // free(mpi_population);
    }

    // if (mpi_den != NULL) {
    if (mpi_real_den != NULL) {
        // printf("test111111\n");
        // fi_den = (double complex *)malloc(Nstate * Nstate * Ngrid * sizeof(double complex));
        // printf("test2222\n");
        // real_rho = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        // printf("test3\n");
        // imag_rho = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        // printf("test4\n");
        // mpi_real_den = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));
        // // printf("test5\n");
        // mpi_imag_den = (double *)malloc(Nstate * Nstate * Ngrid * sizeof(double));

        // printf("tes6666\n");
        // Extract real and imaginary parts
        // for (int i = 0; i < Nstate * Nstate * Ngrid; i++) {
        //     mpi_real_den[i] = creal(mpi_den[i]);
        //     mpi_imag_den[i] = cimag(mpi_den[i]);
        // }

        // printf("mpi_real_den=%p\n",(void*)mpi_real_den);
        // printf("mpi_imag_den=%p\n",(void*)mpi_imag_den);
        // printf("real_rho=%p\n",(void*)real_rho);
        // printf("imag_rho=%p\n",(void*)imag_rho);

        // printf("test7777771\n");
        // MPI_Barrier(MPI_COMM_WORLD);
        // MPI_Reduce(mpi_real_den, mpi_real_den, Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        // MPI_Reduce(mpi_imag_den, mpi_imag_den, Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
            MPI_Reduce(&mpi_real_den[i], &mpi_real_den[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&mpi_imag_den[i], &mpi_imag_den[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // // printf("test9999999\n");
        // // free(mpi_den);
        // // printf("test10\n");
        // free(mpi_real_den);
        // // printf("test11\n");
        // free(mpi_imag_den);
        // // printf("test12s\n");
    }

    // if(mpi_rank==0)printf("den=%18.8E\n",mpi_real_den[0 * Ngrid*Nstate + 0* Ngrid + Ngrid -1 ]); // debug
    

    

    // // if (cfall != NULL) {
    // //     double *mpi_cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    // //     double *real_cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    // //     double *imag_cfall = (double *)malloc(Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    // //     memcpy(mpi_cfall, cfall, Nstate * Nstate * Nstate * Nstate * Ngrid * sizeof(double));
    // //     MPI_Reduce(mpi_cfall, real_cfall, Nstate * Nstate * Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // //     MPI_Reduce(mpi_cfall, imag_cfall, Nstate * Nstate * Nstate * Nstate * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // //     free(mpi_cfall);
    // //     free(real_cfall);
    // //     free(imag_cfall);
    // // }

    // if (P_nuc_mean != NULL) {
    //     mpi_P_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     mpi_R_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     mpi_P2_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     mpi_R2_nuc_mean = (double *)malloc(Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     memcpy(mpi_P_nuc_mean, P_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     memcpy(mpi_R_nuc_mean, R_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     memcpy(mpi_P2_nuc_mean, P2_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     memcpy(mpi_R2_nuc_mean, R2_nuc_mean, Ndof1 * Ndof2 * Ngrid * sizeof(double));
    //     MPI_Reduce(mpi_P_nuc_mean, P_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(mpi_R_nuc_mean, R_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(mpi_P2_nuc_mean, P2_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(mpi_R2_nuc_mean, R2_nuc_mean, Ndof1 * Ndof2 * Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     free(mpi_P_nuc_mean);
    //     free(mpi_R_nuc_mean);
    //     free(mpi_P2_nuc_mean);
    //     free(mpi_R2_nuc_mean);
    // }

    // if (cfeff != NULL) {
    //     mpi_cfeff = (double  complex *)malloc(Ngrid * sizeof(double complex));
    //     mpi_real_cfeff = (double *)malloc(Ngrid * sizeof(double));
    //     mpi_imag_cfeff = (double *)malloc(Ngrid * sizeof(double));
    //     real_cfeff = (double *)malloc(Ngrid * sizeof(double));
    //     imag_cfeff = (double *)malloc(Ngrid * sizeof(double));
    //     memcpy(mpi_cfeff, cfeff, Ngrid * sizeof(double complex));
    //     for (int i = 0; i <  Ngrid; i++) {
    //         mpi_real_cfeff[i] = creal(mpi_cfeff[i]);
    //         mpi_imag_cfeff[i] = cimag(mpi_cfeff[i]);
    //     }

    //     MPI_Reduce(mpi_real_cfeff, real_cfeff, Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(mpi_imag_cfeff, imag_cfeff, Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     free(mpi_cfeff);
    //     free(mpi_real_cfeff);
    //     free(mpi_imag_cfeff);
    // }

    // if (energy_est != NULL) {
    //     mpi_energy_est = (double *)malloc(Ngrid * sizeof(double));
    //     memcpy(mpi_energy_est, energy_est, Ngrid * sizeof(double));
    //     MPI_Reduce(mpi_energy_est, energy_est, Ngrid, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     free(mpi_energy_est);
    // }

    // // if (count_st != NULL) {
    // //     int *mpi_count_st = (int *)malloc(5 * Ngrid * sizeof(int));
    // //     memcpy(mpi_count_st, count_st, 5 * Ngrid * sizeof(int));
    // //     MPI_Reduce(mpi_count_st, count_st, 5 * Ngrid, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    // //     free(mpi_count_st);
    // // }

    // if (if_Pdis == 1) {
    //     mpi_expisP = (double complex *)malloc(s_N * sizeof(double complex));
    //     real_expisP = (double *)malloc(s_N * sizeof(double));
    //     imag_expisP = (double *)malloc(s_N * sizeof(double));
    //     mpi_real_expisP = (double *)malloc(s_N * sizeof(double));
    //     mpi_imag_expisP = (double *)malloc(s_N * sizeof(double));
    //     memcpy(mpi_expisP, expisP, s_N * sizeof(double  complex));
    //     for(int i = 0; i < s_N; i++){
    //         mpi_real_expisP[i] = creal(mpi_expisP[i]);
    //         mpi_imag_expisP[i] = cimag(mpi_expisP[i]);
    //     }
    //     MPI_Reduce(mpi_real_expisP, real_expisP, s_N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(mpi_imag_expisP, imag_expisP, s_N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //     if (mpi_rank == 0) {
    //         for (int i = 0; i < s_N; i++) {
    //             real_expisP[i] /= Ntraj;
    //             imag_expisP[i] /= Ntraj;
    //         }
    //     }
    //     free(mpi_expisP);
    //     free(mpi_real_expisP);
    //     free(mpi_imag_expisP);
    // }

    if (mpi_rank == 0) {
        if (mpi_real_den != NULL) {
            // for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
            //     fi_den[i] = real_rho[i] + I * imag_rho[i];
            // }
            if (if_st_nan == 1) {
                for (int igrid = 0; igrid < Ngrid; igrid++) {
                    for (int i = 0; i < Nstate * Nstate; i++) {
                        mpi_real_den[i * Ngrid + igrid] /= (Ntraj - fi_N_nan_sum[igrid]);
                        mpi_imag_den[i * Ngrid + igrid] /= (Ntraj - fi_N_nan_sum[igrid]);
                    }
                }
            } else {
                for (int i = 0; i < Ngrid * Nstate * Nstate; i++) {
                    mpi_real_den[i] /= Ntraj;
                    mpi_imag_den[i] /= Ntraj;
                }
            }
            
            // if (ifcorreden == 1) {
            //     correct_den();
            // }
        }

        if (mpi_population != NULL) {
            if (if_st_nan == 1) {
                for (int igrid = 0; igrid < Ngrid; igrid++) {
                    for (int i = 0; i < Nstate; i++) {
                        mpi_population[i * Ngrid + igrid] /= (Ntraj - fi_N_nan_sum[igrid]);
                    }
                }
            } else {
                for (int i = 0; i < Ngrid * Nstate; i++) {
                    mpi_population[i] /= Ntraj;
                }
            }
            if (if_st_fb == 1) {
                for (int i = 0; i < 2 * Ngrid * Nstate; i++) {
                    mpi_pop_fb[i] /= Ntraj;
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


    MPI_Barrier(MPI_COMM_WORLD);

    athread_halt();
    exit(-1);
    // MPI_Comm_free(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}








void unsignlong_reduce(unsigned long long *sendbuf, unsigned long long *recvbuf, int count, MPI_Op op, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        // Initialize recvbuf with sendbuf values
        for (int i = 0; i < count; i++) {
            recvbuf[i] = sendbuf[i];
        }

        // Receive data from other processes and perform reduction
        for (int i = 1; i < size; i++) {
            unsigned long long *tempbuf = (unsigned long long *)malloc(count * sizeof(unsigned long long));
            MPI_Recv(tempbuf, count, MPI_UNSIGNED_LONG_LONG, i, 0, comm, MPI_STATUS_IGNORE);
            for (int j = 0; j < count; j++) {
                recvbuf[j] += tempbuf[j]; // Assuming MPI_SUM operation
            }
            free(tempbuf);
        }
    } else {
        // Send data to root process
        MPI_Bsend(sendbuf, count, MPI_UNSIGNED_LONG_LONG, root, 0, comm);
    }
}






void double_reduce(double *sendbuf, double *recvbuf, int count, MPI_Op op, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        // Initialize recvbuf with sendbuf values
        for (int i = 0; i < count; i++) {
            recvbuf[i] = sendbuf[i];
        }

        // Receive data from other processes and perform reduction
        for (int i = 1; i < size; i++) {
            double *tempbuf = (double *)malloc(count * sizeof(double));
            MPI_Recv(tempbuf, count, MPI_DOUBLE, i, 0, comm, MPI_STATUS_IGNORE);
            for (int j = 0; j < count; j++) {
                recvbuf[j] += tempbuf[j]; // Assuming MPI_SUM operation
            }
            free(tempbuf);
        }
    } else {
        // Send data to root process
        MPI_Send(sendbuf, count, MPI_DOUBLE, root, 0, comm);
    }
}


