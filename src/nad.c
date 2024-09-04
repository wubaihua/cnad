

// #include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// 假设所有使用的模块已经转换为C库或函数
#include "def.h"
// #include "def_sse.h"
// #include "constant.h"
// #include "gmath.h"
// #include "multistate_model.h"

int main(int argc, char *argv[]) {
    // int mpi_size, mpi_rank, mpi_ierr;
    char date_now[20], time_now[20];
    double t1, t2, tt1, tt2;

    // MPI_Init(&argc, &mpi_ierr);
    // MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // if (mpi_rank == 0) {
        printf("====================================================================\n");
        printf("*              NAD: Non-Adiabatic Dynamics simulator               *\n");
        printf("*                        Module of Liouville                       *\n");
        printf("*                                                                  *\n");
        printf("*               Author: Baihua Wu (wubaihua@pku.edu.cn)            *\n");
        printf("*                           version: 1.0                           *\n");
        printf("*                           Aug.8, 2021                            *\n");
        printf("====================================================================\n");

        time_t now = time(NULL);
        struct tm *t = localtime(&now);
        strftime(date_now, sizeof(date_now), "%Y%m%d", t);
        strftime(time_now, sizeof(time_now), "%H%M%S", t);
        printf(" NAD start at %c%c%c%c-%c%c-%c%c %c%c:%c%c:%c%c\n",
               date_now[0], date_now[1], date_now[2], date_now[3],
               date_now[4], date_now[5], date_now[6], date_now[7],
               time_now[0], time_now[1], time_now[2], time_now[3], time_now[4], time_now[5]);
    // }

    // t1 = MPI_Wtime();

    // tt1 = MPI_Wtime();
    initial_para();
    // tt2 = MPI_Wtime();
    // if (mpi_rank == 0) 
    printf("Initialization Finish, using time: %f\n", tt2 - tt1);

    // char filepath[256];
    // if (argc > 1) {
    //     strcpy(filepath, argv[1]);
    // } else {
    //     if (mpi_rank == 0) {
    //         printf("File path not provided.\n");
    //     }
    //     MPI_Finalize();
    //     return 1;
    // }

    // FILE *file = fopen(filepath, "r");
    // if (file == NULL) {
    //     if (mpi_rank == 0) printf("Failed to open file: %s\n", filepath);
    //     MPI_Finalize();
    //     return 1;
    // }

    // if (mpi_rank == 0) printf("Reading files: %s\n", filepath);
    // tt1 = MPI_Wtime();
    // readinp_msmodel(file, &Ndof1, &Ndof2, &Nstate);
    // // 此处需要转换读取 NAMELIST 数据的代码
    // tt2 = MPI_Wtime();
    // if (mpi_rank == 0) printf("Reading input file Finish, using time: %f\n", tt2 - tt1);

    // // 以下代码将依次转化为C语言，省略了大部分详细实现

    // MPI_Barrier(MPI_COMM_WORLD);

    // if (mpi_rank == 0) printf("Simulation Finish\n");

    // if (mpi_rank == 0) printf("Outputting data ..\n");
    // tt1 = MPI_Wtime();

    // if (ifoutputmpi == 1) {
    //     if (mpi_rank == 0) printf("ifoutputmpi=1: output data from each mpi process\n");
    //     fileout_mpi(mpi_rank);
    // }

    // // 继续MPI reduce和数据输出的代码转换

    // t2 = MPI_Wtime();
    // if (mpi_rank == 0) {
    //     printf("Total Running time: %f\n", t2 - t1);
    //     now = time(NULL);
    //     t = localtime(&now);
    //     strftime(date_now, sizeof(date_now), "%Y%m%d", t);
    //     strftime(time_now, sizeof(time_now), "%H%M%S", t);
    //     printf(" NAD end at %c%c%c%c-%c%c-%c%c %c%c:%c%c:%c%c\n",
    //            date_now[0], date_now[1], date_now[2], date_now[3],
    //            date_now[4], date_now[5], date_now[6], date_now[7],
    //            time_now[0], time_now[1], time_now[2], time_now[3], time_now[4], time_now[5]);
    // }

    // MPI_Finalize();
    return 0;
}


