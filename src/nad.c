

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
// #include "multistate_model.h"

int main(int argc, char *argv[]) {
    // int mpi_size, mpi_rank, mpi_ierr;
    char date_now[20], time_now[20];
    double t1, t2, tt1, tt2;

    // MPI_Init(&argc, &mpi_ierr);
    // MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    mpi_rank=0;

    // if (mpi_rank == 0) {
    //     printf("====================================================================\n");
    //     printf("*              NAD: Non-Adiabatic Dynamics simulator               *\n");
    //     printf("*                        Module of Liouville                       *\n");
    //     printf("*                                                                  *\n");
    //     printf("*               Author: Baihua Wu (wubaihua@pku.edu.cn)            *\n");
    //     printf("*                           version: 1.0                           *\n");
    //     printf("*                           Aug.8, 2021                            *\n");
    //     printf("====================================================================\n");

    //     time_t now = time(NULL);
    //     struct tm *t = localtime(&now);
    //     strftime(date_now, sizeof(date_now), "%Y%m%d", t);
    //     strftime(time_now, sizeof(time_now), "%H%M%S", t);
    //     printf(" NAD start at %c%c%c%c-%c%c-%c%c %c%c:%c%c:%c%c\n",
    //            date_now[0], date_now[1], date_now[2], date_now[3],
    //            date_now[4], date_now[5], date_now[6], date_now[7],
    //            time_now[0], time_now[1], time_now[2], time_now[3], time_now[4], time_now[5]);
    // // }

    // t1 = MPI_Wtime();

    // tt1 = MPI_Wtime();
    initial_para();
    // tt2 = MPI_Wtime();
    // if (mpi_rank == 0) 
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
        // Ensure msmodelname is the first parameter read
        
        // Read other parameters
       
        // if (cJSON_GetObjectItem(item, "N_bath_SBM") != NULL) {
        //     N_bath_SBM = n_bath_sbm_item->valueint; 
            
        // }
        
        if (NULL !=  cJSON_GetObjectItem(item, "N_bath_SBM")){
            list=cJSON_GetObjectItem(item, "N_bath_SBM");
            N_bath_SBM = list->valueint; 
        }

        if (NULL != cJSON_GetObjectItem(item, "bathtype")) {
        list = cJSON_GetObjectItem(item, "bathtype");
        bathtype = list->valueint; 
        }

        if (NULL != cJSON_GetObjectItem(item, "eps_SBM")) {
            list = cJSON_GetObjectItem(item, "eps_SBM");
            if (list->type == cJSON_Number) {
                eps_SBM = list->valuedouble;
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "delta_SBM")) {
            list = cJSON_GetObjectItem(item, "delta_SBM");
            if (list->type == cJSON_Number) {
                delta_SBM = list->valuedouble; 
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "alpha_SBM")) {
            list = cJSON_GetObjectItem(item, "alpha_SBM");
            if (list->type == cJSON_Number) {
                alpha_SBM = list->valuedouble;
            }
        }

        if (NULL != cJSON_GetObjectItem(item, "omega_c_SBM")) {
            list = cJSON_GetObjectItem(item, "omega_c_SBM");
            if (list->type == cJSON_Number) {
                omega_c_SBM = list->valuedouble; 
            }
        }

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

        // if (NULL != cJSON_GetObjectItem(item, "nproc_sw")) {
        //     list = cJSON_GetObjectItem(item, "nproc_sw");
        //     if (list->type == cJSON_Number) {
        //         nproc_sw = list->valueint; 
        //     }
        // }

        item = item->next;
    }

    cJSON_Delete(json);
    printf("msmodelname: %s\n",msmodelname);

    printf("N_bath_SBM: %d\n", N_bath_SBM);
    printf("bathtype: %d\n", bathtype);
    printf("eps_SBM: %f\n", eps_SBM);
    printf("delta_SBM: %f\n", delta_SBM);
    printf("alpha_SBM: %f\n", alpha_SBM);
    printf("omega_c_SBM: %f\n", omega_c_SBM);


    printf("beta: %f\n", beta);
    printf("Ntraj: %d\n", Ntraj);
    printf("dt: %f\n", dt);
    printf("ttot: %f\n", ttot);
    printf("Nbreak: %d\n", Nbreak);
    printf("method: %s\n", method);
    printf("init_occ: %d\n", init_occ);
    printf("rep: %d\n", rep);
    printf("outputtype: %d\n", outputtype);
    printf("forcetype: %d\n", forcetype);
    printf("calforcetype: %d\n", calforcetype);


    // FILE *file = fopen("sbminp.json", "r");

    
    // // if (file == NULL) {
    // //     if (mpi_rank == 0) printf("Failed to open file: %s\n", filepath);
    // //     // MPI_Finalize();
    // //     return 1;
    // // }
    // // printf(filepath);
    // // printf(file);
    
    // // if (mpi_rank == 0) printf("Reading files: %s\n", filepath);
    // // /;
    
    // // tt1 = MPI_Wtime();
    
    // // // 此处需要转换读取 NAMELIST 数据的代码
    
    // fseek(file, 0, SEEK_END);
    // long length = ftell(file);
    // fseek(file, 0, SEEK_SET);

    // char *content = malloc(length + 1);
    // if (content == NULL) {
    //     if (mpi_rank == 0) printf("Failed to allocate memory\n");
    //     fclose(file);
    //     return 1;
    // }

    // fread(content, 1, length, file);
    // content[length] = '\0';
    // fclose(file);

    

    // cJSON *json = cJSON_Parse(content);
    // if (json == NULL) {
    //     if (mpi_rank == 0) printf("Failed to parse JSON\n");
    //     free(content);
    //     return 1;
    // }

    

    // cJSON *item = cJSON_GetObjectItem(json, "msmodelname");
    // if (item) {
    //     strcpy(msmodelname, item->valuestring);
    // } else {
    //     if (mpi_rank == 0) {
    //         printf("No msmodelname!\n");
    //     }
    //     cJSON_Delete(json);
    //     free(content);
    //     return 1;
    // }
    // cJSON_Delete(json);

    // char *line = strtok(content, "\n");
    // printf("debug677\n");
    // const char *ptr = content;
    //    while ((ptr = strchr(ptr, '{')) != NULL) {
    //     json = cJSON_Parse(ptr);
    //     if (!json) {
    //         printf("Error parsing JSON: %s\n", cJSON_GetErrorPtr());
    //         free((void*)content);
    //         return 1;
    //     }
    //     readinp_msmodel(json, Ndof1, Ndof2, Nstate);
    //     readinp_para(json);
    //     cJSON_Delete(json);
    //     ptr++;
    // }

    // free((void*)content);
    // printf("debug666\n");

    // tt2 = MPI_Wtime();
    if (mpi_rank == 0) printf("Reading input file Finish, using time: %f\n", tt2 - tt1);

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


