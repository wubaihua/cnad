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
    #include <slave.h>
#endif



// void init_slave(){
//     initial_vari();
//     init_msmodel(mass);
// }



void dynamics_slave(struct set_host *seth){
    struct set_slave sets;

    int slavecore_id, itraj;

    int run_size = seth->Ntraj/seth->mpi_size;

  
    #ifdef sunway
    slavecore_id=athread_get_id(-1);
    if (slavecore_id < seth->nproc_sw){
        initial_vari(&sets,seth);
        init_msmodel(sets.mass,seth);
    }
    #elif defined(x86)
    initial_vari(&sets,seth);
    init_msmodel(sets.mass,seth);
    #endif
    

    
    #ifdef sunway
    init_seed(seth->mpi_rank*seth->nproc_sw+slavecore_id);
    for (int itraj = 1; itraj <= run_size; itraj++) {
        if (itraj % seth->nproc_sw == slavecore_id) { 
            sample_msmodel(sets.P_nuc, sets.R_nuc, seth->beta,seth);
            sample_ele(&sets,seth);       
            evo_traj_new(itraj,&sets,seth);
        }
    }
    athread_ssync_array();
    #elif defined(x86)
    init_seed(seth->mpi_rank);
    for (int itraj = 1; itraj <= run_size; itraj++) {
        sample_msmodel(sets.P_nuc, sets.R_nuc, seth->beta,seth);
        sample_ele(&sets,seth);       
        evo_traj_new(itraj,&sets,seth);
    }
    #endif

    
    // printf("%d %18.8E %18.8E %18.8E\n",slavecore_id,sets.population[0 * seth->Ngrid  + seth->Ngrid -2]/run_size*seth->nproc_sw,
    //                                                 sets.population[1 * seth->Ngrid  + seth->Ngrid -2]/run_size*seth->nproc_sw,
    //                                                 sets.population[2 * seth->Ngrid  + seth->Ngrid -2]/run_size*seth->nproc_sw); // debug

    // printf("%d %18.8E %18.8E %18.8E\n",slavecore_id,sets.population[0 * seth->Ngrid  + seth->Ngrid -1]/run_size*seth->nproc_sw,
    //                                                 sets.population[1 * seth->Ngrid  + seth->Ngrid -1]/run_size*seth->nproc_sw,
    //                                                 sets.population[2 * seth->Ngrid  + seth->Ngrid -1]/run_size*seth->nproc_sw); // debug

// printf("%d %18.8E %18.8E\n",slavecore_id,creal(sets.den[0 * seth->Ngrid * seth->Nstate + 0 *seth->Ngrid  + 0]),creal(sets.den[1 * seth->Ngrid * seth->Nstate + 1 *seth->Ngrid + 0])); // debug

    //  for (int i=0;i<seth->Ngrid;i++){
    //         printf("%d %18.8E %18.8E\n",i,seth->mpi_population[0*seth->Ngrid + i],seth->mpi_population[1*seth->Ngrid + i]);
    // }
    

    // for (int i = 0; i < Ngrid; i++){
    //     mpi_N_nan_sum[i] += N_nan_sum[i];
    // }

    // if (den != NULL) {
        
    //     for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
    //         mpi_den[i] += den[i];
    //     }
    // }

    // if (population != NULL) {
    //     for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
    //         mpi_population[i] += population[i];
    //     }
    // }

#ifdef sunway
    athread_ssync_array();
    for(int idcore = 0; idcore < 64; idcore++){
        if(slavecore_id == idcore && slavecore_id < seth->nproc_sw){
            for (int i = 0; i < seth->Ngrid; i++){
                seth->save_N_nan_sum[i*seth->nproc_sw+idcore] += sets.N_nan_sum[i];
            }

            if (seth->outputtype >= 0) {
                
                for (int i = 0; i < seth->Nstate * seth->Nstate * seth->Ngrid; i++){
                    // mpi_den[i] += den[i];
                    // seth->mpi_real_den[i] += creal(sets.den[i]);
                    // seth->mpi_imag_den[i] += cimag(sets.den[i]);
                    seth->save_real_den[i*seth->nproc_sw+idcore] = creal(sets.den[i]);
                    seth->save_imag_den[i*seth->nproc_sw+idcore] = cimag(sets.den[i]);
                }
            }
            if (seth->outputtype != 0) {
                for (int i = 0; i < seth->Nstate * seth->Ngrid; i++){
                    // seth->mpi_population[i] += sets.population[i];
                    seth->save_population[i*seth->nproc_sw+idcore] = sets.population[i];
                    // printf("%d %18.8E %18.8E\n",i,seth->mpi_population[i], sets.population[i]);
                    // seth->mpi_population[i] += 1;
                }
                // for (int i = 0; i < seth->Nstate; i++){
                //     for (int j=0; j< seth->Ngrid; j++){
                //         seth->mpi_population[i* seth->Ngrid + j] += creal(sets.den[i* seth->Nstate * seth->Ngrid + i * seth->Ngrid + j]);
                //     }
                //     // seth->mpi_population[i] += 1;
                // }
                if(seth->if_st_fb == 1){
                    for (int i = 0; i < seth->Nstate * seth->Ngrid * 2; i++){                 
                        seth->save_pop_fb[i*seth->nproc_sw+idcore] = sets.pop_fb[i];
                    }
                }

            }
            // if(id==0){
            //     memcpy(fi_time_grid,timegrid,Ngrid*sizeof(double));
            // }


            // for (int i = 0; i < Ngrid; i++){
            //     fi_time_grid[i] = timegrid[i];
            // }

        }
        athread_ssync_array();
    }
    if(slavecore_id < seth->nproc_sw){
        free_vari(&sets,seth);
    }
    athread_ssync_array();
   
#elif defined(x86)
    
    for (int i = 0; i < seth->Ngrid; i++){
        seth->mpi_N_nan_sum[i] = sets.N_nan_sum[i];
    }
    if (seth->outputtype >= 0) {
        
        for (int i = 0; i < seth->Nstate * seth->Nstate * seth->Ngrid; i++){
            // mpi_den[i] += den[i];
            seth->mpi_real_den[i] = creal(sets.den[i]);
            seth->mpi_imag_den[i] = cimag(sets.den[i]);
            // seth->save_real_den[i*64+idcore] = creal(sets.den[i]);
            // seth->save_imag_den[i*64+idcore] = cimag(sets.den[i]);
        }
    }
    if (seth->outputtype != 0) {
        for (int i = 0; i < seth->Nstate * seth->Ngrid; i++){
            seth->mpi_population[i] = sets.population[i];
            // seth->save_population[i*64+idcore] = sets.population[i];
            // printf("%d %18.8E %18.8E\n",i,seth->mpi_population[i], sets.population[i]);
            // seth->mpi_population[i] += 1;
        }
        // for (int i = 0; i < seth->Nstate; i++){
        //     for (int j=0; j< seth->Ngrid; j++){
        //         seth->mpi_population[i* seth->Ngrid + j] += creal(sets.den[i* seth->Nstate * seth->Ngrid + i * seth->Ngrid + j]);
        //     }
        //     // seth->mpi_population[i] += 1;
        // }
        if (seth->if_st_fb == 1){
            for (int i = 0; i < seth->Nstate * seth->Ngrid * 2; i++){
                seth->mpi_pop_fb[i] = sets.pop_fb[i];
            }
        }
    }
    // if(id==0){
    //     memcpy(fi_time_grid,timegrid,Ngrid*sizeof(double));
    // }
    
    // for (int i = 0; i < Ngrid; i++){
    //     fi_time_grid[i] = timegrid[i];
    // }

    free_vari(&sets,seth);
// exit(-1);
#endif

    // printf("%d %18.8E %18.8E %18.8E\n",slavecore_id,seth->save_population[0 * seth->Ngrid *64  + (seth->Ngrid -1)*64+slavecore_id]/run_size*64,
    //                                          seth->save_population[1 * seth->Ngrid *64  + (seth->Ngrid -1)*64+slavecore_id]/run_size*64,
    //                                          seth->save_population[2 * seth->Ngrid *64  + (seth->Ngrid -1)*64+slavecore_id]/run_size*64); // debug
    // printf("%d %18.8E %18.8E %18.8E\n",slavecore_id,seth->mpi_population[0 * seth->Ngrid  + seth->Ngrid -1]/run_size,
    //                                                 seth->mpi_population[1 * seth->Ngrid  + seth->Ngrid -1]/run_size,
    //                                                 seth->mpi_population[2 * seth->Ngrid  + seth->Ngrid -1]/run_size); // debug

   
    // if (seth->outputtype >= 0) {
                
    //     for (int i = 0; i < seth->Nstate * seth->Nstate * seth->Ngrid; i++){
            
    //         seth->save_real_den[i*64+slavecore_id] = creal(sets.den[i]);
    //         seth->save_imag_den[i*64+slavecore_id] = cimag(sets.den[i]);
    //     }
    // }

    
//     for (int j=0;j<64;j++){
//     printf("%d %18.8E %18.8E\n",j,seth->save_real_den[0 * seth->Ngrid * seth->Nstate + 0 *seth->Ngrid  + 0 + j],seth->save_imag_den[1 * seth->Ngrid * seth->Nstate + 1 *seth->Ngrid + j]); // debug

// }


    
    
    
   

}



// void data_transport(int id){

//     int slavecore_id;
   
//     slavecore_id=athread_get_id(-1);


   

//     if(slavecore_id == id){
//         for (int i = 0; i < Ngrid; i++){
//             mpi_N_nan_sum[i] += N_nan_sum[i];
//         }

//         if (den != NULL) {
            
//             for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
//                 // mpi_den[i] += den[i];
//                 mpi_real_den[i] += creal(den[i]);
//                 mpi_imag_den[i] += cimag(den[i]);

//             }
//         }

//         if (population != NULL) {
//             for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
//                 mpi_population[i] += population[i];
//             }
//         }





//         // if(id==0){
//         //     memcpy(fi_time_grid,timegrid,Ngrid*sizeof(double));
//         // }


//         // for (int i = 0; i < Ngrid; i++){
//         //     fi_time_grid[i] = timegrid[i];
//         // }

        
//     }

    

// }



// void free_slave(){
//     free(population);
//     free(den);
//     free(R_nuc);
//     free(P_nuc);
//     free(dv_adia);
//     free(nac);
    
// }