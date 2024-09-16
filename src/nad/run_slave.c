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
#include <slave.h>



void init_slave(){
    initial_vari();
    init_msmodel(mass);
}



void dynamics_slave(int run_size){

    int slavecore_id, itraj;
   
    slavecore_id=athread_get_id(-1);


    init_seed(mpi_rank*64+slavecore_id);

     
    for (int itraj = 1; itraj <= run_size; itraj++) {
        if (itraj % 64 == slavecore_id) {
            
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


    // printf("%f\n",population[0 * Ngrid  + Ngrid -1 ]); // debug




    

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
    
   

}



void data_transport(int id){

    int slavecore_id;
   
    slavecore_id=athread_get_id(-1);


   

    if(slavecore_id == id){
        for (int i = 0; i < Ngrid; i++){
            mpi_N_nan_sum[i] += N_nan_sum[i];
        }

        if (den != NULL) {
            
            for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
                // mpi_den[i] += den[i];
                mpi_real_den[i] += creal(den[i]);
                mpi_imag_den[i] += cimag(den[i]);

            }
        }

        if (population != NULL) {
            for (int i = 0; i < Nstate * Nstate * Ngrid; i++){
                mpi_population[i] += population[i];
            }
        }





        // if(id==0){
        //     memcpy(fi_time_grid,timegrid,Ngrid*sizeof(double));
        // }


        // for (int i = 0; i < Ngrid; i++){
        //     fi_time_grid[i] = timegrid[i];
        // }

        
    }

    

}