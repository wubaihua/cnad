

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "constant.h"
#include "gmath.h"  // Assuming you have a math.h header for functions like box_muller
// #include "constant.h"
#include <string.h>
#include <stdio.h>
// #include "cJSON.h"
#include "msmodelio.h"
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif




// Initialize model parameters
void parameter_mole(double *mass, struct set_host *setm) {




    char *atom_names[] = {"Bq","H ","He",  // Bq, 1-~2  
    "Li","Be","B ","C ","N ","O ","F ","Ne", // 3~10
    "Na","Mg","Al","Si","P ","S ","Cl","Ar", // !11~18
    "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", //  !19~36
    "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe",//  !37~54
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", //  !55~71
    "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", //  !72~86
    "Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", //  !87~103
    "Rf","Db","Sg","Bh","Hs","Mt" };
    

    double atom_masses[] = { 0.0,	1.00794,	4.002602, // & !Bq; 1-2
        6.941, 9.0121831, 10.811, 12.0107, 14.0067, 15.9994, 18.998403163, 20.1797,// & !3-10
        22.98976928, 24.3050, 26.9815385, 28.0855, 30.973761998, 32.065, 35.453, 39.948, //  !11-18
        39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798, //  !19-36
        85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, 97, 101.07, 102.90550, 106.42, 107.8682, 112.414, 114.818, 118.710, 121.760, 127.60, 126.90, 131.29, // !37-54
        132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 145, 150.36, 151.964, 157.25, 158.92535, 162.500, 164.93033, 167.259, 168.93422, 173.05, 174.967,  // !55-71
        178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.592, 204.38, 207.2, 208.98040, 209, 210, 222,  // !72-86
        223, 226, 227, 232.03806, 231.03588, 238.02891, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0,  257.0, 258.0, 259.0, 262.0, // !87-103
        267.0, 268.0, 269.0, 270.0, 269.0, 277.0 // !104-109
    };

    for (int i = 0; i < 110; i++) {
        setm->atommass_mole[i] = atom_masses[i];
        strcpy(setm->atomname_mole[i], atom_names[i]);
    }




    for (int i = 0; i < setm->Natom_mole; i++){
        for (int j = 0; j < 110; j++){
            if (strcmp(setm->atomname_mole[j], setm->atomlist_mole[i]) == 0){
                setm->atomindexlist_mole[i] = j;
                break;
            }
        }
    }

    for (int i = 0; i < setm->Natom_mole; i++){
        // printf("%s\n",setm->atomlist_mole[i]);
        // printf("%d\n",setm->atomindexlist_mole[i]);
        mass[i * 3 + 0] = setm->atommass_mole[setm->atomindexlist_mole[i]] * amu_2_au;
        mass[i * 3 + 1] = setm->atommass_mole[setm->atomindexlist_mole[i]] * amu_2_au;
        mass[i * 3 + 2] = setm->atommass_mole[setm->atomindexlist_mole[i]] * amu_2_au;
        // printf("mass=%f\n",mass[i * 3 + 0]/amu_2_au);
    }


    



    
}

// Sample the initial conditionals for trajectories of the model
void sample_mole(double *P, double *R, double beta, struct set_host *setm) {
    
    // double x2;
    // for (int j = 0; j < setm->N_bath_mole; j++) {
    //     if (beta > 99999) {
    //         if(setm->if_classical == 0){
    //             box_muller(&P[j], &x2, sqrt(0.5 * hbar * setm->omega_mole[j]), 0.0);
    //             box_muller(&R[j], &x2, sqrt(0.5 * hbar / setm->omega_mole[j]), 0.0);
    //         } else {
    //             P[j] = 0;
    //             R[j] = 0;
    //         }
            
    //     } else {
    //         if(setm->if_classical == 0){
    //             box_muller(&P[j], &x2, sqrt(0.5 * hbar * setm->omega_mole[j] / tanh(0.5 * beta * hbar * setm->omega_mole[j])), 0.0);
    //             box_muller(&R[j], &x2, sqrt(0.5 * hbar / (tanh(0.5 * beta * hbar * setm->omega_mole[j]) * setm->omega_mole[j])), 0.0);
    //         // P[j]=1,R[j]=1;   //debug 
    //         } else {
    //             box_muller(&P[j], &x2, sqrt(1.0 / beta), 0.0);
    //             box_muller(&R[j], &x2, sqrt(1.0 / (beta * setm->omega_mole[j] * setm->omega_mole[j])), 0.0);
    //         }
        
    //     }
    // }


    for (int i = 0; i < 3 * setm->Natom_mole; i++){
        R[i] = setm->R0_nuc_mole[i];
        P[i] = setm->P0_nuc_mole[i];
    }

    // exit(-1);
}

// Build the diabatic potential matrix of the model
void V_mole(double *R, double *H, int forcetype, struct set_host *setm) {
    // double sum1=0, sum2=0;
    // for(int i=0;i<setm->N_bath_mole;i++){
    //     sum1 += setm->c_mole[i] * R[i];
    //     sum2 += setm->omega_mole[i] * setm->omega_mole[i] * R[i] * R[i];
    // }
    // switch (forcetype) {
    //     case 0:
    //         H[0] = setm->eps_mole + sum1 + 0.5 * sum2;
    //         H[1] = setm->delta_mole;
    //         H[3] = -setm->eps_mole - sum1 + 0.5 * sum2;
    //         H[2] = setm->delta_mole;
    //         break;
    //     case 1:
    //         H[0] = setm->eps_mole + sum1;
    //         H[1] = setm->delta_mole;
    //         H[3] = -setm->eps_mole - sum1;
    //         H[2] = setm->delta_mole;
    //         break;
    //     default:
    //         fprintf(stderr, "ERROR: unknown forcetype\n");
    //         exit(1);
    // }

    // int slavecore_id=athread_get_id(-1);
//     if(slavecore_id == 0){
//         printf("c=%18.8E %18.8E w=%18.8E %18.8E\n", setm->c_mole[0],setm->c_mole[1],setm->omega_mole[0],setm->omega_mole[1]);
//         printf("R=%18.8E %18.8E V=%18.8E %18.8E %18.8E %18.8E\n", R[0],R[1],H[0],H[1],H[2],H[3]);
       
//     }

}

// Build the first-order derivative matrix of the model
void dV_mole(double *R, double complex *dH, int forcetype, struct set_host *setm) {
    FILE *inp = fopen("bdfinp", "w");
    fprintf(inp, "%s\n", "[GEOM]");
    fprintf(inp, "%s%d\n", "xyz=\"\"\"",setm->Natom_mole);
    fprintf(inp, "%s\n", "Comment: generated by CNAD");
    for (int i = 0; i < setm->Natom_mole; i++){
        fprintf(inp, "%s  %f  %f  %f\n",setm->atomlist_mole[i], R[i * 3 + 0] * au_2_angstrom, R[i * 3 + 1] * au_2_angstrom, R[i * 3 + 2] * au_2_angstrom);
    }
    
    fprintf(inp, "\"\"\"\n");
    fprintf(inp, "%s\n", setm->qmkeyword_mole);
    
    fclose(inp);

    // char cmd[256];
    // int ret;
    // snprintf(cmd, sizeof(cmd),"mkdir ~/scratch");
    // // 执行命令
    // ret = system(cmd);
    
    // snprintf(cmd, sizeof(cmd),
    //     "python /mnt/f/GitHub/psinad_wbh/scripts/QM.py -d run_bdf -i bdfinp -qm \"bdf\" -t 0 > QMlog");
    // // 执行命令
    // ret = system(cmd);


    FILE *qmlog= fopen("QMlog","r");
    if (!qmlog) {
        printf("Error: Cannot open QMlog\n");
        exit(1);
    }

    char line[512];
    char last_line[512] = {0};
    while (fgets(line, sizeof(line), qmlog)) {
        strcpy(last_line, line); // 保存最后一行
    }
    fclose(qmlog);

    if (strstr(last_line, "--- calculation terminated normally ---") == NULL) {
        printf("Error: QM finishes unnormally!\n");
        exit(-1);
    };// else {
    //    printf("yayyy\n");// debug
    //}

    FILE *qmdata= fopen("run_bdf/interface.ds","r");
    if (!qmdata) {
        printf("Error: Cannot open interface.ds\n");
        exit(1);
    }
    for (int i = 0; i < 6; i++){
        if (fgets(line, sizeof(line), qmdata) != NULL) {
            // printf("%s\n",line);
        }
    }

    double PES[setm->Nstate];
    for (int i = 0; i < setm->Nstate; i++) {
        if (fscanf(qmdata, "%lf", &PES[i]) != 1) {
            
        }
        printf("%f\n",PES[i]);
    }

    for (int i = 0; i < 4; i++){
        if (fgets(line, sizeof(line), qmdata) != NULL) {
            printf("%s\n",line);
        }
    }

    double dPES[setm->Nstate * setm->Natom_mole * 3];
    for (int j = 0; j < setm->Natom_mole * 3; j++) {
        for (int i = 0; i < setm->Nstate; i++) {
            if (fscanf(qmdata, "%lf", &dPES[i * setm->Natom_mole * 3 + j]) != 1) {
                  
            }
            printf("%18.8E  ",dPES[i * setm->Natom_mole * 3 + j]); 
        }
        
        printf("\n");
    }



    

    exit(-1);
}

// Calculate the nuclear force of the model
void nucforce_mole(double *R, double *nf, struct set_host *setm) {
    // for (int j = 0; j < setm->N_bath_mole; j++) {
    //     nf[j] = setm->omega_mole[j]* setm->omega_mole[j] * R[j];
    // }
    
}

