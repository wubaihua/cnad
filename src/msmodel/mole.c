

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
#ifdef x86
    #include "def.h"
#endif
// #ifdef sunway
//     #include <slave.h>
//     #include <athread.h>
// #endif






// Initialize model parameters
void parameter_mole(double *mass, struct set_host *setm) {

    FILE *qmpath= fopen(setm->path_qmpath_mole, "rb");
    if (!qmpath) {
        perror("qmpath_mole opening failed");
        exit(1);
    }
    
    if (fgets(setm->qmsoft_mole, sizeof(setm->qmsoft_mole), qmpath) != NULL) {
        // printf("%s\n",line);
    }
    if (fgets(setm->qmscript_mole, sizeof(setm->qmscript_mole), qmpath) != NULL) {
        // printf("%s\n",line);
    }
    if (fgets(setm->qmscratch_mole, sizeof(setm->qmscratch_mole), qmpath) != NULL) {
        // printf("%s\n",line);
    }

    setm->qmscript_mole[strcspn(setm->qmscript_mole, "\r\n")] = 0;
    setm->qmsoft_mole[strcspn(setm->qmsoft_mole, "\r\n")] = 0;
    setm->qmscratch_mole[strcspn(setm->qmscratch_mole, "\r\n")] = 0;
    
    // printf("%s\n",setm->qmsoft_mole);
    // printf("%s\n",setm->qmscript_mole);
    // printf("%s\n",setm->qmscratch_mole);
    // exit(-1);



    // if (strcmp(setm->msmodelname, "mole") == 0) {
    FILE *qmkeywordfile = fopen(setm->path_qmkeyword_mole, "rb");
    if (!qmkeywordfile) {
        perror("qmkeyword_mole opening failed");
        // return NULL;
    }
    char buffer[2000];  // 读取缓冲区
    size_t total_length = 0;
    setm->qmkeyword_mole = malloc(1);  // 动态字符串初始化
    setm->qmkeyword_mole[0] = '\0';
    while (fgets(buffer, sizeof(buffer), qmkeywordfile)) {
        total_length += strlen(buffer);
        setm->qmkeyword_mole = realloc(setm->qmkeyword_mole, total_length + 1);
        strcat(setm->qmkeyword_mole, buffer);  // 追加到动态字符串
    }
    // printf("=============================\n");
    // printf("%s", setm->qmkeyword_mole);
    fclose(qmkeywordfile);
    setm->R0_nuc_mole = (double *)malloc(setm->Natom_mole * 3 * sizeof(double));
    setm->P0_nuc_mole = (double *)malloc(setm->Natom_mole * 3 * sizeof(double));
    setm->atomlist_mole = malloc(setm->Natom_mole * sizeof(char*));
    setm->atomindexlist_mole = malloc(setm->Natom_mole * sizeof(int));
    FILE *R0_nuc_file = fopen(setm->path_R0_nuc_mole, "rb");
    if (!R0_nuc_file) {
        perror("R0_nuc_mole opening failed");
        // return NULL;
    }
    int natom_check;
    int istat = fscanf(R0_nuc_file, "%d\n", &natom_check);  // 读取原子数
    // printf("natomcheck= %d \n",natom_check);
    char comment_line[256];
    if (!fgets(comment_line, sizeof(comment_line), R0_nuc_file)) {
        fprintf(stderr, "警告: 读取注释行失败\n");
    }
    char c2[3];
    double x, y, z;
    for (int i = 0; i < setm->Natom_mole; i++) {
        int read_count = fscanf(R0_nuc_file, "%2s %lf %lf %lf", c2, &x, &y, &z);
        // printf("read_count: %d\n", read_count);
        if (read_count != 4) {
            fprintf(stderr, "Error reading R0_nuc_mole file at line %d\n", i + 1);
            exit(EXIT_FAILURE);
        }
        // printf("c2: %s, x: %lf, y: %lf, z: %lf\n", c2, x, y, z);
        c2[2] = '\0';
        if (strlen(c2) == 1) {
            c2[1] = ' ';
            c2[2] = '\0';
        }
        setm->atomlist_mole[i] = strdup(c2);
        setm->R0_nuc_mole[i * 3] = x / au_2_angstrom;
        setm->R0_nuc_mole[i * 3 + 1] = y / au_2_angstrom;
        setm->R0_nuc_mole[i * 3 + 2] = z / au_2_angstrom;
    }
    // for (int i = 0; i < setm->Natom_mole; i++) {
    //     printf("%s %18.8E %18.8E %18.8E\n", setm->atomlist_mole[i], setm->R0_nuc_mole[i * 3], setm->R0_nuc_mole[i * 3 + 1], setm->R0_nuc_mole[i * 3 + 2]);
    // }
    fclose(R0_nuc_file);
    FILE *P0_nuc_file = fopen(setm->path_P0_nuc_mole, "rb");
    if (!P0_nuc_file) {
        perror("P0_nuc_mole opening failed");
        // return NULL;
    }

    istat = fscanf(P0_nuc_file, "%d\n", &natom_check);  // 读取原子数
    // printf("natomcheck= %d \n",natom_check);
    // char comment_line[256];
    if (!fgets(comment_line, sizeof(comment_line), P0_nuc_file)) {
        fprintf(stderr, "警告: 读取注释行失败\n");
    }
    for (int i = 0; i < setm->Natom_mole; i++) {
        int read_count = fscanf(P0_nuc_file, "%2s %lf %lf %lf", c2, &x, &y, &z);
        // printf("read_count: %d\n", read_count);
        if (read_count != 4) {
            fprintf(stderr, "Error reading R0_nuc_mole file at line %d\n", i + 1);
            exit(EXIT_FAILURE);
        }
        // printf("c2: %s, x: %lf, y: %lf, z: %lf\n", c2, x, y, z);
        c2[2] = '\0';
        // setm->atomlist_mole[i] = strdup(c2);
        setm->P0_nuc_mole[i * 3] = x;
        setm->P0_nuc_mole[i * 3 + 1] = y;
        setm->P0_nuc_mole[i * 3 + 2] = z;
    }
    // for (int i = 0; i < setm->Natom_mole; i++) {
    //     printf("%s %18.8E %18.8E %18.8E\n", setm->atomlist_mole[i], setm->P0_nuc_mole[i * 3], setm->P0_nuc_mole[i * 3 + 1], setm->P0_nuc_mole[i * 3 + 2]);
    // }
    fclose(P0_nuc_file);



        // exit(0);
    // }


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
    // FILE *inp = fopen("bdfinp", "w");
    // fprintf(inp, "%s\n", "[GEOM]");
    // fprintf(inp, "%s%d\n", "xyz=\"\"\"",setm->Natom_mole);
    // fprintf(inp, "%s\n", "Comment: generated by CNAD");
    // for (int i = 0; i < setm->Natom_mole; i++){
    //     fprintf(inp, "%s  %f  %f  %f\n",setm->atomlist_mole[i], R[i * 3 + 0] * au_2_angstrom, R[i * 3 + 1] * au_2_angstrom, R[i * 3 + 2] * au_2_angstrom);
    // }
    
    // fprintf(inp, "\"\"\"\n");
    // fprintf(inp, "%s\n", setm->qmkeyword_mole);
    
    // fclose(inp);

    // // char cmd[256];
    // // int ret;
    // // snprintf(cmd, sizeof(cmd),"mkdir ~/scratch");
    // // // 执行命令
    // // ret = system(cmd);
    
    // // snprintf(cmd, sizeof(cmd),
    // //     "python /mnt/f/GitHub/psinad_wbh/scripts/QM.py -d run_bdf -i bdfinp -qm \"bdf\" -t 0 > QMlog");
    // // // 执行命令
    // // ret = system(cmd);


    // FILE *qmlog= fopen("QMlog","r");
    // if (!qmlog) {
    //     printf("Error: Cannot open QMlog\n");
    //     exit(1);
    // }

    // char line[512];
    // char last_line[512] = {0};
    // while (fgets(line, sizeof(line), qmlog)) {
    //     strcpy(last_line, line); // 保存最后一行
    // }
    // fclose(qmlog);

    // if (strstr(last_line, "--- calculation terminated normally ---") == NULL) {
    //     printf("Error: QM finishes unnormally!\n");
    //     exit(-1);
    // };// else {
    // //    printf("yayyy\n");// debug
    // //}

    // FILE *qmdata= fopen("run_bdf/interface.ds","r");
    // if (!qmdata) {
    //     printf("Error: Cannot open interface.ds\n");
    //     exit(1);
    // }
    // for (int i = 0; i < 6; i++){
    //     if (fgets(line, sizeof(line), qmdata) != NULL) {
    //         // printf("%s\n",line);
    //     }
    // }

    // double PES[setm->Nstate];
    // for (int i = 0; i < setm->Nstate; i++) {
    //     if (fscanf(qmdata, "%lf", &PES[i]) != 1) {
            
    //     }
    //     printf("%f\n",PES[i]);
    // }

    // for (int i = 0; i < 4; i++){
    //     if (fgets(line, sizeof(line), qmdata) != NULL) {
    //         printf("%s\n",line);
    //     }
    // }

    // double dPES[setm->Nstate * setm->Natom_mole * 3];
    // for (int j = 0; j < setm->Natom_mole * 3; j++) {
    //     for (int i = 0; i < setm->Nstate; i++) {
    //         if (fscanf(qmdata, "%lf", &dPES[i * setm->Natom_mole * 3 + j]) != 1) {
                  
    //         }
    //         printf("%18.8E  ",dPES[i * setm->Natom_mole * 3 + j]); 
    //     }
        
    //     printf("\n");
    // }


    // for (int i = 0; i < 4; i++){
    //     if (fgets(line, sizeof(line), qmdata) != NULL) {
    //         printf("%s\n",line);
    //     }
    // }



    

    // exit(-1);
}



void qm_mole(double *R, struct set_host *setm, struct set_slave *sets) {
    if (sets->if_recal_qm == 1) {
        FILE *inp = fopen("qminp", "w");
        fprintf(inp, "%s\n", "[GEOM]");
        fprintf(inp, "%s%d\n", "xyz=\"\"\"",setm->Natom_mole);
        fprintf(inp, "%s\n", "Comment: generated by CNAD");
        for (int i = 0; i < setm->Natom_mole; i++){
            fprintf(inp, "%s  %f  %f  %f\n",setm->atomlist_mole[i], R[i * 3 + 0] * au_2_angstrom, R[i * 3 + 1] * au_2_angstrom, R[i * 3 + 2] * au_2_angstrom);
        }

        fprintf(inp, "\"\"\"\n");
        fprintf(inp, "%s\n", setm->qmkeyword_mole);

        fclose(inp);

        char cmd[256];
        int ret;
        snprintf(cmd, sizeof(cmd),"mkdir %s",setm->qmscratch_mole);
        // 执行命令
        ret = system(cmd);

        // printf("python %s -d run_%s -i qminp -qm %s%s%s -t 0 > QMlog\n",
        //     setm->qmscript_mole,setm->qmsoft_mole,"\"",setm->qmsoft_mole,"\"");
        //     exit(-1);


        int ncheck = snprintf(cmd, sizeof(cmd),
            "python %s -d run_qm -i qminp -qm \"%s\" -t 0 > QMlog",
            setm->qmscript_mole,setm->qmsoft_mole);
        // 执行命令
        if (ncheck >= sizeof(cmd)) {
            fprintf(stderr, "Warning: command string truncated!\n");
        }
        ret = system(cmd);
    }

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

    FILE *qmdata= fopen("run_qm/interface.ds","r");
    if (!qmdata) {
        printf("Error: Cannot open interface.ds\n");
        exit(1);
    }
    for (int i = 0; i < 6; i++){
        if (fgets(line, sizeof(line), qmdata) != NULL) {
            // printf("%s\n",line);
        }
    }

    double PES_read[setm->Nstate];
    for (int i = 0; i < setm->Nstate; i++) {
        if (fscanf(qmdata, "%lf", &PES_read[i]) != 1) {
            
        }
        // printf("%f\n",PES[i]);
    }

    for (int i = 0; i < 4; i++){
        if (fgets(line, sizeof(line), qmdata) != NULL) {
            // printf("%s\n",line);
        }
    }

    double dPES_read[setm->Nstate * setm->Natom_mole * 3];
    for (int j = 0; j < setm->Natom_mole * 3; j++) {
        for (int i = 0; i < setm->Nstate; i++) {
            if (fscanf(qmdata, "%lf", &dPES_read[i * setm->Natom_mole * 3 + j]) != 1) {
                  
            }
            // printf("%18.8E  ",dPES[i * setm->Natom_mole * 3 + j]); 
        }
        
        // printf("\n");
    }


    for (int i = 0; i < 4; i++){
        if (fgets(line, sizeof(line), qmdata) != NULL) {
            // printf("%d, %s\n",i,line);
        }
    }

    double nacread[setm->Nstate * setm->Nstate * setm->Natom_mole * 3];
    for (int i = 0; i < setm->Natom_mole * 3; i++) {
        for (int j = 0; j < setm->Nstate * setm->Nstate; j++) {
            if (fscanf(qmdata, "%lf", &nacread[j * setm->Natom_mole * 3 + i]) != 1) {
            }
            // printf("%18.8E  ",nacread[j * setm->Natom_mole * 3 + i]); 
        }
        
        // printf("\n");
    }


    for (int i = 0; i < 4; i++){
        if (fgets(line, sizeof(line), qmdata) != NULL) {
            // printf("%d, %s\n",i,line);
        }
    }


    double complex socread[setm->Nstate * setm->Nstate];
    double re_read, im_read;
    for (int i = 0; i < setm->Nstate; i++) {
        for (int j = 0; j < setm->Nstate; j++) {
            if (fscanf(qmdata, "(%lf, %lf) ", &re_read, &im_read) == 2) {
                socread[i * setm->Nstate + j] = re_read + im_read * I;
            }
            // printf("(%18.8E, %18.8E) ",creal(socread[i * setm->Nstate + j]),cimag(socread[i * setm->Nstate + j])); 
            
        }
        // printf("\n");
        
    }



    if (setm->rep == 1) {
        memset(sets->dv_adia, 0, setm->Nstate * setm->Nstate * setm->Natom_mole * 3 *  sizeof(double complex));
        memset(sets->nac, 0, setm->Nstate * setm->Nstate * setm->Natom_mole * 3 *  sizeof(double complex));
        for (int i = 0; i < setm->Nstate; i++){
            for (int j = 0; j < setm->Nstate; j++){
                if (i == j) {
                    sets->E_adia[i] = PES_read[i];
                    for (int k = 0; k < setm->Natom_mole * 3; k++){
                        sets->dv_adia[i * setm->Nstate * setm->Natom_mole * 3 + i * setm->Natom_mole * 3 + k] = dPES_read[i * setm->Natom_mole * 3 + k];
                    }
                } else {
                    for (int k = 0; k < setm->Natom_mole * 3; k++){
                        sets->nac[i * setm->Nstate * setm->Natom_mole * 3 + j * setm->Natom_mole * 3 + k] = nacread[i * setm->Nstate * setm->Natom_mole * 3 + j * setm->Natom_mole * 3 + k];
                    }
                }
            }
        }
    } else if (setm->rep == 2 || setm->rep == 3) {
        memset(sets->dV, 0, setm->Nstate * setm->Nstate * setm->Natom_mole * 3 *  sizeof(double complex));
        memset(sets->nac, 0, setm->Nstate * setm->Nstate * setm->Natom_mole * 3 *  sizeof(double complex));
        memset(sets->V, 0, setm->Nstate * setm->Nstate * sizeof(double complex));
        for (int i = 0; i < setm->Nstate; i++){
            for (int j = 0; j < setm->Nstate; j++){
                if (i == j) {
                    sets->V[i * setm->Nstate + i] = PES_read[i];
                    for (int k = 0; k < setm->Natom_mole * 3; k++){
                        sets->dV[i * setm->Nstate * setm->Natom_mole * 3 + i * setm->Natom_mole * 3 + k] = dPES_read[i * setm->Natom_mole * 3 + k];
                    }
                } else {
                    for (int k = 0; k < setm->Natom_mole * 3; k++){
                        sets->nac[i * setm->Nstate * setm->Natom_mole * 3 + j * setm->Natom_mole * 3 + k] = nacread[i * setm->Nstate * setm->Natom_mole * 3 + j * setm->Natom_mole * 3 + k];
                        sets->V[i * setm->Nstate + j] = socread[i * setm->Nstate + j];
                    }
                }
            }
        }
    }


    if (setm->rep == 2 || setm->rep == 3) {
        for (int i = 0; i < setm->Nstate * setm->Nstate * setm->Natom_mole * 3; ++i) {
            if (isnan(creal(sets->dV[i]))) {
                printf("Error: read NAN!\n");
                exit(-1);
            }
        }
    }


    // printf("V=\n");
    // for (int i = 0; i < setm->Nstate; i++){
    //     for (int j = 0; j < setm->Nstate; j++){
    //         printf("(%18.8E, %18.8E) ",creal(sets->V[i * setm->Nstate + j]),cimag(sets->V[i * setm->Nstate + j]));
    //     }
    //     printf("\n");
    // }
    
    // printf("dV=\n");
    // for (int i = 0; i < setm->Nstate; i++){
    //     for (int j = 0; j < setm->Nstate; j++){
    //         for (int k = 0; k < setm->Natom_mole; k++){
    //             for (int l = 0; l < 3; l++){
    //                 printf("%18.8E ",creal(sets->dV[i * setm->Nstate * setm->Natom_mole * 3 + j * setm->Natom_mole * 3 + k * 3 + l]));
    //             }
    //             printf("\n");
    //         }
    //     }
    // }




    // printf("nac=\n");
    // for (int i = 0; i < setm->Nstate; i++){
    //     for (int j = 0; j < setm->Nstate; j++){
    //         for (int k = 0; k < setm->Natom_mole; k++){
    //             for (int l = 0; l < 3; l++){
    //                 printf("%18.8E ",creal(sets->nac[i * setm->Nstate * setm->Natom_mole * 3 + j * setm->Natom_mole * 3 + k * 3 + l]));
    //             }
    //             printf("\n");
    //         }
    //     }
    // }


    

    // exit(-1);
}

// Calculate the nuclear force of the model
void nucforce_mole(double *R, double *nf, struct set_host *setm) {
    // for (int j = 0; j < setm->N_bath_mole; j++) {
    //     nf[j] = setm->omega_mole[j]* setm->omega_mole[j] * R[j];
    // }
    
}

