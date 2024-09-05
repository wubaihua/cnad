#include <stdio.h>
#include <stdlib.h>
#include "cJSON.h"

// Function to read the JSON file
char *read_file(const char *filename) {
    FILE *file = fopen(filename, "rb");
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
    return data;
}

// Function to parse the JSON data
void parse_json(const char *json_data) {
    cJSON *json = cJSON_Parse(json_data);
    if (!json) {
        printf("Error parsing JSON: %s\n", cJSON_GetErrorPtr());
        return;
    }

    cJSON *item = json;
    while (item) {
        // Ensure msmodelname is the first parameter read
        cJSON *msmodelname = cJSON_GetObjectItem(item, "msmodelname");
        if (msmodelname && msmodelname->type == cJSON_String) {
            printf("msmodelname: %s\n", msmodelname->valuestring);
        }

        // Read other parameters
        cJSON *N_bath_SBM = cJSON_GetObjectItem(item, "N_bath_SBM");
        if (N_bath_SBM && N_bath_SBM->type == cJSON_Number) {
            printf("N_bath_SBM: %d\n", N_bath_SBM->valueint);
        }

        cJSON *bathtype = cJSON_GetObjectItem(item, "bathtype");
        if (bathtype && bathtype->type == cJSON_Number) {
            printf("bathtype: %d\n", bathtype->valueint);
        }

        cJSON *eps_SBM = cJSON_GetObjectItem(item, "eps_SBM");
        if (eps_SBM && eps_SBM->type == cJSON_Number) {
            printf("eps_SBM: %d\n", eps_SBM->valueint);
        }

        cJSON *delta_SBM = cJSON_GetObjectItem(item, "delta_SBM");
        if (delta_SBM && delta_SBM->type == cJSON_Number) {
            printf("delta_SBM: %d\n", delta_SBM->valueint);
        }

        cJSON *alpha_SBM = cJSON_GetObjectItem(item, "alpha_SBM");
        if (alpha_SBM && alpha_SBM->type == cJSON_Number) {
            printf("alpha_SBM: %f\n", alpha_SBM->valuedouble);
        }

        cJSON *omega_c_SBM = cJSON_GetObjectItem(item, "omega_c_SBM");
        if (omega_c_SBM && omega_c_SBM->type == cJSON_Number) {
            printf("omega_c_SBM: %f\n", omega_c_SBM->valuedouble);
        }

        cJSON *beta = cJSON_GetObjectItem(item, "beta");
        if (beta && beta->type == cJSON_Number) {
            printf("beta: %f\n", beta->valuedouble);
        }

        cJSON *ntraj = cJSON_GetObjectItem(item, "ntraj");
        if (ntraj && ntraj->type == cJSON_Number) {
            printf("ntraj: %d\n", ntraj->valueint);
        }

        cJSON *dt = cJSON_GetObjectItem(item, "dt");
        if (dt && dt->type == cJSON_Number) {
            printf("dt: %f\n", dt->valuedouble);
        }

        cJSON *ttot = cJSON_GetObjectItem(item, "ttot");
        if (ttot && ttot->type == cJSON_Number) {
            printf("ttot: %f\n", ttot->valuedouble);
        }

        cJSON *nbreak = cJSON_GetObjectItem(item, "nbreak");
        if (nbreak && nbreak->type == cJSON_Number) {
            printf("nbreak: %d\n", nbreak->valueint);
        }

        cJSON *method = cJSON_GetObjectItem(item, "method");
        if (method && method->type == cJSON_String) {
            printf("method: %s\n", method->valuestring);
        }

        cJSON *init_occ = cJSON_GetObjectItem(item, "init_occ");
        if (init_occ && init_occ->type == cJSON_Number) {
            printf("init_occ: %d\n", init_occ->valueint);
        }

        cJSON *rep = cJSON_GetObjectItem(item, "rep");
        if (rep && rep->type == cJSON_Number) {
            printf("rep: %d\n", rep->valueint);
        }

        cJSON *outputtype = cJSON_GetObjectItem(item, "outputtype");
        if (outputtype && outputtype->type == cJSON_Number) {
            printf("outputtype: %d\n", outputtype->valueint);
        }

        cJSON *forcetype = cJSON_GetObjectItem(item, "forcetype");
        if (forcetype && forcetype->type == cJSON_Number) {
            printf("forcetype: %d\n", forcetype->valueint);
        }

        cJSON *calforcetype = cJSON_GetObjectItem(item, "calforcetype");
        if (calforcetype && calforcetype->type == cJSON_Number) {
            printf("calforcetype: %d\n", calforcetype->valueint);
        }

        cJSON *nproc_sw = cJSON_GetObjectItem(item, "nproc_sw");
        if (nproc_sw && nproc_sw->type == cJSON_Number) {
            printf("nproc_sw: %d\n", nproc_sw->valueint);
        }

        item = item->next;
    }

    cJSON_Delete(json);
}

int main() {
    const char *filename = "input.json";
    char *json_data = read_file(filename);
    if (json_data) {
        parse_json(json_data);
        free(json_data);
    }
    return 0;
}
