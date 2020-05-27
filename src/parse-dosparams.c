#include "cjson/cJSON.h"
#include "structs.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int read_file(const char *filename, char **content) {
    long length = 0;

    /* open in read binary mode */
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        return 1;
    }

    /* get the length */
    if (fseek(file, 0, SEEK_END) != 0) {
        fclose(file);
        return 1;
    }
    length = ftell(file);
    if (length < 0) {
        fclose(file);
        return 1;
    }
    if (fseek(file, 0, SEEK_SET) != 0) {
        fclose(file);
        return 1;
    }

    /* allocate content buffer */
    *content = (char *)malloc((size_t)length + 1);
    if (*content == NULL) {
        fclose(file);
        return 1;
    }

    /* read the file into memory */
    long read_chars = fread(*content, sizeof(char), (size_t)length, file);
    if (read_chars != length) {
        free(*content);
        fclose(file);
        return 1;
    }
    (*content)[read_chars] = '\0';
    fclose(file);

    return 0;
}

void json_parse_file(const char *file, cJSON **json) {
    char *file_content = NULL;
    if (read_file(file, &file_content) != 0) {
        fprintf(stderr, "ERROR: file '%s' could not be read in.\n", file);
        exit(1);
    }

    // parse
    *json = cJSON_Parse(file_content);
    if (*json == NULL) {
        const char *error_ptr = cJSON_GetErrorPtr();
        fprintf(stderr, "ERROR: Could not parse json. Error before: %s.\n",
                error_ptr);
        cJSON_Delete(*json);
        exit(1);
    }
    free(file_content);
}

int json_parse_int(cJSON *json, char *key) {
    const cJSON *json_number = NULL;
    json_number = cJSON_GetObjectItemCaseSensitive(json, key);
    if (!cJSON_IsNumber(json_number)) {
        fprintf(stderr, "ERROR: %s could not be parsed.", key);
        cJSON_Delete(json);
        exit(1);
    }
    return json_number->valueint;
}

int json_parse_float(cJSON *json, char *key) {
    const cJSON *json_number = NULL;
    json_number = cJSON_GetObjectItemCaseSensitive(json, key);
    if (!cJSON_IsNumber(json_number)) {
        fprintf(stderr, "ERROR: %s could not be parsed.", key);
        cJSON_Delete(json);
        exit(1);
    }
    return (float)json_number->valuedouble;
}

char json_parse_char(cJSON *json, char *key) {
    const cJSON *json_string = NULL;
    json_string = cJSON_GetObjectItemCaseSensitive(json, key);
    if (!cJSON_IsString(json_string)) {
        fprintf(stderr, "ERROR: %s could not be parsed.", key);
        cJSON_Delete(json);
        exit(1);
    }
    if (strlen(json_string->valuestring) != 1) {
        fprintf(stderr, "ERROR: %s is not a single char.", key);
        cJSON_Delete(json);
        exit(1);
    }
    return json_string->valuestring[0];
}

char *json_parse_string(cJSON *json, char *key) {
    const cJSON *json_string = NULL;
    json_string = cJSON_GetObjectItemCaseSensitive(json, key);
    if (!cJSON_IsString(json_string)) {
        fprintf(stderr, "ERROR: %s could not be parsed.", key);
        cJSON_Delete(json);
        exit(1);
    }
    if (strlen(json_string->valuestring) > 80) {
        fprintf(stderr, "ERROR: %s is longer than 80 chars.", key);
        cJSON_Delete(json);
        exit(1);
    }
    return json_string->valuestring;
}

float json_parse_float_nokey(cJSON *json_number) {
    if (!cJSON_IsNumber(json_number)) {
        fprintf(stderr, "ERROR: float could not be parsed.");
        cJSON_Delete(json_number);
        exit(1);
    }
    return json_number->valuedouble;
}

int json_parse_int_nokey(cJSON *json_number) {
    if (!cJSON_IsNumber(json_number)) {
        fprintf(stderr, "ERROR: int could not be parsed.");
        cJSON_Delete(json_number);
        exit(1);
    }
    return json_number->valueint;
}

void json_parse_json_array(cJSON *json, char *key, cJSON **array) {
    *array = cJSON_GetObjectItemCaseSensitive(json, key);
    if (!cJSON_IsArray(*array)) {
        fprintf(stderr, "ERROR: '%s' array could not be parsed.", key);
        cJSON_Delete(json);
        exit(1);
    }
}

int parse_dosparams(const char *dosparams_file,
                    size_t *nsamples, // output
                    size_t *nblocks, unsigned long *nblocksteps,
                    size_t *nmoltypes, size_t **moltypes_nmols,
                    size_t **moltypes_natomspermol,
                    float ***moltypes_atommasses, char **moltypes_rot_treat,
                    int ***moltypes_abc_indicators, size_t *ncross_spectra,
                    cross_spectrum_def **cross_spectra_def) {
    // tokenize file
    cJSON *dosparams_json;
    json_parse_file(dosparams_file, &dosparams_json);

    // parse numbers
    *nsamples = (size_t)json_parse_int(dosparams_json, "nsamples");
    *nblocks = (size_t)json_parse_int(dosparams_json, "nblocks");
    *nblocksteps = (size_t)json_parse_int(dosparams_json, "nblocksteps");

    // parse moltypes array
    cJSON *moltypes_json = NULL;
    json_parse_json_array(dosparams_json, "moltypes", &moltypes_json);
    *nmoltypes = cJSON_GetArraySize(moltypes_json);

    // allocate moltypes arrays
    *moltypes_nmols = calloc(*nmoltypes, sizeof(size_t));
    *moltypes_natomspermol = calloc(*nmoltypes, sizeof(size_t));
    *moltypes_atommasses = (float **)malloc(*nmoltypes * sizeof(float *));
    *moltypes_rot_treat = calloc(*nmoltypes, sizeof(char));
    *moltypes_abc_indicators = (int **)malloc(*nmoltypes * sizeof(int *));

    // parse each moltype
    for (size_t h = 0; h < *nmoltypes; h++) {
        // moltype
        cJSON *moltype_json = NULL;
        moltype_json = cJSON_GetArrayItem(moltypes_json, h);
        // nmols
        (*moltypes_nmols)[h] = (size_t)json_parse_int(moltype_json, "nmols");
        // atommasses
        cJSON *atommasses_json = NULL;
        json_parse_json_array(moltype_json, "atom_masses", &atommasses_json);
        (*moltypes_natomspermol)[h] = cJSON_GetArraySize(atommasses_json);
        (*moltypes_atommasses)[h] =
            (float *)malloc((*moltypes_natomspermol)[h] * sizeof(float));
        for (size_t j = 0; j < (*moltypes_natomspermol)[h]; j++) {
            cJSON *atommass = cJSON_GetArrayItem(atommasses_json, j);
            (*moltypes_atommasses)[h][j] = json_parse_float_nokey(atommass);
        }
        // rot_treat
        (*moltypes_rot_treat)[h] = json_parse_char(moltype_json, "rot_treat");
        // abc_indicators
        cJSON *abc_indicators_json = NULL;
        json_parse_json_array(moltype_json, "abc_indicators",
                              &abc_indicators_json);
        (*moltypes_abc_indicators)[h] = (int *)malloc(4 * sizeof(int));
        for (size_t j = 0; j < 4; j++) {
            cJSON *abc_indicator = cJSON_GetArrayItem(abc_indicators_json, j);
            (*moltypes_abc_indicators)[h][j] =
                (size_t)json_parse_int_nokey(abc_indicator);
            if (((*moltypes_abc_indicators)[h][j] == -1) &&
                (j == 0 || j == 2)) {
                fprintf(stderr, "ERROR: only the second and fourth abc "
                                "indicator can be -1.\n");
                exit(1);
            }
        }
    }

    // parse cross_spectra array
    cJSON *cross_spectra = NULL;
    json_parse_json_array(dosparams_json, "cross_spectra", &cross_spectra);
    *ncross_spectra = cJSON_GetArraySize(cross_spectra);

    // allocate cross_spectra arrays
    *cross_spectra_def = (cross_spectrum_def *)malloc(
        *ncross_spectra * sizeof(**cross_spectra_def));

    // parse each cross_spectrum
    for (size_t d = 0; d < *ncross_spectra; d++) {
        // cross_spectrum
        cJSON *cross_spectrum = NULL;
        cross_spectrum = cJSON_GetArrayItem(cross_spectra, d);
        // name
        strcpy((*cross_spectra_def)[d].name,
               json_parse_string(cross_spectrum, "name"));
        // type
        (*cross_spectra_def)[d].type = json_parse_char(cross_spectrum, "type");
        // dof_pairs
        cJSON *dof_pairs = NULL;
        json_parse_json_array(cross_spectrum, "dof_pairs", &dof_pairs);
        (*cross_spectra_def)[d].ndof_pair_defs = cJSON_GetArraySize(dof_pairs);
        (*cross_spectra_def)[d].dof_pair_defs = (dof_pair_def *)malloc(
            (*cross_spectra_def)[d].ndof_pair_defs * sizeof(dof_pair_def));

        for (size_t p = 0; p < (*cross_spectra_def)[d].ndof_pair_defs; p++) {
            cJSON *dof_pair = cJSON_GetArrayItem(dof_pairs, p);
            cJSON *dofA = cJSON_GetArrayItem(dof_pair, 0);
            cJSON *dofB = cJSON_GetArrayItem(dof_pair, 1);
            // dof{A,B}_moletype
            (*cross_spectra_def)[d].dof_pair_defs[p].dofA_moltype =
                (size_t)json_parse_int(dofA, "moltype");
            (*cross_spectra_def)[d].dof_pair_defs[p].dofB_moltype =
                (size_t)json_parse_int(dofB, "moltype");
            if ((*cross_spectra_def)[d].dof_pair_defs[p].dofA_moltype >=
                    *nmoltypes ||
                (*cross_spectra_def)[d].dof_pair_defs[p].dofB_moltype >=
                    *nmoltypes) {
                fprintf(stderr, "ERROR: moltype in cross spectrum higher than "
                                "nmoltypes.\n");
                exit(1);
            }
            if ((*cross_spectra_def)[d].type == 'i' &&
                ((*cross_spectra_def)[d].dof_pair_defs[p].dofA_moltype !=
                 (*cross_spectra_def)[d].dof_pair_defs[p].dofB_moltype)) {
                fprintf(stderr,
                        "ERROR: For cross_spectrum type 'i' both moltypes have "
                        "to match.\n");
                exit(1);
            }
            // dof{A,B}_type
            (*cross_spectra_def)[d].dof_pair_defs[p].dofA_type =
                json_parse_char(dofA, "dof_type");
            (*cross_spectra_def)[d].dof_pair_defs[p].dofB_type =
                json_parse_char(dofB, "dof_type");
            // dof{A,B}_list und ndof{A,B}
            cJSON *dofA_list = NULL;
            cJSON *dofB_list = NULL;
            json_parse_json_array(dofA, "dof_list", &dofA_list);
            json_parse_json_array(dofB, "dof_list", &dofB_list);
            size_t ndofA = cJSON_GetArraySize(dofA_list);
            size_t ndofB = cJSON_GetArraySize(dofB_list);
            (*cross_spectra_def)[d].dof_pair_defs[p].ndofA = ndofA;
            (*cross_spectra_def)[d].dof_pair_defs[p].ndofB = ndofB;
            (*cross_spectra_def)[d].dof_pair_defs[p].dofA_list =
                (size_t *)malloc(ndofA * sizeof(size_t));
            (*cross_spectra_def)[d].dof_pair_defs[p].dofB_list =
                (size_t *)malloc(ndofB * sizeof(size_t));
            for (size_t l = 0; l < ndofA; l++) {
                cJSON *dof_indicator = cJSON_GetArrayItem(dofA_list, l);
                (*cross_spectra_def)[d].dof_pair_defs[p].dofA_list[l] =
                    (size_t)json_parse_int_nokey(dof_indicator);
            }
            for (size_t l = 0; l < ndofB; l++) {
                cJSON *dof_indicator = cJSON_GetArrayItem(dofB_list, l);
                (*cross_spectra_def)[d].dof_pair_defs[p].dofB_list[l] =
                    (size_t)json_parse_int_nokey(dof_indicator);
            }
        }
    }
    cJSON_Delete(dosparams_json);
    return 0;
}

void print_dosparams(size_t nsamples, size_t nblocks, unsigned long nblocksteps,
                     size_t nmoltypes, size_t *moltypes_nmols,
                     size_t *moltypes_natomspermol, float **moltypes_atommasses,
                     char *moltypes_rot_treat, int **moltypes_abc_indicators) {
    printf("nsamples: %zu\n", nsamples);
    printf("nblocks: %zu\n", nblocks);
    printf("nblocksteps: %lu\n", nblocksteps);
    printf("nmoltypes: %zu\n", nmoltypes);
    for (size_t h = 0; h < nmoltypes; h++) {
        printf("moltype %zu nmols: %zu\n", h, moltypes_nmols[h]);
        printf("moltype %zu natomspermol: %zu\n", h, moltypes_natomspermol[h]);
        printf("moltype %zu atommasses: ", h);
        for (size_t j = 0; j < moltypes_natomspermol[h]; j++)
            printf("%f ", moltypes_atommasses[h][j]);
        printf("\n");
        printf("moltype %zu rot_treat: %c\n", h, moltypes_rot_treat[h]);
        printf("moltype %zu abc_indicators: ", h);
        for (size_t j = 0; j < 4; j++)
            printf("%d ", moltypes_abc_indicators[h][j]);
        printf("\n");
    }
}

void free_dosparams_arrays(size_t nmoltypes, size_t **moltypes_nmols,
                           size_t **moltypes_natomspermol,
                           float ***moltypes_atommasses,
                           char **moltypes_rot_treat,
                           int ***moltypes_abc_indicators,
                           size_t *ncross_spectra,
                           cross_spectrum_def **cross_spectra_def) {
    // free arrays
    free(*moltypes_nmols);
    free(*moltypes_natomspermol);

    for (size_t h = 0; h < nmoltypes; h++) {
        free((*moltypes_atommasses)[h]);
        free((*moltypes_abc_indicators)[h]);
    }
    free(*moltypes_atommasses);
    free(*moltypes_abc_indicators);
    free(*moltypes_rot_treat);

    for (size_t d = 0; d < *ncross_spectra; d++) {
        for (size_t p = 0; p < (*cross_spectra_def)[d].ndof_pair_defs; p++) {
            free((*cross_spectra_def)[d].dof_pair_defs[p].dofA_list);
            free((*cross_spectra_def)[d].dof_pair_defs[p].dofB_list);
        }
        free((*cross_spectra_def)[d].dof_pair_defs);
    }
    free(*cross_spectra_def);
}

void calc_convenience_variables(size_t nmoltypes, size_t *moltypes_nmols,
                                size_t *moltypes_natomspermol,
                                float **moltypes_atommasses,
                                size_t *natoms, // output
                                size_t *nmols, size_t *moltypes_firstmol,
                                size_t *moltypes_firstatom,
                                size_t **mols_moltypenr, size_t **mols_natoms,
                                float **mols_mass, size_t **mols_firstatom) {
    for (size_t h = 0; h < nmoltypes; h++) {
        moltypes_firstmol[h] = *nmols;
        moltypes_firstatom[h] = *natoms;
        *nmols += moltypes_nmols[h];
        *natoms += moltypes_nmols[h] * moltypes_natomspermol[h];
    }
    *mols_moltypenr = calloc(*nmols, sizeof(size_t));
    *mols_natoms = calloc(*nmols, sizeof(size_t));
    *mols_mass = calloc(*nmols, sizeof(float));
    *mols_firstatom = calloc(*nmols, sizeof(size_t));
    for (size_t i = 0; i < *nmols; i++) {
        for (size_t h = 0; h < nmoltypes; h++) {
            if (i >= moltypes_firstmol[h] &&
                i < moltypes_firstmol[h] + moltypes_nmols[h])
                (*mols_moltypenr)[i] = h;
        }
        (*mols_natoms)[i] = moltypes_natomspermol[(*mols_moltypenr)[i]];
        (*mols_mass)[i] = 0;
        for (size_t j = 0; j < (*mols_natoms)[i]; j++)
            (*mols_mass)[i] += moltypes_atommasses[(*mols_moltypenr)[i]][j];
    }
    for (size_t h = 0; h < nmoltypes; h++) {
        for (size_t ii = 0; ii < moltypes_nmols[h]; ii++) {
            (*mols_firstatom)[moltypes_firstmol[h] + ii] =
                moltypes_firstatom[h] + moltypes_natomspermol[h] * ii;
        }
    }
}
