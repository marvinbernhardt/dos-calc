#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "cjson/cJSON.h"

int write_to_file(const char *filename, char *content)
{
    FILE *file = fopen(filename, "wb");
    if (file == NULL) return 1;
    fputs(content, file);
    fclose(file);
    return 0;
}

int write_dos(const char *dos_file,
        size_t nsamples,
        unsigned long nblocksteps,
        unsigned long nfrequencies,
        float framelength,
        size_t ndos,
        size_t ncross_spectra,
        const char **dos_names,
        size_t nmoltypes,
        float *moltypes_dos_samples,
        float *moltypes_dos_cross_samples,
        float *moltypes_samples_moments_of_inertia,
        cross_spectrum_def *cross_spectra_def)
{

    cJSON *root = cJSON_CreateObject();
    if (root == NULL) return 1;

    // frequencies
    cJSON *frequencies = cJSON_AddArrayToObject(root, "frequencies");
    if (frequencies == NULL) return 1;
    for (unsigned long t=0; t<nfrequencies; t++)
    {
        cJSON *number = cJSON_CreateNumber(t / (framelength * (float)nblocksteps));
        if (number == NULL) return 1;
        cJSON_AddItemToArray(frequencies, number);
    }

    // moltypes array
    cJSON *moltypes = cJSON_AddArrayToObject(root, "moltypes");
    if (moltypes == NULL) return 1;
    for (size_t h=0; h<nmoltypes; h++)
    {
        // moltype object
        cJSON *moltype = cJSON_CreateObject();
        if (moltype == NULL) return 1;
        cJSON_AddItemToArray(moltypes, moltype);

        // dos array
        cJSON *dos = cJSON_AddArrayToObject(moltype, "dos");
        if (dos == NULL) return 1;

        // moi array
        cJSON *moi = cJSON_AddArrayToObject(moltype, "moments_of_inertia");
        if (dos == NULL) return 1;

        // fill dos array
        for (size_t d=0; d<ndos; d++)
        {
            // this_dos object
            cJSON *this_dos = cJSON_CreateObject();
            if (this_dos == NULL) return 1;
            cJSON_AddItemToArray(dos, this_dos);

            // this_dos.name
            if (cJSON_AddStringToObject(this_dos, "name", dos_names[d]) == NULL) return 1;

            // this_dos.data
            cJSON *dos_data = cJSON_AddArrayToObject(this_dos, "data");
            if (dos_data == NULL) return 1;

            // fill dos_data array
            for (size_t sample=0; sample<nsamples; sample++)
            {
                // dos_data_sample array
                cJSON *dos_data_sample = cJSON_CreateArray();
                if (dos_data_sample == NULL) return 1;
                cJSON_AddItemToArray(dos_data, dos_data_sample);

                // fill dos_data_sample array
                for (unsigned long t=0; t<nfrequencies; t++)
                {
                    int index =
                        h*ndos*nsamples*nfrequencies
                        +    d*nsamples*nfrequencies
                        +        sample*nfrequencies
                        +                       t;
                    // single number
                    cJSON *number = cJSON_CreateNumber(moltypes_dos_samples[index]);
                    if (number == NULL) return 1;
                    cJSON_AddItemToArray(dos_data_sample, number);
                }
            }
        }

        // fill moi array
        for (size_t sample=0; sample<nsamples; sample++)
        {
            // this_moi array
            cJSON *this_moi = cJSON_CreateArray();
            if (this_moi == NULL) return 1;
            cJSON_AddItemToArray(moi, this_moi);
            
            // fill this_moi array
            for (size_t abc=0; abc<3; abc++)
            {
                size_t index =
                      h*nsamples*  3
                    +     sample*  3
                    +            abc;
                // single number
                cJSON *number = cJSON_CreateNumber(moltypes_samples_moments_of_inertia[index]);
                if (number == NULL) return 1;
                cJSON_AddItemToArray(this_moi, number);
            }
        }
    }

    // dos_cross array
    cJSON *dos_cross = cJSON_AddArrayToObject(root, "cross_spectra");
    if (dos_cross == NULL) return 1;

    // fill dos_cross array
    for (size_t d=0; d<ncross_spectra; d++)
    {
        // this_dos object
        cJSON *this_dos = cJSON_CreateObject();
        if (this_dos == NULL) return 1;
        cJSON_AddItemToArray(dos_cross, this_dos);

        // this_dos.name
        if (cJSON_AddStringToObject(this_dos, "name", cross_spectra_def[d].name) == NULL) return 1;

        // this_dos.data
        cJSON *dos_data = cJSON_AddArrayToObject(this_dos, "data");
        if (dos_data == NULL) return 1;

        // fill dos_data array
        for (size_t sample=0; sample<nsamples; sample++)
        {
            // dos_data_sample array
            cJSON *dos_data_sample = cJSON_CreateArray();
            if (dos_data_sample == NULL) return 1;
            cJSON_AddItemToArray(dos_data, dos_data_sample);

            // fill dos_data_sample array
            for (unsigned long t=0; t<nfrequencies; t++)
            {
                size_t index =
                    +          d*nsamples*nfrequencies
                    +              sample*nfrequencies
                    +                             t;
                // single number
                cJSON *number = cJSON_CreateNumber(moltypes_dos_cross_samples[index]);
                if (number == NULL) return 1;
                cJSON_AddItemToArray(dos_data_sample, number);
            }
        }
    }

    // generate json string
    char *string = cJSON_Print(root);
    if (string == NULL) {
        fprintf(stderr, "ERROR: Failed generate json.\n");
        return 1;
    }
    cJSON_Delete(root);

    // write string to file
    if (write_to_file(dos_file, string) != 0) {
        fprintf(stderr, "ERROR: Failed write json to file.\n");
        return 1;
    }
    free(string);

    return 0;
}
