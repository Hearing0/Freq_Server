#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>      // FFT transform library
// #include <cblas.h>      // Matrix Multi library
// #include <utils/phasing_utils.c>
#include <string.h>
#include "../utils/ini_parser/ini.c"
#include "../utils/read_config.c"
#include <time.h>


#ifndef TEST_SAMPLES
#define TEST_SAMPLES 0
#endif

#ifndef TEST_CLR_RANGE
#define TEST_CLR_RANGE 0
#endif

#ifndef CLR_BANDS_MAX
#define CLR_BANDS_MAX 6
#endif


typedef struct sample_meta_data {
    int *antenna_list;
    int num_antennas;
    int number_of_samples;
    double x_spacing;
    int usrp_rf_rate;
    int usrp_fcenter;
} sample_meta_data;

typedef struct freq_data {
    double *restricted_freq;
    double *clear_freq_range;
} freq_data;

typedef struct freq_band {
    int f_start;
    int f_end;
    double noise;
} freq_band;


typedef struct clear_freq {
    double noise;
    double tfreq;
} clear_freq;


/**
 * @brief  Writes a complex Frequency Spectrum to csv file to be plotted in python.
 * @note   By DF
 * @param  *filename:       The filepath for the saved CSV file
 * @param  *spectrum:       Frequency Spectrum in complex form
 * @param  *freq_vector:    Array for Frequency (steps of frequency per sample)
 * @param  num_samples:     Number of four_spectrums collected per antenna
 * @retval None
 * @deprecated Replaced by write_spectrum_mag_csv. No longer using complex form to similify calculations
 */
void write_spectrum_csv(char *filename, fftw_complex *spectrum, double *freq_vector, int num_samples) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Frequency,Power\n");
    for (int i = 0; i < num_samples; i++) {
        double magnitude = sqrt(creal(spectrum[i]) * creal(spectrum[i]) + cimag(spectrum[i]) * cimag(spectrum[i]));
        fprintf(file, "%f,%f\n", freq_vector[i], magnitude);
    }

    fclose(file);
}

/**
 * @brief  Writes a magnitude Frequency Spectrum to csv file to be plotted in python.
 * @note   By DF
 * @param  *filename:       The filepath for the saved CSV file
 * @param  *spectrum:       Frequency Spectrum in magnitude form
 * @param  *freq_vector:    Array for Frequency (steps of frequency per sample)
 * @param  num_samples:     Number of samples in spectrum
 * @retval None
 */
void write_spectrum_mag_csv(char *filename, double *spectrum, double *freq_vector, int num_samples) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Frequency,Power\n");
    for (int i = 0; i < num_samples; i++) {
        fprintf(file, "%f,%f\n", freq_vector[i], spectrum[i]);
    }

    fclose(file);
}


void write_clr_freq_csv(char *filename, freq_band *clr_bands) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    // Find Start and End of Clear Freq Range to Write later
    int clr_start = RAND_MAX;
    int clr_end = 0;
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        if (clr_bands[i].f_start < clr_start && clr_bands[i].noise < RAND_MAX) clr_start = clr_bands[i].f_start;
        if (clr_bands[i].f_end > clr_end && clr_bands[i].noise < RAND_MAX) clr_end = clr_bands[i].f_end;
    }    

    fprintf(file, "Start Frequency,End Frequency,Noise,Clear Freq Start,Clear Freq End\n");
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        // Special: Print Clear Freq Range on Line 0
        if (i == 0) fprintf(file, "%d,%d,%f,%d,%d\n", clr_bands[i].f_start, clr_bands[i].f_end, clr_bands[i].noise,clr_start,clr_end);
        else fprintf(file, "%d,%d,%f\n", clr_bands[i].f_start, clr_bands[i].f_end, clr_bands[i].noise);
    }

    fclose(file);
}

/**
 * @brief  Writes the Real/Imaginary magnitude to csv file to be plotted in python.
 * @note   By DF
 * @param  *filename:           The filepath for the saved CSV file
 * @param  *raw_samples_mag:    Int array of the Real/Imaginary magnitude
 * @param  *freq_vector:        Array for Frequency (steps of frequency per sample)
 * @param  num_samples:         Number of four_spectrums collected per antenna
 * @retval None
 */
void write_sample_mag_csv(char *filename, int **raw_samples_mag, double *freq_vector, sample_meta_data *meta_data) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Samples,Power\n");
    for (int j = 0; j < meta_data->num_antennas; j++) {
        for (int i = 0; i < meta_data->number_of_samples; i++) {
            fprintf(file, "%f,%d\n", freq_vector[i], raw_samples_mag[j][i]);
        }
    }
    fclose(file);
}


/**
 * @brief  Reads intial Clear Freq parameters from a converted txt file (see utils/pickle_text_convert.py)
 * @note   
 * @param  *filename:           Filepath of converted txt file
 * @param  *meta_data:          Meta data for current sample batch
 * @param  **clear_freq_range:  Range for the Clear Freq
 * @param  ***raw_samples:      14x2500 complex sample array
 * @retval None
 */
void read_input_data(const char *filename, sample_meta_data *meta_data, double **clear_freq_range, fftw_complex ***raw_samples) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line[256];
    int antenna_list_size = 0;

    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "number_of_samples: %d", &meta_data->number_of_samples) == 1) continue;
        if (sscanf(line, "usrp_rf_rate: %d", &meta_data->usrp_rf_rate) == 1) continue;
        if (sscanf(line, "usrp_fcenter: %d", &meta_data->usrp_fcenter) == 1) continue;
        // if (sscanf(line, "x_spacing: %lf", &meta_data->x_spacing) == 1) continue;
        if (strncmp(line, "clear_freq_range:", 15) == 0 && TEST_CLR_RANGE) {
            clear_freq_range = realloc(clear_freq_range, 2 * sizeof(double));
            int i = 0;
            char *token = strtok(line + 16, ",");
            while (token != NULL) {
                (*clear_freq_range)[i] = atof(token);
                i++;
            }
            continue;
        }
        
        // Antenna List
        if (strncmp(line, "antenna_list:", 13) == 0) {
            char *token = strtok(line + 14, ",");
            while (token != NULL) {
                // Remove any leading or trailing whitespace from token
                while (isspace(*token)) token++;
                char *end = token + strlen(token) - 1;
                while (end > token && isspace(*end)) end--;
                *(end + 1) = '\0';

                meta_data->antenna_list = realloc(meta_data->antenna_list, (++antenna_list_size) * sizeof(int));
                meta_data->antenna_list[antenna_list_size - 1] = atoi(token);

                token = strtok(NULL, ",");
            }
            meta_data->num_antennas = 2; //antenna_list_size;
            continue;
        }

        // Raw Samples
        if (strncmp(line, "raw_samples:", 12) == 0 && TEST_SAMPLES) {
            printf("[Clear Freq Search] Aquiring test four_spectrums from pickle files...\n");
            // Allocate mem
            *raw_samples = (fftw_complex **)fftw_malloc(meta_data->num_antennas * sizeof(fftw_complex *));
            for (int i = 0; i < meta_data->num_antennas; i++) {
                (*raw_samples)[i] = (fftw_complex *)fftw_malloc(meta_data->number_of_samples * sizeof(fftw_complex));
            }
            if (*raw_samples == NULL) {
                perror("Error allocating memory for raw four_spectrums");
                exit(EXIT_FAILURE);
            }
            
            // Store data
            for (int i = 0; i < meta_data->num_antennas; i++) {
                fftw_complex *ant_samples = (*raw_samples)[i];

                for (int j = 0; j < meta_data->number_of_samples; j++) {
                    double real, imag;
                    fgets(line, sizeof(line), file);
                    sscanf(line, "%lf,%lf", &real, &imag);
                    
                    ant_samples[j] = real + I * imag;
                }
            }
            break;
        }
    }

    fclose(file);
}
