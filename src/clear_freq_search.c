#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>      // FFT transform library
#include <cblas.h>      // Matrix Multi library
// #include <utils/phasing_utils.c>
#include <string.h>
#include <math.h>
#include "../utils/ini_parser/ini.c"

#define PI 3.14159265358979323846

//
int C = (int)pow(3.0, 8.0);
bool verbose = false;

// #define RESTRICT_FILE = '/home/radar/repos/SuperDARN_MSI_ROS/linux/home/radar/ros.3.6/tables/superdarn/site/site.sps/restrict.dat.inst'


/*
* NOTE: 
* When compiling make sure to add "-lm" and "-lfftw3" to link the libraries 
* to the compilier. 
*/


typedef struct {
    int *antenna_list;
    int num_antennas;
    int number_of_samples;
    double x_spacing;
    int usrp_rf_rate;
    int usrp_fcenter;
} SampleMetaData;

void read_input_data(const char *filename, SampleMetaData *meta_data, fftw_complex **raw_samples) {
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
        if (sscanf(line, "x_spacing: %lf", &meta_data->x_spacing) == 1) continue;

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
            meta_data->num_antennas = antenna_list_size;
            continue;
        }

        if (strncmp(line, "raw_samples:", 12) == 0) {
            *raw_samples = malloc((meta_data->number_of_samples) * sizeof(fftw_complex));
            if (*raw_samples == NULL) {
                perror("Error allocating memory for raw samples");
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < meta_data->number_of_samples; i++) {
                double real, imag;
                fgets(line, sizeof(line), file);
                sscanf(line, "%lf,%lf", &real, &imag);
                (*raw_samples)[i] = real + I * imag;
            }
            break;
        }
    }

    fclose(file);
}

/**
 * @brief  Calculates phase increment between antennas to produce a mainlobe sterring 
 *         of beam_angle at center_frequency.
 * 
 * @param  beam_angle: Distrubution of the beam (in radians)
 * @param  center_frequency: Frequency at the center of the beam; 
 * * used to phase-shift allign the other frequencys (in Hz)
 * @param  x_spacing: Spacing in-between antennas
 *         
 * @retval phase_shift: Phase shift (in degrees)
 */
double calc_phase_increment(double beam_angle, double center_frequency, double x_spacing) {
    double wavelength = C / center_frequency;
    double phase_shift = (2 * PI * x_spacing * sin(beam_angle)) / wavelength;
    if (verbose) {
        printf("center_freq: %lf x_spacing: %lf \nphase_shift: %lf degree", center_frequency, x_spacing, phase_shift * 180 / PI);
    }
    return phase_shift; 
}

/**
 * Converts radians to rectangular form
 */
double complex rad_to_rect(double phase) {
    return cexp(I * phase);
}

/**
 * @brief  Fast Fourier Transform (FFT) conversion from beamformed_samples 
 * * into a spectrum.
 * @note   // HACK Check time vs regular fft; plans are expensive to create. 
 * @param  *beamformed_samples: Array of samples that have already been 
 * * beamformed
 * @param  number_of_samples: Number of samples
 * @param  *spectrum_power: Output array for the resultant spectrum
 * @retval None
 */
void fft_clrfreq_samples(fftw_complex *beamformed_samples, int num_samples, fftw_complex *spectrum_power) {
    fftw_plan plan = fftw_plan_dft_1d(num_samples, beamformed_samples, spectrum_power, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

void find_clear_freq(fftw_complex *spectrum_power, double *freq_vector, double start_freq, double end_freq, double clear_bw, double *tfreq, double *noise) {

}

// HACK apply efficient matrix multi via cblas_dgemm
void calc_clear_freq_on_raw_samples(fftw_complex *raw_samples, SampleMetaData *meta_data, double *restricted_frequencies, double *clear_freq_range, double beam_angle, double smsep) {
    // Extract meta data
    int num_samples = meta_data->number_of_samples;
    int *antennas = meta_data->antenna_list;

    // Ensure inputs exist
    if (!raw_samples || !meta_data || !antennas) {
        fprintf(stderr, "Error: Null input detected.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for Variables
    fftw_complex *phasing_vector = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_complex *beamformed_samples = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_complex *spectrum_power = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    double *freq_vector = (double*) malloc(sizeof(double) * num_samples);
    
    if (!phasing_vector || !beamformed_samples || !spectrum_power || !freq_vector) {
        perror("Error allocating memory");
        exit(EXIT_FAILURE);
    }

    // Calculate and Apply phasing vector
    double phase_increment = calc_phase_increment(beam_angle, (clear_freq_range[0] + clear_freq_range[1]) / 2 * 1000, meta_data->x_spacing);
    printf("phase_increment: %lf\n", phase_increment);

    printf("antenna size: %d\n", meta_data->num_antennas);

    for (int i = 0; i < meta_data->num_antennas; i++) {
        if (i <= 15 || i >= 20) {
            if (i < meta_data->num_antennas) {
                printf("antenna[%d]: %d\n", i, antennas[i]);
                phasing_vector[i] = rad_to_rect(antennas[i] * phase_increment);
                printf("phase_vector[%d]: %f + %fi\n", i, creal(phasing_vector[i]), cimag(phasing_vector[i]));
            } else {
                fprintf(stderr, "Error: Accessing antennas out of bounds at index %d\n", i);
                exit(EXIT_FAILURE);
            }
        }
    }

    // // Apply beamforming
    // for (int i = 0; i < num_samples; i++) {
    //     beamformed_samples[i] = phasing_vector[i] * raw_samples[i];
    //     printf("beamformed_samples[%d]: %f + %fi\n", i, creal(beamformed_samples[i]), cimag(beamformed_samples[i]));
    // }

    // // Spectral Estimation
    // fft_clrfreq_samples(beamformed_samples, num_samples, spectrum_power);

    // // Print the spectrum power
    // for (int i = 0; i < num_samples; i++) {
    //     printf("Frequency bin %d: %f + %fi\n", i, creal(spectrum_power[i]), cimag(spectrum_power[i]));
    // }


    /// END of Spectrum Calc


    // // Frequency Vector Calculation
    // double delta_f = meta_data->usrp_rf_rate / num_samples;
    // for (int i = 0; i < num_samples; i++) {
    //     freq_vector[i] = i * delta_f - (meta_data->usrp_rf_rate / 2) + meta_data->usrp_fcenter * 1000;
    // }

    // // Mask restricted frequencies
    // // 

    // // Find clear frequency
    // double clear_bw = 2e6 / smsep;
    // double tfreq, noise;
    // find_clear_freq(spectrum_power, freq_vector, clear_freq_range[0] * 1e3, clear_freq_range[1] * 1e3, clear_bw, &tfreq, &noise);

    // // Output results
    // printf("Clear Frequency: %f, Noise: %f\n", tfreq, noise);

    printf("Finished Clear Freq!");
    
    fftw_free(phasing_vector);
    fftw_free(beamformed_samples);
    fftw_free(spectrum_power);
    free(freq_vector);

}

int main() {
    // Stopwatch Start
    // double t1 = dsecnd();

    // HACK: Setup file_path environment variable
    const char *input_file_path = "../Freq_Server/utils/clear_freq_input/clrfreq_dump.1.txt";
    const char *config_file_path = "../Freq_Server/utils/clear_freq_input/array_config.ini";
    // const char *output_file_path = "../utils/txt_output/result.txt";

    // Initial Data Variables
    SampleMetaData meta_data = {0};
    fftw_complex *raw_samples = NULL;

    // Load raw_samples and meta_data
    read_input_data(input_file_path, &meta_data, &raw_samples);
       
    printf("num_samples: %d\nx_spacing: %lf\nusrp_rf_rate: %d\nusrp_fcenter: %d\n",
        meta_data.number_of_samples,
        meta_data.x_spacing,
        meta_data.usrp_rf_rate,
        meta_data.usrp_fcenter
    );

    


    // XXX: Define other parameters
    double restricted_frequencies[] = { 0,0 };
    double clear_freq_range[] = { 25 * pow(10,3), 40 * pow(10,3) };
    double beam_angle = 0.08482300164692443; //0.08482300164692443
    double smsep = 1 / (2 * 250 * pow(10, 3)); // ~4 ms


    // Call calc_clear_freq_on_raw_samples
    calc_clear_freq_on_raw_samples(
        raw_samples, &meta_data, restricted_frequencies, 
        clear_freq_range, beam_angle, smsep);
    
    // Free allocated memory
    fftw_free(raw_samples);
    free(meta_data.antenna_list);
    
    // Print processing time; Stopwatch End
    // double t2 = dsecnd();
    // printf("Time for Clear Freq Search: %lf", (t2-t1));

    return 0;
};
