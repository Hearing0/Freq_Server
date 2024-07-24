#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>      // FFT transform library
#include <cblas.h>      // Matrix Multi library
// #include <utils/phasing_utils.c>
#include <string.h>
#include "../utils/ini_parser/ini.c"
#include "../utils/read_config.c"


// Define Constants
#define PI 3.14159265358979323846
double C = 3e8;//* pow(10, 8.0);
bool verbose = true;
bool test_samples = false;
// #define RESTRICT_FILE = '/home/radar/repos/SuperDARN_MSI_ROS/linux/home/radar/ros.3.6/tables/superdarn/site/site.sps/restrict.dat.inst'


/*
* NOTE: 
* When compiling make sure to add "-lm" and "-lfftw3" to link the libraries 
* to the compilier. 
* 
* TODO: Find GCC optimize flags
*/

// TODO: Find GCC optimization flags



typedef struct {
    int *antenna_list;
    int num_antennas;
    int number_of_samples;
    double x_spacing;
    int usrp_rf_rate;
    int usrp_fcenter;
} sample_meta_data;

typedef struct {
    double *restricted_freq;
    double *clear_freq_range;
} freq_data;

typedef struct {
    double noise;
    double tfreq;
} clear_freq;

/**
 * @brief  Reads intial Clear Freq parameters from a converted txt file (see utils/pickle_text_convert.py)
 * @note   
 * @param  *filename: Filepath of converted txt file
 * @param  *meta_data: Meta data for current sample batch
 * @param  **clear_freq_range: Range for the Clear Freq
 * @param  ***raw_samples: 14x2500 complex sample array
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
        if (strncmp(line, "clear_freq_range:", 15) == 0) {
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
            meta_data->num_antennas = antenna_list_size;
            continue;
        }

        // Raw Samples
        if (strncmp(line, "raw_samples:", 12) == 0 && test_samples) {
            printf("[Clear Freq Search] Aquiring test samples from pickle files...");
            // Allocate mem
            *raw_samples = (fftw_complex **)fftw_malloc(meta_data->num_antennas * sizeof(fftw_complex *));
            for (int i = 0; i < meta_data->num_antennas; i++) {
                (*raw_samples)[i] = (fftw_complex *)fftw_malloc(meta_data->number_of_samples * sizeof(fftw_complex));
            }
            if (*raw_samples == NULL) {
                perror("Error allocating memory for raw samples");
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

/**
 * @brief  Loads in the beam configuration from array_config.ini. 
 * @note   By DF
 * @param  *n_beams:    Number of beams
 * @param  *beam_sep:   Angle Offset between beams (in degrees)
 * @retval None
 */
void load_beam_config(double *x_spacing, int *n_beams, double *beam_sep){
    const char *config_path = "../Freq_Server/utils/clear_freq_input/array_config.ini";
    Config config;

    if (ini_parse(config_path, config_ini_handler, &config) < 0) {
        printf("Can't load 'config.ini'\n");
        return;
    }

    *x_spacing = config.array_info.x_spacing;
    *n_beams = config.array_info.nbeams;
    *beam_sep = config.array_info.beam_sep;
}

/**
 * @brief  Writes Frequency Spectrum to csv file to be plotted in python.
 * @note   By DF
 * @param  *filename: 
 * @param  *spectrum: 
 * @param  num_samples: 
 * @retval None
 */
void write_spectrum_to_csv(const char *filename, fftw_complex *spectrum, double *freq_vector, int num_samples) {
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
 * @brief  Calculates Beam Azimuth Angle.
 * @note   By DF
 * @param  n_beams:     Total number of beams
 * @param  beam_num:    Beam number for calculating Beam Azimuth
 * @param  beam_sep:    Beam Seperation (in degrees)
 * @retval Returns the Beam Azumuth Angle (in radians)
 */
double calc_beam_angle(int n_beams, int beam_num, double beam_sep) {
    // Calculate Beamforming shift
    double center_beam = (n_beams - 1) / 2;

    // Calculate Beam Azimuth
    double b_azi = ((beam_num - center_beam) * beam_sep) * (PI / 180);
    if (verbose){
        printf("n_beams: %d, beam_num: %d, beam_sep: %lf\n", n_beams, beam_num, beam_sep);
        printf("    beam = %lf degree", (b_azi * 180 / PI));
    }
    return b_azi;
}

/**
 * @brief  Calculates phase increment between antennas to produce a mainlobe sterring 
 *         of beam_angle at center_frequency.
 * @note   By DF
 * @param  beam_angle:          Distrubution of the beam (in radians)
 * @param  center_frequency:    Frequency at the center of the beam; 
 * * used to phase-shift allign the other frequencys (in Hz)
 * @param  x_spacing:           Spacing in-between antennas
 *         
 * @retval Returns the Phase Shift (in degrees)
 */
float calc_phase_increment(float beam_angle, double center_frequency, double x_spacing) {
    double wavelength = C / center_frequency;
    double phase_shift = (2 * PI * x_spacing * sin(beam_angle)) / wavelength;
    if (verbose) {
        printf("center_freq: %lf\nx_spacing: %lf\nphase_shift: %lf degree\n", center_frequency, x_spacing, phase_shift * 180 / PI);
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


clear_freq find_clear_freq(fftw_complex *spectrum_power, double *freq_vector, double start_freq, double end_freq, double clear_bw, double *tfreq, double *noise) {

}

// HACK apply efficient matrix multi via cblas_dgemm
clear_freq calc_clear_freq_on_raw_samples(fftw_complex **raw_samples, sample_meta_data *meta_data, double *restricted_frequencies, double *clear_freq_range, double beam_angle, double smsep) {
    char *spectrum_file = "spectrum_output.csv";
    
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
        perror("Error allocating memory.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate and Apply phasing vector
    float phase_increment = calc_phase_increment(beam_angle, (clear_freq_range[0] + clear_freq_range[1]) / 2, meta_data->x_spacing);
    printf("phase_increment: %lf\n", phase_increment);
    // phase_increment = .000738
    
    for (int i = 0; i < meta_data->num_antennas; i++) {
        if (i <= 15 || i >= 20) {
            if (i < meta_data->num_antennas) {
                // printf("antenna[%d]: %d\n", i, antennas[i]);
                phasing_vector[i] = rad_to_rect(antennas[i] * phase_increment);
                // printf("phase_vector[%d]: %f + %fi\n", i, creal(phasing_vector[i]), cimag(phasing_vector[i]));
            } else {
                fprintf(stderr, "Error: Accessing antennas out of bounds at index %d\n", i);
                exit(EXIT_FAILURE);
            }
        }
    }
    printf("antenna[13]: %d\n", antennas[13]);
    // Antennas[13] = 18
    printf("phasing_vector[13]: %f + %fi\n", creal(phasing_vector[13]), cimag(phasing_vector[13]));
    // phasing_vector[13] = 0.999912 + 0.013288i

    // Apply beamforming
    for (int i = 0; i < num_samples; i++) {
        double real_sum = 0.0;
        double imag_sum = 0.0;
        
        for (int aidx = 0; aidx < meta_data->num_antennas; aidx++) {
            double real_sample = creal(raw_samples[aidx][i]);       
            double imag_sample = cimag(raw_samples[aidx][i]);
            double real_phase = creal(phasing_vector[aidx]);
            double imag_phase = cimag(phasing_vector[aidx]);
            if (i == 2499) {
                printf("sample[%d][2499]    = %f + %fi\n", aidx, real_sample, imag_sample);
                printf("phase[%d]           = %f + %fi\n", aidx, real_phase, imag_phase);
            }

            real_sum += real_sample * real_phase - imag_sample * imag_phase;
            imag_sum += real_sample * imag_phase + imag_sample * real_phase;
            // printf("beamformed_samples[%d]: %f + %fi\n", i, creal(beamformed_samples[i]), cimag(beamformed_samples[i]));
        }

        beamformed_samples[i] = real_sum + I * imag_sum;
    }
    printf("beamformed[2499]    = %f + %fi\n", creal(beamformed_samples[2499]), cimag(beamformed_samples[2499]));
    // beamformed[2499] =  70.466100 + -169.215264i
    // TODO: Verify answer


    // Spectral Estimation
    fft_clrfreq_samples(beamformed_samples, num_samples, spectrum_power);

    // Print the spectrum power
    // for (int i = 0; i < num_samples; i++) {
    //     printf("Frequency bin %d: %f + %fi\n", i, creal(spectrum_power[i]), cimag(spectrum_power[i]));
    // }
    


    // Frequency Vector Calculation
    double delta_f = meta_data->usrp_rf_rate / num_samples;

    for (int i = 0; i < num_samples; i++) {
        freq_vector[i] = i * delta_f - (meta_data->usrp_rf_rate / 2) + meta_data->usrp_fcenter * 1000;
    }

    write_spectrum_to_csv(spectrum_file, spectrum_power, freq_vector, num_samples);


    /// END of Spectrum Calc

    // // Mask restricted frequencies
    // // 

    // // Find clear frequency
    // double clear_bw = 2e6 / smsep;
    // double tfreq, noise;
    // clear_freq cfreq = find_clear_freq(spectrum_power, freq_vector, clear_freq_range[0] * 1e3, clear_freq_range[1] * 1e3, clear_bw, &tfreq, &noise);

    // // Output results
    // printf("Clear Frequency: %f, Noise: %f\n", tfreq, noise);

    printf("Finished Clear Freq!\n");
    
    fftw_free(phasing_vector);
    fftw_free(beamformed_samples);
    fftw_free(spectrum_power);
    free(freq_vector);

    // return cfreq;
}

void init() {
    
}

// int main() {
clear_freq clear_freq_search(fftw_complex **raw_samples) {
    // Stopwatch Start
    // double t1 = dsecnd();

    // HACK: Setup file_path environment variable
    const char *input_file_path = "../Freq_Server/utils/clear_freq_input/clrfreq_dump.1.txt";
    // const char *output_file_path = "../utils/txt_output/result.txt";

    // Initial Data Variables
    sample_meta_data meta_data = {0};
    // fftw_complex **raw_samples = NULL;
    int n_beams, beam_num;
    double beam_sep;
    freq_data freq_data;

    // Load in data for Clear Freq Calculation
    read_input_data(input_file_path, &meta_data, &freq_data.clear_freq_range, &raw_samples);

    // Load in data for Beam Angle Calculation
    load_beam_config(&meta_data.x_spacing, &n_beams, &beam_sep);
    beam_num = 1;

    printf("num_samples: %d\nx_spacing: %lf\nusrp_rf_rate: %d\nusrp_fcenter: %d\n",
        meta_data.number_of_samples,
        meta_data.x_spacing,
        meta_data.usrp_rf_rate,
        meta_data.usrp_fcenter
    );       
    printf("n_beams: %d\nbeam_sep: %f\nbeam_num: %d\n", n_beams, beam_sep, beam_num);

    // Check last sample
    // fftw_complex sample = raw_samples[meta_data.num_antennas - 1][meta_data.number_of_samples - 1];
    // printf("raw_samples[%d][%d]: %f + %fi\n", meta_data.num_antennas - 1, meta_data.number_of_samples - 1, creal(sample), cimag(sample));
    // Should be raw_samples[13][2499]: -134.000000 + 168.000000i


    // XXX: Define other parameters
    double restricted_frequencies[] = { 0,0 };
    double clear_freq_range[] = { 25 * pow(10,3), 40 * pow(10,3) };
    // double beam_angle = calc_beam_angle(n_beams, beam_num, beam_sep);    // 
    double beam_angle = 0.08482300164692443;        // in radians
    double smsep = 1 / (2 * 250 * pow(10, 3));      // ~4 ms

    // Find Clear Frequency
    clear_freq cfreq = calc_clear_freq_on_raw_samples(
        raw_samples, &meta_data, restricted_frequencies, 
        clear_freq_range, beam_angle, smsep);
    
    // Free allocated memory
    free(meta_data.antenna_list);  

    // Print processing time; Stopwatch End
    // double t2 = dsecnd();
    // printf("Time for Clear Freq Search: %lf", (t2-t1));

    // return 1;
    return cfreq;
};
