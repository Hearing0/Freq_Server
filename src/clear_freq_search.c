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


#ifndef CLK_TCK
#define CLK_TCK 60
#endif


// Define Constants
#define PI 3.14159265358979323846
double C = 3e8;//* pow(10, 8.0);
bool verbose = true;
bool test_samples = true;
int CLR_BANDS_MAX = 6;
int CLRFREQ_RES = 2e3;                      // XXX: Is CLRFREQ_RES just delta_f? Both are 2kHz.

// #define RESTRICT_FILE = '/home/radar/repos/SuperDARN_MSI_ROS/linux/home/radar/ros.3.6/tables/superdarn/site/site.sps/restrict.dat.inst'


/*
* NOTE: 
* When compiling make sure to add "-lm" and "-lfftw3" to link the libraries 
* to the compilier. 
*/

// TODO: Find GCC optimization flags



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
            meta_data->num_antennas = 2; //antenna_list_size;
            continue;
        }

        // Raw Samples
        if (strncmp(line, "raw_samples:", 12) == 0 && test_samples) {
            printf("[Clear Freq Search] Aquiring test samples from pickle files...\n");
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
 * @param  *filename:       The filepath for the saved CSV file
 * @param  *spectrum:       Frequency Spectrum
 * @param  *freq_vector:    Array for Frequency (steps of frequency per sample)
 * @param  num_samples:     Number of samples collected per antenna
 * @retval None
 */
void write_spectrum_to_csv(char *filename, fftw_complex *spectrum, double *freq_vector, int num_samples) {
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
 * @brief  Writes the Real/Imaginary magnitude to csv file to be plotted in python.
 * @note   By DF
 * @param  *filename:           The filepath for the saved CSV file
 * @param  *raw_samples_mag:    Int array of the Real/Imaginary magnitude
 * @param  *freq_vector:        Array for Frequency (steps of frequency per sample)
 * @param  num_samples:         Number of samples collected per antenna
 * @retval None
 */
void write_sample_mag_to_csv(char *filename, int **raw_samples_mag, double *freq_vector, sample_meta_data *meta_data) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Samples,Power\n");
    for (int j = 0; j < meta_data->num_antennas; j++)
        for (int i = 0; i < meta_data->number_of_samples; i++) {
            fprintf(file, "%f,%d\n", freq_vector[i], raw_samples_mag[j][i]);
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
        printf("search_center_freq: %lf\nx_spacing: %lf\nphase_shift: %lf degree\n", center_frequency, x_spacing, phase_shift * 180 / PI);
    }
    return phase_shift; 
}

/**
 * Converts radians to rectangular form
 */
double complex rad_to_rect(double phase) {
    return cexp(I * phase);
}

// HACK Check time vs regular fft; plans are expensive to create.
/**
 * @brief  Fast Fourier Transform (FFT) conversion from beamformed_samples 
 * * into a spectrum.
 * @note    By DF
 * @param  *samples: Array of samples that have already been 
 * * prepared for FFT
 * @param  number_of_samples: Number of samples
 * @param  *spectrum_power: Output array for the resultant spectrum
 * @retval None
 */
void fft_clrfreq_samples(fftw_complex *samples, int num_samples, fftw_complex *spectrum_power) {
    fftw_plan plan = fftw_plan_dft_1d(num_samples, samples, spectrum_power, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * @brief  Convolves 'u' array with 'v' array.
 * @note   
 * @param  u: 1st Array
 * @param  u_size: Size of array 'u' 
 * @param  v: 2nd Array 
 * @param  v_size: Size of array 'v'
 * @param  *result: Array for Convolution Result
 * @retval N/A
 */
void convolve(double* u, int u_size, int* v, int v_size, double* result) {  
    for (int i = 0; i < u_size - v_size; i++) {
        result[i] = 0;
        for (int j = 0; j < v_size; j++) {
            result[i] += v[j] * u[i + j];
        }
    }
}

// TODO: Reformat finding Insertion
// TODO: Test by Compile
// TODO: Parse radar_config_constants.py for CLRFREQ_RES
// TODO: Strike balance b/w speed and time using clear_sample_bw
// TODO: Test with -20dBm 10MHz signal
// TODO: Compare Time Complexity w/ Py Version
// TODO: Obtain f_start, f_end from RESTRICT file
// Try convolve with filter then find min of convolve
//      gnu or intel scientific library
/**
 * @brief  Processes Spectrum data to find the lowest noise frequency bands 
 * * (returned as clr_freq_bands). Steps are as follows: 
 * * (1) Setup Freq Search parameter data,
 * * (2) Scan Search Range with Bandpass Filter by way of convolution, 
 * * (3) Find Noise of BPF and compare w/ current Clear Freq Bands, and
 * * (4) If appropriate spot found, insert BPF's Range as New Clear Freq Band.
 * * * Overwriting any pre-existing, Intersecting Clear Freq Bands as necesary.  
 * @note   By DF
 * @param  *spectrum_power: Spectrum Data (Power per Sample)
 * @param  meta_data: Misc info on operating Radar parameters  
 * @param  delta_f: Frequency step per Sample 
 * @param  f_start: Clear Freq Search Bound start
 * @param  f_end: Clear Freq Search Bound end
 * @param  clear_bw: Bandwidth of the Clear Frequency Bands. If 0, default to 40kHz.
 * @param  *lowest_freq_bands: Passed by reference; Overwritten with an array of the 
 * * lowest noise freq_bands.  
 * @retval None
 */
void find_clear_freqs(fftw_complex *spectrum_power, sample_meta_data meta_data, double delta_f, double f_start, double f_end, int clear_bw, freq_band *clr_freq_bands) {
    printf("[find_clear_freqs()] Entered find_clear_freqs()...\n");
    if (clear_bw == 0) clear_bw = 40e3;
    int clear_sample_bw = clear_bw / CLRFREQ_RES; 
    // avg spectrum_power using ratio of delta_f to CLRFREQ_RES 

    // Define Range of Clear Freq Search 
    int spectrum_sample_start = (int) ((meta_data.usrp_fcenter * 1000 - meta_data.usrp_rf_rate / 2) / delta_f);
    int clr_search_sample_start = (int) (f_start / delta_f) - spectrum_sample_start;
    int clr_search_sample_end = (int) (f_end / delta_f) - spectrum_sample_start;

    // Trim Spectrum Data to only Clear Search Range (Used for convolving)
    int clr_search_sample_bw = clr_search_sample_end - clr_search_sample_start;
    double clr_search_band[clr_search_sample_bw];
    memset(clr_search_band, 0, clr_search_sample_bw * sizeof(double));
    for (int i = clr_search_sample_start; i < clr_search_sample_end; i++) {
        clr_search_band[i - clr_search_sample_start] = sqrt(creal(spectrum_power[i]) * creal(spectrum_power[i]) + cimag(spectrum_power[i]) * cimag(spectrum_power[i]));
        if (i < 2 == 0 || i > clr_search_sample_end - 2) {
            printf("clr_search_band[%d]: %f\n", i - clr_search_sample_start, clr_search_band[i - clr_search_sample_start]);
            printf("                   : %f + i%f\n", creal(spectrum_power[i]), cimag(spectrum_power[i]));
        }
    }

    // Scan Search range w/ Bandpass Filter (BPF) to find Clear Freq Band
    printf("[find_clear_freqs()] Scanning Search Range w/ Bandpass...\n");
    int *bpf = (int *) malloc(sizeof(int) * clear_sample_bw);
    for (int band_i = 0; band_i < clear_sample_bw; band_i++) {
        bpf[band_i] = 1;
    }
    printf("    clr_search_sample_start: %d\n    clr_search_sample_end: %d\n    clr_search_sample_bw: %d\n", 
        clr_search_sample_start, clr_search_sample_end, clr_search_sample_bw);

    // Convolve BPF with Search Range
    double *convolve_result = calloc(clr_search_sample_bw, sizeof(double));
    convolve(clr_search_band, clr_search_sample_bw, bpf, clear_sample_bw, convolve_result);
    printf("    Convolved Scan Band and BPF...\n");

    // for (int i = 0; i < clr_search_sample_bw; i++)
    // {
    //     printf("convolve[%d]: %f", i, convolve_result[i]);
    // }
    
    // Initialize Clear Freq Bands
    freq_band *clr_bands = clr_freq_bands;
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        clr_bands[i].f_start = clr_search_sample_start * delta_f - (meta_data.usrp_rf_rate / 2) + meta_data.usrp_fcenter * 1000;
        clr_bands[i].f_end = clr_search_sample_end * delta_f - (meta_data.usrp_rf_rate / 2) + meta_data.usrp_fcenter * 1000;
        clr_bands[i].noise = RAND_MAX;
    };
    int min_idx[CLR_BANDS_MAX];
    
    // Identify lowest noise bands from convolve results...
    freq_band curr_band;
    for (int i = 0; i < clr_search_sample_bw; i++) {
        curr_band.f_start = (spectrum_sample_start + i) * delta_f;
        curr_band.f_end = (spectrum_sample_start + i + clear_sample_bw) * delta_f;
        curr_band.noise = convolve_result[i];
        printf("[%d] | %d -- %f -- %d|\n", i, curr_band.f_start, curr_band.noise, curr_band.f_end);
        // printf("curr noise[%d]: %f\n", i, convolve_result[i]);
        
        int insert_idx = -1;
        int intersect_idx = -1;
        // Find Insertion spot in clr_freq_bands
        // Compare curr power with min_powers...
        for (int j = CLR_BANDS_MAX - 1; j >= 0 ; j--) {
            // Update Insert Index; maintaining ascending order 
            if (curr_band.noise < clr_bands[j].noise && curr_band.noise > 0) {
                insert_idx = j;
            }
            // Check for Intersecting Band; get intersecting clr_band index
            if ( intersect_idx == -1 && 
                ((clr_bands[j].f_start < curr_band.f_start && curr_band.f_start < clr_bands[j].f_end) ||
                    (clr_bands[j].f_start < curr_band.f_end && curr_band.f_end < clr_bands[j].f_end))) {
                intersect_idx = j;
            }
            // Continue Intersection Search 
        }
        printf("    Intersection Search finished...\n");
        printf("    intersect_idx: %d\n   insert_idx: %d\n", intersect_idx, insert_idx);

        // Insertion Point was Found...
        if (insert_idx != -1) {
            // Intersection w/ curr_band was also Found...
            if (intersect_idx != -1) {
                // Special: If Intersect is less noisy, do not place 
                if (insert_idx > intersect_idx) continue;
                printf("    Intersecting Insertion found w/...\n");
                freq_band inter_band = clr_bands[intersect_idx];
                printf("        i-band = | %d -- %f -- %d|\n", inter_band.f_start, inter_band.noise, inter_band.f_end);

                // Special: Shift right till the Intersecting band is overwritten 
                if (insert_idx < intersect_idx) {
                    // Debug: verify shifting @ clr sample 210
                    if (i == 210) for (int j = 0; j < CLR_BANDS_MAX; j++) {
                        printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", j, clr_bands[j].f_start, clr_bands[j].noise, clr_bands[j].f_end);
                    }
                    
                    printf("        shifting clr_bands for intersect...\n");
                    for (int j = intersect_idx - 1; j >= insert_idx; j--) {
                        if (j + 1 < CLR_BANDS_MAX) {
                            clr_bands[j + 1] = clr_bands[j];
                            min_idx[j + 1] = min_idx[j];
                        }
                    }
                }
            } 
            // Only Insertion Point Found...
            else {
                printf("    Insertion found...\n");
                // Special: Keep pre-existing bands by shifting them to right
                for (int j = CLR_BANDS_MAX - 2; j >= insert_idx; j--) {
                    printf("        shifting clr_bands for insert...\n");
                    if (j + 1 < CLR_BANDS_MAX) {
                        clr_bands[j + 1] = clr_bands[j];
                        min_idx[j + 1] = min_idx[j];
                    }
                }
            }

            // Insert curr_band and store its index
            clr_bands[insert_idx] = curr_band;
            min_idx[insert_idx] = i;

            // Debug: verify shifting @ clr sample 210  
            if (i == 210) for (int j = 0; j < CLR_BANDS_MAX; j++) {
                printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", j, clr_bands[j].f_start, clr_bands[j].noise, clr_bands[j].f_end);
            }
        }
    }
    printf("    Convolution Results packed up...\n");

    // Debug: Output results
    // for (int i = 0; i < CLR_BANDS_MAX; i++)
    //     printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", i, clr_bands[i].f_start, clr_bands[i].noise, clr_bands[i].f_end);

    // Free allocated memory
    free(convolve_result);
    free(bpf);

    printf("[find_clear_freqs()] Exiting find_clear_freqs()...\n");
}


// HACK apply efficient matrix multi via cblas_dgemm
void calc_clear_freq_on_raw_samples(fftw_complex **raw_samples, sample_meta_data *meta_data, double *restricted_frequencies, double *clear_freq_range, double beam_angle, double smsep, freq_band *clr_bands) {
    char *spectrum_file = "../Freq_Server/utils/csv_dump/spectrum_output.csv";
    char *sample_re_file = "../Freq_Server/utils/csv_dump/samples/sample_re_output.csv";
    char *sample_im_file = "../Freq_Server/utils/csv_dump/samples/sample_im_output.csv";
    int **sample_re = NULL;
    int **sample_im = NULL;
    
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
    sample_im = (int**) malloc(sizeof(int*) * meta_data->num_antennas);
    sample_re = (int**) malloc(sizeof(int*) * meta_data->num_antennas);
    for (int i = 0; i < meta_data->num_antennas; i++) {
        sample_im[i] = (int *)malloc(meta_data->number_of_samples * sizeof(int));
        sample_re[i] = (int *)malloc(meta_data->number_of_samples * sizeof(int));
    }   
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
    printf("phasing_vector[13]: %f + %fi\n", creal(phasing_vector[13]), cimag(phasing_vector[13]));
    // Antennas[13] = 18
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
            if (verbose && i == 2499) {
                printf("sample[%d][2499]    = %f + %fi\n", aidx, real_sample, imag_sample);
                printf("phase[%d]           = %f + %fi\n", aidx, real_phase, imag_phase);
            }

            real_sum += real_sample * real_phase - imag_sample * imag_phase;
            imag_sum += real_sample * imag_phase + imag_sample * real_phase;
            // printf("beamformed_samples[%d]: %f + %fi\n", i, creal(beamformed_samples[i]), cimag(beamformed_samples[i]));

            sample_im[aidx][i] = cimag(raw_samples[aidx][i]);
            sample_re[aidx][i] = creal(raw_samples[aidx][i]);
        }
        beamformed_samples[i] = real_sum + I * imag_sum;
    }
    printf("beamformed[2499]    = %f + %fi\n", creal(beamformed_samples[2499]), cimag(beamformed_samples[2499]));
    // beamformed[2499] =  70.466100 + -169.215264i

    // TODO: Plot usrp/antenna im and re separately (4 plots) 
    //          to confirm that one of them is correctly all zeros for im

    // Spectrum Averaging; delinate Transmitters and filter out noise
    //      Split time series into blocks, then average them to filter out noise
    // int num_of_avg = 500;
    // if ((num_samples % num_of_avg) != 0) num_of_avg = 100;
    // int samples_per_avg = num_samples / num_of_avg; 
    // fftw_complex *avg_samples = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples_per_avg);
    // printf("[Spectral Avg] Number of averages: %d\n", num_of_avg);
    // for (int i = 0; i < num_samples; i++) {
    //     avg_samples[i % samples_per_avg] += beamformed_samples[i];
    // }

    // Spectrum Calculation
    fft_clrfreq_samples(beamformed_samples, num_samples, spectrum_power);

    // Display Spectrum Power
    // for (int i = 0; i < meta_data.number_of_samples; i++) {
    //     printf("Frequency bin %d: %f + %fi\n", i, creal(spectrum_power[i]), cimag(spectrum_power[i]));
    // }


    // Frequency Vector Calculation
    double delta_f = meta_data->usrp_rf_rate / num_samples;
    double f_start = meta_data->usrp_fcenter * 1000 - (meta_data->usrp_rf_rate / 2);
    for (int i = 0; i < num_samples; i++) {
        freq_vector[i] = i * delta_f + f_start;
    }

    printf("delta_f: %f\nnum_samples: %d\nfcenter: %d\n", delta_f, num_samples, meta_data->usrp_fcenter * 1000);

    // Save data to csv
    write_spectrum_to_csv(spectrum_file, spectrum_power, freq_vector, num_samples);
    write_sample_mag_to_csv(sample_im_file, sample_im, freq_vector, meta_data);
    write_sample_mag_to_csv(sample_re_file, sample_re, freq_vector, meta_data);

    /// END of Spectrum Calc


    // Mask restricted frequencies

    // int clr_search_sample_bw = (int) (round(clear_freq_range[1] / delta_f) - round(clear_freq_range[0] / delta_f));
    int clear_sample_start = (int) round((clear_freq_range[0] - f_start) / delta_f);
    int clear_sample_end = (int) round((clear_freq_range[1] - f_start) / delta_f);
    for (int i = clear_sample_start; i < clear_sample_end; i++) {
        if (i < 2 || i > clear_sample_end - 2) printf("spectrum_pow[%d]: %f + i%f\n", i, creal(spectrum_power[i]), cimag(spectrum_power[i]));   
    }


    // Find clear frequency
    double clear_bw = 2e6 / smsep;
    clear_bw = 0;
    clock_t t1, t2;
    t1 = clock();
    find_clear_freqs(spectrum_power, *meta_data, delta_f, clear_freq_range[0], clear_freq_range[1], clear_bw, clr_bands);
    t2 = clock();
    printf("clear_freq_search (mili sec): %lf\n", ((double) (t2 - t1)) / (CLOCKS_PER_SEC));

    // Debug: Output results
    for (int i = 0; i < CLR_BANDS_MAX; i++)
        printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", i, clr_bands[i].f_start, clr_bands[i].noise, clr_bands[i].f_end);
    
    printf("Finished Clear Freq Search!\n");
    
    fftw_free(phasing_vector);
    fftw_free(beamformed_samples);
    fftw_free(spectrum_power);
    free(freq_vector);
}

// XXX: Add a initialization???

// int main() {
clear_freq clear_freq_search(fftw_complex **raw_samples, freq_band *clr_bands) {

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
    double clear_freq_range[] = { 12 * pow(10,6), 12.5 * pow(10,6) };
    // double beam_angle = calc_beam_angle(n_beams, beam_num, beam_sep);  
    double beam_angle = 0.08482300164692443;        // in radians
    double smsep = .004; // 1 / (2 * 250 * pow(10, 3));      // ~4 ms

    // Stopwatch Start
    double t1,t2;
    t1 = clock();

    // Find Clear Frequency Bands
    calc_clear_freq_on_raw_samples(
        raw_samples, &meta_data, restricted_frequencies, 
        clear_freq_range, beam_angle, smsep, clr_bands);
    
    // Print processing time; Stopwatch End
    t2 = clock();
    printf("clear_freq_search (mili sec): %lf", ((double) (t2 - t1)) / (CLOCKS_PER_SEC * 1000));
    
    // Free allocated memory
    free(meta_data.antenna_list);  


};