#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>      // FFT transform library
// #include <cblas.h>      // Matrix Multi library
// #include <utils/phasing_utils.c>
#include <string.h>
#include <time.h>


#ifndef CLK_TCK
#define CLK_TCK 60
#endif


// Define Constants
#define IDX_LAST_IA 19     // Last Interferrometer Array
#define IDX_LAST_MA 15     // Last Main Array
#define PI 3.14159265358979323846
#define C  3e8
#define CLR_BANDS_MAX 6
#define CLRFREQ_RES 2e3          

// Debug Flags
#define VERBOSE 1
#define SPECTRAL_AVGING 1
#define TEST_SAMPLES 1
#define TEST_CLR_RANGE 1

#include "../utils/misc_read_writes.c" // Pulls Debug Testing Data

// TODO: Pass in clr_freq_range via restrict actual file
// #define RESTRICT_FILE = '/home/radar/repos/SuperDARN_MSI_ROS/linux/home/radar/ros.3.6/tables/superdarn/site/site.sps/restrict.dat.inst'

// TODO: Pass in x_spacing, etc. via config actual file

// TODO: Store Clear Freq Bands for each Beam 

/*
* NOTE: 
* When compiling make sure to add "-lm" and "-lfftw3" to link the libraries 
* to the compilier. 
*/

// TODO: Find GCC optimization flags


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
    if (VERBOSE){
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
    if (VERBOSE) {
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

/**
 * @brief  Fast Fourier Transform (FFT) conversion from Time four_spectrums  
 * * into Frequency spectrum.
 * @note    By DF
 * @param  *four_spectrums: Array of four_spectrums that have already been 
 * * prepared for FFT
 * @param  number_of_samples: Number of four_spectrums
 * @param  *spectrum: Output array for the resultant spectrum
 * @retval None
 * 
 * @deprecated At the time of writing, each component of the function is used 
 * * separately for fft averaging
 */
void fft_samples(fftw_complex *four_spectrums, int num_samples, fftw_complex *spectrum) {
    fftw_plan plan = fftw_plan_dft_1d(num_samples, four_spectrums, spectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * @brief  Convolves 'u' array with 'v' array.
 * @note   By DF and B. Bristrow
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
            result[i] += v[j] * u[i + j]; // result[i + (int) v_size] += ...
        }
    }
}

// TODO: Reconfigure logic to find most probable signal
// TODO: Test with -20dBm 10MHz signal
// TODO: Obtain f_start, f_end from RESTRICT file
// TODO: Parse radar_config_constants.py for CLRFREQ_RES
// TODO: Strike balance b/w speed and time using clear_sample_bw
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
 * @param  *spectrum: Spectrum Data (Power per Sample)
 * @param  meta_data: Misc info on operating Radar parameters  
 * @param  delta_f: Frequency step per Sample 
 * @param  f_start: Clear Freq Search Bound start
 * @param  f_end: Clear Freq Search Bound end
 * @param  clear_bw: Bandwidth of the Clear Frequency Bands. If 0, default to 40kHz.
 * @param  *lowest_freq_bands: Passed by reference; Overwritten with an array of the 
 * * lowest noise freq_bands.  
 * @retval None
 */
void find_clear_freqs(double *spectrum, sample_meta_data meta_data, double delta_f, double f_start, double f_end, int clear_bw, freq_band *clr_freq_bands) {
    printf("[find_clear_freqs()] Entered find_clear_freqs()...\n");
    if (clear_bw == 0) clear_bw = 10e3;
    int clear_sample_bw = clear_bw / CLRFREQ_RES; 

    // int guard_band = ; // twice the bw of transmitted signal

    // Define Range of Clear Freq Search 
    int spectrum_sample_start = (int) ((meta_data.usrp_fcenter * 1000 - meta_data.usrp_rf_rate / 2) / delta_f);
    int clr_search_sample_start = (int) (f_start / delta_f) - spectrum_sample_start;
    int clr_search_sample_end = (int) (f_end / delta_f) - spectrum_sample_start;

    // Trim Spectrum Data to only Clear Search Range (Used for convolving)
    int clr_search_sample_bw = clr_search_sample_end - clr_search_sample_start;
    double clr_search_band[clr_search_sample_bw];
    memset(clr_search_band, 0, clr_search_sample_bw * sizeof(double));
    for (int i = 0; i < clr_search_sample_bw; i++) {
        clr_search_band[i] = spectrum[i + clr_search_sample_start];

        // Debug: Check the Clear Search Range 
        if (i < 2  || i > clr_search_sample_bw - 2) {
            printf("clr_search_band[%d]: %f\n", i, clr_search_band[i]);
            printf("                   : %f \n", spectrum[i + clr_search_sample_start]);
        }
    }

    // Scan Search range w/ Bandpass Filter (BPF) to find Clear Freq Band
    printf("[find_clear_freqs()] Scanning Search Range w/ Bandpass...\n");
    int *bpf = (int *) malloc(sizeof(int) * clear_sample_bw);
    for (int band_i = 0; band_i < clear_sample_bw; band_i++) {
        bpf[band_i] = 1;
    }
    if (VERBOSE) printf("    clr_search_sample_start: %d\n    clr_search_sample_end: %d\n    clr_search_sample_bw: %d\n", 
        clr_search_sample_start, clr_search_sample_end, clr_search_sample_bw);

    // Convolve BPF with Search Range
    int convolve_bw = clr_search_sample_bw - clear_sample_bw;
    double *convolve_result = calloc(convolve_bw, sizeof(double));
    convolve(clr_search_band, clr_search_sample_bw, bpf, clear_sample_bw, convolve_result);
    printf("    Convolved Scan Band and BPF...\n");

    // Debug: Check convolve result
    // for (int i = 0; i < clr_search_sample_bw; i++)
    // {
    //     printf("convolve[%d]: %f\n", i, convolve_result[i]);
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
    for (int i = 0; i < convolve_bw; i++) {
        curr_band.f_start = (spectrum_sample_start + i) * delta_f;
        curr_band.f_end = (spectrum_sample_start + i + clear_sample_bw) * delta_f;
        curr_band.noise = convolve_result[i];
        printf("[%d] | %d -- %f -- %d|\n", i, curr_band.f_start, curr_band.noise, curr_band.f_end);
        
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
                // Special: If Intersect is less noisy, do not place/skip
                if (insert_idx > intersect_idx) continue;
                printf("    Intersecting Insertion found w/...\n");
                freq_band inter_band = clr_bands[intersect_idx];

                printf("        i-band = | %d -- %f -- %d|\n", inter_band.f_start, inter_band.noise, inter_band.f_end);

                // Special: Shift right till the Intersecting band is overwritten 
                if (insert_idx < intersect_idx) {
                    // Debug: verify bands shift properly @ clr sample 210
                    // if (i == 210) for (int j = 0; j < CLR_BANDS_MAX; j++) {
                    //     printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", j, clr_bands[j].f_start, clr_bands[j].noise, clr_bands[j].f_end);
                    // }
                    
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

            // Insert curr_band and store its sample index
            clr_bands[insert_idx] = curr_band;
            min_idx[insert_idx] = i;

            // Debug: verify shifting @ clr sample 210  
            // if (i == 210) for (int j = 0; j < CLR_BANDS_MAX; j++) {
            //     printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", j, clr_bands[j].f_start, clr_bands[j].noise, clr_bands[j].f_end);
            // }
        }
    }

    printf("    Convolution Results packed up...\n");

    // Debug: Output Final Clear Freq Bands
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
    char *clr_freq_file = "../Freq_Server/utils/csv_dump/clr_freq_output.csv";
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
    double *freq_vector = (double*) malloc(sizeof(double) * num_samples);
    sample_im = (int**) malloc(sizeof(int*) * meta_data->num_antennas);
    sample_re = (int**) malloc(sizeof(int*) * meta_data->num_antennas);
    for (int i = 0; i < meta_data->num_antennas; i++) {
        sample_im[i] = (int *)malloc(meta_data->number_of_samples * sizeof(int));
        sample_re[i] = (int *)malloc(meta_data->number_of_samples * sizeof(int));
    }   
    if (!phasing_vector || !beamformed_samples || !freq_vector) {
        perror("Error allocating memory.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate and Apply phasing vector
    float phase_increment = calc_phase_increment(beam_angle, (clear_freq_range[0] + clear_freq_range[1]) / 2, meta_data->x_spacing);
    if (VERBOSE) printf("phase_increment: %lf\n", phase_increment);
    // phase_increment = .000738
    
    for (int i = 0; i < meta_data->num_antennas; i++) {
        if (i <= IDX_LAST_MA) {
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
    if (VERBOSE) {
        printf("antenna[13]: %d\n", antennas[13]);
        printf("phasing_vector[13]: %f + %fi\n", creal(phasing_vector[13]), cimag(phasing_vector[13]));
    }

    // Apply beamforming
    for (int i = 0; i < num_samples; i++) {
        double real_sum = 0.0;
        double imag_sum = 0.0;
        
        for (int aidx = 0; aidx < meta_data->num_antennas; aidx++) {
            double real_sample = creal(raw_samples[aidx][i]);       
            double imag_sample = cimag(raw_samples[aidx][i]);
            double real_phase = creal(phasing_vector[aidx]);
            double imag_phase = cimag(phasing_vector[aidx]);
            if (VERBOSE && i == 2499) {
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
    if (VERBOSE) printf("beamformed[625]    = %f + %fi\n", creal(beamformed_samples[625]), cimag(beamformed_samples[625]));

    // Frequency Vector Calculation
    double delta_f = meta_data->usrp_rf_rate / num_samples;
    double f_start = meta_data->usrp_fcenter * 1000 - (meta_data->usrp_rf_rate / 2);
    for (int i = 0; i < num_samples; i++) {
        freq_vector[i] = i * delta_f + f_start;
    }

    /// Spectrum Calculation and Averaging; delinate Transmitters and filter out noise
    // Spectrum Avg (avg of 4 fft)
    // if (SPECTRAL_AVGING) {
    printf("=----Starting Spectral Average----=\n");
    clock_t t_avg_curr, t_avg;
    int avg_freq_ratio = 4;     //(int) delta_f / CLRFREQ_RES;
    int num_avg_samples = num_samples / avg_freq_ratio; 

    // Determine (Avg version of) Freq vector
    double *freq_vector_avg = (double*) malloc(sizeof(double) * num_avg_samples);
    double delta_f_avg = delta_f * avg_freq_ratio;
    delta_f = delta_f_avg;
    for (int i = 0; i < num_avg_samples; i++) freq_vector_avg[i] = i * delta_f_avg + f_start;

    double *avg_spectrum = (double*) fftw_malloc(sizeof(double) * num_avg_samples);
    fftw_complex *four_spectrums = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    if (VERBOSE) printf("num_avg_sample: %d\navg_freq_ratio: %d\n", num_avg_samples, avg_freq_ratio);

        
    // FFT Beamformed Samples
    t_avg_curr = clock();
    // for (int j = 0; j < num_samples; j++) {
    //     four_spectrums[j] = beamformed_samples[j];
    // }
    

    // TODO: Optimize plan usage; Long-term plan storage (create and store in samples_server.c)
    fftw_plan plan = fftw_plan_dft_1d(num_samples, beamformed_samples, four_spectrums, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);    

    // For every avg element,
    for (int k = 0; k < num_avg_samples; k++) {
        // Avg the magnitude of four sequential spectrum samples 
        for (int j = 0; j < avg_freq_ratio; j++) {
            double re = creal(four_spectrums[k * avg_freq_ratio + j]) * creal(four_spectrums[k * avg_freq_ratio + j]);
            double im = cimag(four_spectrums[k * avg_freq_ratio + j]) * cimag(four_spectrums[k * avg_freq_ratio + j]);
            if (j == 1) avg_spectrum[k] = sqrt(re + im);
            else avg_spectrum[k] += sqrt(re + im);

            if (k == 9) {
                printf("sample[%d][%d][%d]: %f + j%f\n", k, j, k + j, creal(four_spectrums[k + j]), cimag(four_spectrums[k + j]));
            }
        }
        avg_spectrum[k] /= avg_freq_ratio;

        if (k == 9) printf("avg_spectrum[%d]: %f\n", k, avg_spectrum[k]);
    }
    t_avg = clock() - t_avg_curr;
    if (VERBOSE) printf("====> Spectral Avg took (s): %lf\n", ((double) (t_avg)) / (CLOCKS_PER_SEC));

    if (VERBOSE) printf("avg_spectrum[0]: %f\n", avg_spectrum[0]);


    // Debug: Store 4 time-domain samples (in magnitude) and Not Averaged Freq Spectrum (complex to magnitude)
    // if (VERBOSE) {
    //     write_spectrum_csv("../Freq_Server/utils/csv_dump/spectrum_output.NoAvg.csv", four_spectrums, freq_vector, num_samples);
    //     double sample_set[num_avg_samples];
    //     int i_set = 0;  
    //     for (int i = 0; i < num_avg_samples; i++) {
    //         double re = creal(four_spectrums[i + i_set]) * creal(four_spectrums[i + i_set]);
    //         double im = cimag(four_spectrums[i + i_set]) * cimag(four_spectrums[i + i_set]);
    //         sample_set[i] = sqrt(re + im);
    //     }
    //     write_spectrum_mag_csv("../Freq_Server/utils/csv_dump/avg_spect.0.csv", sample_set, freq_vector_avg, num_avg_samples);
    //     i_set+= num_avg_samples;
    //     for (int i = 0; i < num_avg_samples; i++) {
    //         double re = creal(four_spectrums[i + i_set]) * creal(four_spectrums[i + i_set]);
    //         double im = cimag(four_spectrums[i + i_set]) * cimag(four_spectrums[i + i_set]);
    //         sample_set[i] = sqrt(re + im);
    //     }
    //     write_spectrum_mag_csv("../Freq_Server/utils/csv_dump/avg_spect.1.csv", sample_set, freq_vector_avg, num_avg_samples);
    //     i_set+= num_avg_samples;
    //     for (int i = 0; i < num_avg_samples; i++) {
    //         double re = creal(four_spectrums[i + i_set]) * creal(four_spectrums[i + i_set]);
    //         double im = cimag(four_spectrums[i + i_set]) * cimag(four_spectrums[i + i_set]);
    //         sample_set[i] = sqrt(re + im);
    //     }
    //     write_spectrum_mag_csv("../Freq_Server/utils/csv_dump/avg_spect.2.csv", sample_set, freq_vector_avg, num_avg_samples);
    //     i_set+= num_avg_samples;
    //     for (int i = 0; i < num_avg_samples; i++) {
    //         double re = creal(four_spectrums[i + i_set]) * creal(four_spectrums[i + i_set]);
    //         double im = cimag(four_spectrums[i + i_set]) * cimag(four_spectrums[i + i_set]);
    //         sample_set[i] = sqrt(re + im);
    //     }
    //     write_spectrum_mag_csv("../Freq_Server/utils/csv_dump/avg_spect.3.csv", sample_set, freq_vector_avg, num_avg_samples);
    // }
    

    // Dispose of temp variables
    fftw_destroy_plan(plan);
    fftw_free(four_spectrums);

    /// END of Spectrum Calc
    // TODO: Plot usrp/antenna im and re separately (4 plots) 
    //          to confirm that one of them is correctly all zeros for im

    if (VERBOSE) printf("delta_f: %f\nnum_samples: %d\nfcenter: %d\n", delta_f, num_samples, meta_data->usrp_fcenter * 1000);


    // Mask restricted frequencies
    

    int clear_sample_start = (int) round((clear_freq_range[0] - f_start) / delta_f);
    int clear_sample_end = (int) round((clear_freq_range[1] - f_start) / delta_f);
    for (int i = clear_sample_start; i < clear_sample_end; i++) {
        if ((i < 2 + clear_sample_start || i > clear_sample_end - 3) && VERBOSE) printf("spectrum_pow[%d]: %f\n", i, avg_spectrum[i]);   
    }

    // Find clear frequency
    double clear_bw = 4e6 / smsep; // ~ 300 us
    clear_bw = 0;
    clock_t t1, t2;
    t1 = clock();
    find_clear_freqs(avg_spectrum, *meta_data, delta_f, clear_freq_range[0], clear_freq_range[1], clear_bw, clr_bands);
    t2 = clock();
    if (VERBOSE) printf("clear_freq_search (ms): %lf\n", ((double) (t2 - t1)) / (CLOCKS_PER_SEC * 1000));

    // Debug: Output results
    if (VERBOSE) for (int i = 0; i < CLR_BANDS_MAX; i++)
        printf("Clear Freq Band[%d]: | %dMHz -- Noise: %f -- %dMHz |\n", i, clr_bands[i].f_start, clr_bands[i].noise, clr_bands[i].f_end);
    

    // Debug: Save data to csv
    write_sample_mag_csv(sample_im_file, sample_im, freq_vector, meta_data);          // Used to check complex Samples after Beamforming
    // write_sample_mag_csv(sample_re_file, sample_re, freq_vector, meta_data);          // Plot w/ sample_plot.py
    write_spectrum_mag_csv(spectrum_file, avg_spectrum, freq_vector, num_avg_samples);  // Spectrum after Spectrum FFT averaging; plot w/ spectrum_plot.py
    write_clr_freq_csv(clr_freq_file, clr_bands);                                       // Used to plot Clear Freq Bands w/ spectrum_plot.clr_freq.py

    printf("Finished Clear Freq Search!\n");
    
    fftw_free(phasing_vector);
    fftw_free(beamformed_samples);
    fftw_free(avg_spectrum);
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
    double smsep = .0003; // 1 / (2 * 250 * pow(10, 3));      // ~4 ms

    // Stopwatch Start
    double t1,t2;
    t1 = clock();

    // Find Clear Frequency Bands
    calc_clear_freq_on_raw_samples(
        raw_samples, &meta_data, restricted_frequencies, 
        clear_freq_range, beam_angle, smsep, clr_bands);
    
    // Print processing time; Stopwatch End
    t2 = clock();
    printf("clear_freq_search (ms): %lf\n", ((double) (t2 - t1)) / (CLOCKS_PER_SEC * 1000));
    
    // Free allocated memory
    free(meta_data.antenna_list);  


};