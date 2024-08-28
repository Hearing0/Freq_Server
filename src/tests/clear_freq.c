#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define CLRFREQ_RES 1.0
#define CLR_BANDS_MAX 10

typedef struct {
    double f_start;
    double f_end;
    double noise;
} freq_band;

typedef struct {
    double usrp_rf_rate;
    double usrp_fcenter;
} sample_meta_data;

// Function to perform convolution
int* convolve(int* u, int u_size, int* v, int v_size) {
    int result_size = u_size + v_size - 1;
    int* result = (int*)calloc(result_size, sizeof(int));
    for (int i = 0; i < result_size; i++) {
        for (int j = 0; j < v_size; j++) {
            if (i - j >= 0 && i - j < u_size) {
                result[i] += u[i - j] * v[j];
            }
        }
    }
    return result;
}

// Function to find clear frequency bands
void find_clear_freqs(fftw_complex *spectrum_power, sample_meta_data meta_data, double delta_f, double f_start, double f_end, int clear_bw, freq_band *clr_freq_bands) {
    printf("[find_clear_freqs()] Entered find_clear_freqs()...\n");
    if (clear_bw == 0) clear_bw = 40e3;
    int clear_sample_bw = clear_bw / CLRFREQ_RES;  

    // Define Clear Freq Search Range
    int start_sample = (int) (f_start / delta_f);
    int end_sample = (int) (f_end / delta_f);

    // Trim Spectrum Data to Search Range (Used for convolving)
    int search_bw = end_sample - start_sample;
    int* search_band = (int*)malloc(sizeof(int) * search_bw);
    for (int i = 0; i < search_bw; i++) {
        search_band[i] = (int) round(sqrt(creal(spectrum_power[start_sample + i]) * creal(spectrum_power[start_sample + i]) +
                                           cimag(spectrum_power[start_sample + i]) * cimag(spectrum_power[start_sample + i])));
        if (i % 10 == 0) printf("search_band[%d]: %d\n", i, search_band[i]);
    }
    
    // Preconfig Clear Freq Bands
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        clr_freq_bands[i].f_start = f_start;
        clr_freq_bands[i].f_end = f_end;
        clr_freq_bands[i].noise = 500;
    }

    printf("[find_clear_freqs()] Scanning Search Range w/ Bandpass...\n");

    int* bpf_result;
    int* bpf = (int*)malloc(sizeof(int) * clear_sample_bw);
    for (int band_i = 0; band_i < clear_sample_bw; band_i++) {
        bpf[band_i] = 1;
    }

    for (int i = 0; i < search_bw - clear_sample_bw; i++) {
        printf("[find_clear_freqs()] BPF @ f_start = %d...\n", i);

        // Convolve BPF with Search Range
        bpf_result = convolve(&search_band[i], clear_sample_bw, bpf, clear_sample_bw);
        printf("    Convolved Search Band and BPF...\n");

        // Calculate noise
        int noise_sum = 0;
        for (int j = 0; j < clear_sample_bw; j++) {
            noise_sum += bpf_result[j];
        }
        double avg_noise = (double)noise_sum / clear_sample_bw;

        // Update clear frequency bands if necessary
        if (avg_noise < clr_freq_bands[CLR_BANDS_MAX - 1].noise) {
            clr_freq_bands[CLR_BANDS_MAX - 1].f_start = (start_sample + i) * delta_f;
            clr_freq_bands[CLR_BANDS_MAX - 1].f_end = (start_sample + i + clear_sample_bw) * delta_f;
            clr_freq_bands[CLR_BANDS_MAX - 1].noise = avg_noise;

            // Sort bands by noise
            for (int j = CLR_BANDS_MAX - 1; j > 0; j--) {
                if (clr_freq_bands[j].noise < clr_freq_bands[j - 1].noise) {
                    freq_band temp = clr_freq_bands[j];
                    clr_freq_bands[j] = clr_freq_bands[j - 1];
                    clr_freq_bands[j - 1] = temp;
                }
            }
        }

        free(bpf_result);
    }

    free(bpf);
    free(search_band);
    printf("[find_clear_freqs()] Exiting find_clear_freqs()...\n");
}

int main() {
    // Initialize sample metadata
    sample_meta_data meta_data;
    meta_data.usrp_rf_rate = 1e6; // 1 MHz
    meta_data.usrp_fcenter = 2.4e9; // 2.4 GHz

    // Frequency step per sample
    double delta_f = 1e3; // 1 kHz

    // Clear frequency search bounds
    double f_start = 0.0;
    double f_end = 1e6; // 1 MHz

    // Bandwidth of the clear frequency bands
    int clear_bw = 40e3; // 40 kHz

    // Allocate memory for spectrum power
    fftw_complex* spectrum_power = (fftw_complex*)malloc(sizeof(fftw_complex) * SAMPLE_COUNT);
    for (int i = 0; i < SAMPLE_COUNT; i++) {
        spectrum_power[i] = (cos(2 * PI * i / SAMPLE_COUNT) + I * sin(2 * PI * i / SAMPLE_COUNT));
    }

    // Allocate memory for clear frequency bands
    freq_band* clr_freq_bands = (freq_band*)malloc(sizeof(freq_band) * CLR_BANDS_MAX);

    // Call the function to find clear frequency bands
    find_clear_freqs(spectrum_power, meta_data, delta_f, f_start, f_end, clear_bw, clr_freq_bands);

    // Print the results
    printf("Clear Frequency Bands:\n");
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        printf("Band %d: f_start = %f Hz, f_end = %f Hz, noise = %f\n", i, clr_freq_bands[i].f_start, clr_freq_bands[i].f_end, clr_freq_bands[i].noise);
    }

    // Free allocated memory
    free(spectrum_power);
    free(clr_freq_bands);

    return 0;
}
