#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>      // FFT transform library

typedef struct {
    int band_start;
    int band_end;
    double noise;
} clr_band;

void calculate_noise_bands(fftw_complex *spectrum, int num_samples, int band_size, clr_band *noise_bands, int *num_bands) {
    int num_total_bands = num_samples / band_size;
    *num_bands = num_total_bands;

    for (int i = 0; i < num_total_bands; i++) {
        noise_bands[i].band_start = i * band_size;
        noise_bands[i].band_end = (i + 1) * band_size - 1;
        noise_bands[i].noise = 0.0;

        for (int j = noise_bands[i].band_start; j <= noise_bands[i].band_end; j++) {
            noise_bands[i].noise += cabs(spectrum[j]);
        }

        noise_bands[i].noise /= band_size; // Average noise
    }
}


/**
 * @brief  Compares noise between two different bands
 * @note   
 * @param  *a: 
 * @param  *b: 
 * @retval 
 */
int compare_noise(const void *a, const void *b) {
    clr_band *bandA = (clr_band *)a;
    clr_band *bandB = (clr_band *)b;
    return (bandA->noise > bandB->noise) - (bandA->noise < bandB->noise);
}

void find_lowest_noise_bands(clr_band *noise_bands, int num_bands, clr_band *lowest_noise_bands, int num_lowest) {
    qsort(noise_bands, num_bands, sizeof(clr_band), compare_noise);

    for (int i = 0; i < num_lowest; i++) {
        lowest_noise_bands[i] = noise_bands[i];
    }
}

