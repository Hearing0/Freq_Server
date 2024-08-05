#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define SPECTRUM_SIZE 2500
#define FILTER_SIZE 100
#define NUM_BANDS 6
#define RAND_MAX 500

// Function to perform convolution
void convolve(const double* spectrum, const double* filter, double* result, int spectrum_size, int filter_size) {
    int i, j;
    for (i = 0; i < spectrum_size; i++) {
        result[i] = 0;
        for (j = 0; j < filter_size; j++) {
            if (i - j >= 0) {
                result[i] += spectrum[i - j] * filter[j];
            }
        }
    }
}

// Function to find the smallest non-intersecting bands
void find_smallest_bands(const double* convolved, int* bands, int spectrum_size, int num_bands) {
    int i, j, band_start, band_end;
    double band_sum;
    int smallest_bands[NUM_BANDS][2]; // Store start and end indices of smallest bands

    // Setup Smallest Bands
    for (i = 0; i < num_bands; i++) {
        smallest_bands[i][0] = RAND_MAX;
        smallest_bands[i][1] = -1;
    }

    for (i = 0; i < spectrum_size; i++) {
        for (j = i + 1; j < spectrum_size; j++) {
            band_start = i;
            band_end = j;
            band_sum = 0;
            for (int k = band_start; k <= band_end; k++) {
                band_sum += convolved[k];
            }

            // Check if the band is smaller than the largest band in smallest_bands
            double max_band_sum = -1;
            int max_band_index = -1;
            for (int k = 0; k < num_bands; k++) {
                if (smallest_bands[k][0] != -1 && smallest_bands[k][1] != -1) {
                    double current_band_sum = 0;
                    for (int l = smallest_bands[k][0]; l <= smallest_bands[k][1]; l++) {
                        current_band_sum += convolved[l];
                    }
                    if (current_band_sum > max_band_sum) {
                        max_band_sum = current_band_sum;
                        max_band_index = k;
                    }
                } else {
                    max_band_index = k;
                    break;
                }
            }

            if (band_sum < max_band_sum || max_band_index == -1) {
                // Check if the band intersects with any existing smallest bands
                int intersects = 0;
                for (int k = 0; k < num_bands; k++) {
                    if (smallest_bands[k][0] != -1 && smallest_bands[k][1] != -1) {
                        if ((band_start >= smallest_bands[k][0] && band_start <= smallest_bands[k][1]) ||
                            (band_end >= smallest_bands[k][0] && band_end <= smallest_bands[k][1])) {
                            intersects = 1;
                            break;
                        }
                    }
                }
                if (!intersects) {
                    smallest_bands[max_band_index][0] = band_start;
                    smallest_bands[max_band_index][1] = band_end;
                }
            }
        }
    }

    for (i = 0; i < num_bands; i++) {
        bands[i * 2] = smallest_bands[i][0];
        bands[i * 2 + 1] = smallest_bands[i][1];
    }
}

int main() {
    int spectrum[SPECTRUM_SIZE];
    int filter[FILTER_SIZE];
    int convolved[SPECTRUM_SIZE];
    int bands[NUM_BANDS * 2]; // Store start and end indices of the bands

    // Initialize spectrum and filter with some values
    for (int i = 0; i < SPECTRUM_SIZE; i++) {
        spectrum[i] = (int)rand() / RAND_MAX;
    }
    for (int i = 0; i < FILTER_SIZE; i++) {
        filter[i] = RAND_MAX;
    }

    // Perform convolution
    convolve(spectrum, filter, convolved, SPECTRUM_SIZE, FILTER_SIZE);

    // Find smallest non-intersecting bands
    // find_smallest_bands(convolved, bands, SPECTRUM_SIZE, NUM_BANDS);

    // Print the smallest bands
    printf("Smallest bands:\n");
    for (int i = 0; i < NUM_BANDS; i++) {
        printf("Band %d: Start = %d, End = %d\n", i + 1, bands[i * 2], bands[i * 2 + 1]);
    }

    return 0;
}
