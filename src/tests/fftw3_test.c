#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>

// Function to calculate phase increment
double calc_phase_increment(double beam_angle, double frequency, double x_spacing) {
    
    return 0; 
}

// Function to convert radians to rectangular coordinates
double complex rad_to_rect(double phase) {
    return cexp(I * phase);
}

// Function to beamform samples
void beamform_samples(fftw_complex *raw_samples, fftw_complex *phasing_vector, int num_samples, fftw_complex *beamformed_samples) {
    for (int i = 0; i < num_samples; i++) {
        beamformed_samples[i] = raw_samples[i] * phasing_vector[i];
    }
}

// Function to calculate the FFT
void fft_clrfreq_samples(fftw_complex *beamformed_samples, int num_samples, fftw_complex *spectrum_power) {
    fftw_plan plan = fftw_plan_dft_1d(num_samples, beamformed_samples, spectrum_power, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

// Main function
int main() {
    int num_samples = 1024; // Example number of samples
    fftw_complex *raw_samples = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_complex *phasing_vector = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_complex *beamformed_samples = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_complex *spectrum_power = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);

    // Example initialization (replace with actual data)
    for (int i = 0; i < num_samples; i++) {
        raw_samples[i] = i + I * i;
        phasing_vector[i] = rad_to_rect(i * 0.1);
    }

    // Beamform samples
    beamform_samples(raw_samples, phasing_vector, num_samples, beamformed_samples);

    // Calculate FFT
    fft_clrfreq_samples(beamformed_samples, num_samples, spectrum_power);

    // Print the spectrum power
    for (int i = 0; i < num_samples; i++) {
        printf("Frequency bin %d: %f + %fi\n", i, creal(spectrum_power[i]), cimag(spectrum_power[i]));
    }

    // Free allocated memory
    fftw_free(raw_samples);
    fftw_free(phasing_vector);
    fftw_free(beamformed_samples);
    fftw_free(spectrum_power);

    return 0;
}
