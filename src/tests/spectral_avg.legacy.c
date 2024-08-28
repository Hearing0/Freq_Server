// Legacy test code for determining the average computation time for the fft spectral avgeraging protocal
// 
// To use simply copy all in this main and replace the spectral averaging section in clear_freq_search.c.
// Run the server-client model and the server should output the average time for this section.

int main(int argc, char const *argv[])
{
    /// Spectrum Averaging; delinate Transmitters and filter out noise
    // Spectrum Avg (4 fft into avg)
    clock_t t_avg_curr, t_avg;
    int avg_freq_ratio = 4;     //(int) delta_f / CLRFREQ_RES;
    int num_avg_samples = num_samples / avg_freq_ratio; 


    /// Computes one avg in ~1.100 s
    printf("=----------------=\n");
    printf("starting avg...\n");

    // Avg version of freq vector
    double *freq_vector_avg = (double*) malloc(sizeof(double) * num_avg_samples);
    double delta_f_avg = delta_f * avg_freq_ratio;
    for (int i = 0; i < num_avg_samples; i++) freq_vector_avg[i] = i * delta_f_avg + f_start;
    
    double *avg_spect_samples = (double*) fftw_malloc(sizeof(double) * num_avg_samples);
    fftw_complex *samples = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    printf("num_avg_sample: %d\navg_freq_ratio: %d\n", num_avg_samples, avg_freq_ratio);

    


    // TODO: Optimize plan usage; Long-term plan storage (create and store in samples_server.c)
    fftw_plan plan = fftw_plan_dft_1d(num_samples, samples, samples, FFTW_FORWARD, FFTW_ESTIMATE);

    // Find avg compute time of 100 fft avgs 
    t_avg_curr = clock();
    // for (int i = 0; i < 100; i++) {
        
        // Copy beamformed_samples for FFT
        for (int j = 0; j < num_samples; j++) samples[j] = beamformed_samples[j];
        
        fftw_execute(plan);
        // fft_samples(samples, num_samples, samples);

        // For each avg_sample element,
        for (int k = 0; k < num_avg_samples; k++) {
            // Avg the corresponding elements from the split chunks 
            for (int j = 0; j < num_samples; j += num_avg_samples) {
                double re = creal(samples[k + j * avg_freq_ratio]) * creal(samples[k + j * avg_freq_ratio]);
                double im = cimag(samples[k + j * avg_freq_ratio]) * cimag(samples[k + j * avg_freq_ratio]);
                if (j == 1) avg_spect_samples[k] = sqrt(re + im);
                else avg_spect_samples[k] += sqrt(re + im);
            }
            avg_spect_samples[k] /= avg_freq_ratio;

            if (k == 4) printf("avg_spect_samples[%d]: %f\n", k, avg_spect_samples[k]);
        }

        // Store spectrum to verify plotting
        // if (i == 0)    
    // }
    t_avg = clock() - t_avg_curr;
    // t_avg /= 100;
    printf("====> avg fft v1 (s): %lf\n", ((double) (t_avg)) / (CLOCKS_PER_SEC));

    write_spect_mag_to_csv("../Freq_Server/utils/csv_dump/spectrum_output.avg1.csv", avg_spect_samples, freq_vector_avg, num_avg_samples);
    printf("avg_spect_samples[0]: %f\n", avg_spect_samples[0]);

    fftw_destroy_plan(plan);
    fftw_free(samples);
    fftw_free(avg_spect_samples);
    return 0;
}