#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>      // FFT transform library
#include <signal.h>
#include "../../clear_freq_search.c"

#define SAMPLES_NUM 2500
#define ANTENNAS_NUM 14
// #define SHM_SIZE (2 * SAMPLES_NUM * ANTENNAS_NUM * sizeof(int))  // Size for a 2x14x2500 array of integers
#define SHM_SIZE (ANTENNAS_NUM * SAMPLES_NUM * sizeof(int)) // TODO: Finish 2x2500 test
#define SHM_NAME "/shared_memory"
#define SEM_SERVER "/sem_server"
#define SEM_CLIENT "/sem_client"


// Initiate with flags:
// -lrt -pthread -lfftw3 -lm



sem_t *sem_server = NULL;
sem_t *sem_client = NULL;
int *shm_ptr = NULL;
int shm_fd = -1;





// TODO: store sample_meta_data *meta_data, double **clear_freq_range
void storeData(int **temp_arr, fftw_complex ***raw_samples) {
    
    // Allocate mem
    *raw_samples = (fftw_complex **)fftw_malloc(ANTENNAS_NUM * sizeof(fftw_complex *));
    for (int i = 0; i < ANTENNAS_NUM; i++) {
        (*raw_samples)[i] = (fftw_complex *)fftw_malloc(SAMPLES_NUM * sizeof(fftw_complex));
    }
    if (*raw_samples == NULL) {
        perror("Error allocating memory for raw samples");
        exit(EXIT_FAILURE);
    }

    // Store data
    for (int i = 0; i < ANTENNAS_NUM; i++) {
        fftw_complex *ant_samples = (*raw_samples)[i];

        for (int j = 0; j < SAMPLES_NUM; j++) {
            double real, imag;
            ant_samples[j] = real + I * imag;
        }
    }
}

/**
 * @brief  Deallocates all service semaphores and SHM pointers.
 * @note   
 * @retval None
 */
void cleanup() {
    if (sem_server) sem_close(sem_server);
    if (sem_client) sem_close(sem_client);
    if (shm_ptr) munmap(shm_ptr, SHM_SIZE);
    if (shm_fd >= 0) close(shm_fd);
    sem_unlink(SHM_NAME);
    sem_unlink(SEM_SERVER);
    sem_unlink(SEM_CLIENT);
}

/**
 * @brief  Catch signals and exit gracefully.
 * @note   
 * @param  sig: Caught signal
 * @retval None
 */
void handle_sigint(int sig) {
    printf("\n[Frequency Server] Caught signal %d, cleaning up and exiting...\n", sig);
    cleanup();

    // Prompt exit to terminal  
    printf("[Frequency Server] Main processes and communication terminated.\n"
           "Goodbye.\n");
           
    exit(0);
}



int main() {
    // Setup Signal Handler
    signal(SIGINT, handle_sigint);

    // Open Shared Memory Object
    printf("[Frequency Server] Initializing Shared Memory Object...\n");
    int shm_fd = shm_open(SHM_NAME, O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1) {
        perror("[Frequency Server] shm_open failed.\n");
        exit(EXIT_FAILURE);
    } else {
        printf("[Frequency Server] Created Shared Memory Object...\n");
    }

    // Set Size of Shared Memory Object
    if (ftruncate(shm_fd, SHM_SIZE) == -1) {
        perror("[Frequency Server] ftruncate failed\n");
        exit(EXIT_FAILURE);
    }
    
    // Request Block of Memory
    printf("[Frequency Server] Requesting Shared Memory Cache...\n");    
    int *shm_ptr = mmap(0, SHM_SIZE, PROT_WRITE | PROT_READ, MAP_SHARED, shm_fd, 0);
    if (shm_ptr == MAP_FAILED) {
        perror("[Frequency Server] Memory Mapping failed.\n");
        exit(EXIT_FAILURE);
    } else {
        printf("[Frequency Server] Memory successfully cached...\n");
    }

    // Open Semaphores for synchronization
    printf("[Frequency Server] Opening Communication Semaphores...\n");
    sem_t *sem_server = sem_open(SEM_SERVER, O_CREAT, 0666, 0);
    sem_t *sem_client = sem_open(SEM_CLIENT, O_CREAT, 0666, 0);
    if (sem_server == SEM_FAILED || sem_client == SEM_FAILED) {
        perror("[Frequency Server] One or more sem_open failed.\n");
        exit(EXIT_FAILURE);
    } else {
        printf("[Frequency Server] All semaphores have been established...\n");
    }

    printf("[Frequency Server] Done Initializing...\n\n");

    // Allocate mem
    fftw_complex **temp_samples = NULL;
    temp_samples = (fftw_complex **)fftw_malloc(ANTENNAS_NUM * sizeof(fftw_complex *));
    if (temp_samples == NULL) {
        perror("Error allocating memory for temp_samples pointers");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < ANTENNAS_NUM; i++) {
        temp_samples[i] = (fftw_complex *)fftw_malloc(SAMPLES_NUM * sizeof(fftw_complex));
        if (temp_samples[i] == NULL) {
            perror("Error allocating memory for temp_samples elements");
            exit(EXIT_FAILURE);
        }
    }
    freq_band *clr_bands = NULL;
    clr_bands = (freq_band *)malloc(CLR_BANDS_MAX * sizeof(freq_band));
    if (clr_bands == NULL) {
        perror("Error allocating memory for clr_bands elements");
        exit(EXIT_FAILURE);
    }
            
    // Continuously Send and Receive messages via Shared Memory
    while (1) {
        // Write data to Shared Memory
        printf("[Frequency Server] Writing frequency data to Shared Memory...\n");
        for (int i = 0; i < ANTENNAS_NUM; i++) {
            for (int i = 0; i < SAMPLES_NUM; i++) {
                shm_ptr[i] = i; 
            }
        }
        printf("[Frequency Server] Requesting new client to respond...\n\n");
        sem_post(sem_client); 
        
        printf("[Frequency Server] Awaiting client response...\n");
        sem_wait(sem_server);   

        printf("[Frequency Server] Processing client frequency data...\n");
        printf("[Frequency Server; from Client] ");
        // for (int i = 0; i < 20; i += 2) { 
        //     printf("    %d + j%d\n", shm_ptr[i], shm_ptr[i+1]);
        // }
        // Store data into complex form
        for (int i = 0; i < ANTENNAS_NUM; i++) {    
            for (int j = 0; j < SAMPLES_NUM; j+=2) {
                temp_samples[i][j] = shm_ptr[i * SAMPLES_NUM + j] + I * shm_ptr[i * SAMPLES_NUM + j + 1];

                // Print first complex of each antenna sample batch for verification
                // if (j < 20) {
                //     printf("shm[%d]      =   %d + i%d\n", i * SAMPLES_NUM + j, (int)shm_ptr[i * SAMPLES_NUM + j], (int)shm_ptr[i * SAMPLES_NUM + j + 1]);   
                //     printf("vs\n");
                //     printf("temp_samples[%d][%d] =  %f + i%f\n\n", i, j, creal(temp_samples[i][j]), cimag(temp_samples[i][j]));
                // }
            }
        }
        printf("\n");

        // Process data for clear freqs and store result
        clear_freq_search(temp_samples, clr_bands);
        printf("[Frequency Server] Processed Client response successfully...\n");
    }

    return 0;
}