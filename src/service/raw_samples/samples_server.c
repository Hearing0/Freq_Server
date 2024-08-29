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
#define RESTRICT_NUM 15             // Number of restricted freq bands in the restrict.dat.inst
#ifndef CLR_BANDS_MAX
#define CLR_BANDS_MAX 6
#endif
#define RESTRICT_SHM_SIZE (RESTRICT_NUM * 2 * sizeof(int))          // 2 = start and end freqs
#define SAMPLES_SHM_SIZE (ANTENNAS_NUM * SAMPLES_NUM * sizeof(int)) // TODO: Finish 2x2500 test
#define CLR_BANDS_SHM_SIZE (CLR_BANDS_MAX * sizeof(double))         // TODO: Round to convert freqs to int again 
// #define SAMPLES_SHM_SIZE (2 * SAMPLES_NUM * ANTENNAS_NUM * sizeof(int))  // Size for a 2x14x2500 array of integers
#define SAMPLES_SHM_NAME "/samples"
#define RESTRICT_SHM_NAME "/restricted_freq"
#define SEM_SERVER "/sem_server"
#define SEM_CLIENT "/sem_client"




// Initiate with flags:
// -lrt -pthread -lfftw3 -lm

// TODO: Read/Write of specific data types (only write samples and restrict vs read clr freq)
// TODO: Test Clr Freq by flipping logic at lab
// TODO: Calc all beam directions at recieve of samples
// TODO: Num of Radar {has table of all beam dirc {has Clr freq in beam direction}}} while sharing restricting assigned freq 
// TODO: Rewrite of usrp sample send/clr freq request timing logic 

sem_t *sem_server = NULL;
sem_t *sem_client = NULL;
int *restrict_shm_ptr = NULL;
int restrict_shm_fd = -1;
int *samples_shm_ptr = NULL;
int samples_shm_fd = -1;




/**
 * @brief  Deallocates all service semaphores and SHM pointers.
 * @note   
 * @retval None
 */
void cleanup() {
    if (sem_server) sem_close(sem_server);
    if (sem_client) sem_close(sem_client);
    if (samples_shm_ptr) munmap(samples_shm_ptr, SAMPLES_SHM_SIZE);
    if (samples_shm_fd >= 0) close(samples_shm_fd);
    if (restrict_shm_ptr) munmap(restrict_shm_ptr, RESTRICT_SHM_SIZE);
    if (restrict_shm_fd >= 0) close(restrict_shm_fd);
    sem_unlink(SAMPLES_SHM_NAME);
    sem_unlink(RESTRICT_SHM_NAME);
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
    samples_shm_fd = shm_open(SAMPLES_SHM_NAME, O_CREAT | O_RDWR, 0666);
    restrict_shm_fd = shm_open(RESTRICT_SHM_NAME, O_CREAT | O_RDWR, 0666);
    
    if (samples_shm_fd == -1) {
        perror("[Frequency Server] shm_open failed for samples_shm");
        exit(EXIT_FAILURE);
    } if (restrict_shm_fd == -1) {
        perror("[Frequency Server] shm_open failed for restrict_shm");
        exit(EXIT_FAILURE);
    } else {
        printf("[Frequency Server] Created Shared Memory Objects...\n");
    }

    // Set Size of Shared Memory Object
    if (ftruncate(samples_shm_fd, SAMPLES_SHM_SIZE) == -1 || ftruncate(restrict_shm_fd, RESTRICT_SHM_SIZE) == -1) {
        perror("[Frequency Server] ftruncate failed\n");
        exit(EXIT_FAILURE);
    }
    
    // Request Block of Memory
    printf("[Frequency Server] Requesting Shared Memory Cache...\n");    
    int *samples_shm_ptr = mmap(0, SAMPLES_SHM_SIZE, PROT_WRITE | PROT_READ, MAP_SHARED, samples_shm_fd, 0);
    int *restrict_shm_ptr = mmap(0, RESTRICT_SHM_SIZE, PROT_WRITE | PROT_READ, MAP_SHARED, restrict_shm_fd, 0);
    if (samples_shm_ptr == MAP_FAILED || restrict_shm_ptr == MAP_FAILED) {
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

    // Allocate temp mem for shm varibles
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
    // TODO: could restricted freq array size be dynamic?
    freq_band *restricted_freq = NULL;
    restricted_freq = (freq_band *)malloc(RESTRICT_SHM_SIZE);
    if (restricted_freq == NULL) {
        perror("Error allocating memory for restricted_freq elements");
        exit(EXIT_FAILURE);
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
                samples_shm_ptr[i] = i; 
            }
        }
        printf("[Frequency Server] Requesting new client to respond...\n\n");
        sem_post(sem_client); 
        
        // TODO: Client writes to param shm
        printf("[Frequency Server] Awaiting client response...\n");
        sem_wait(sem_server);   

        // TODO: Add Param Buffer for meta data and clear freq range
        printf("[Frequency Server] Processing client frequency data...\n");
        printf("[Frequency Server; from Client] \n");
        // Debug: View 
        // for (int i = 0; i < 20; i += 2) printf("    %d + j%d\n", samples_shm_ptr[i], samples_shm_ptr[i+1]);
        // for (int i = 0; i < RESTRICT_NUM; i++) printf("    restrict[%d]: %d\n", i, restrict_shm_ptr[i]);
        // Store sample data into complex form
        for (int i = 0; i < ANTENNAS_NUM; i++) {    
            for (int j = 0; j < SAMPLES_NUM; j+=2) {
                temp_samples[i][j] = samples_shm_ptr[i * SAMPLES_NUM + j] + I * samples_shm_ptr[i * SAMPLES_NUM + j + 1];

                // Print first complex of each antenna sample batch for verification
                // if (j < 20) {
                //     printf("shm[%d]      =   %d + i%d\n", i * SAMPLES_NUM + j, (int)samples_shm_ptr[i * SAMPLES_NUM + j], (int)samples_shm_ptr[i * SAMPLES_NUM + j + 1]);   
                //     printf("vs\n");
                //     printf("temp_samples[%d][%d] =  %f + i%f\n\n", i, j, creal(temp_samples[i][j]), cimag(temp_samples[i][j]));
                // }
            }
        }
        /// Transfer param data into referencable forms
        // Transfer restricted frequencies
        for (int i = 0; i < RESTRICT_NUM; i++) {                
            // Store Restricted Freq
            // TODO: Verify this storing logic
            restricted_freq[i].f_start = restrict_shm_ptr[i * 2];
            restricted_freq[i].f_end = restrict_shm_ptr[i * 2 + 1];
            printf("restricted_freq[%d]: |%d -- %d|\n", i, restricted_freq[i].f_start, restricted_freq[i].f_end);
        }
        // TODO: Add the rest of meta data storing (clear_freq, etc.)
        printf("\n");

        // Process data for clear freqs and store result
        // clear_freq_search(temp_samples, clr_bands, restricted_freq);
        printf("[Frequency Server] Processed Client response successfully...\n");
    }

    return 0;
}