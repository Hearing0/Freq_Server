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
#define ANTENNAS_NUM 2 //14
#define RESTRICT_NUM 15             // Number of restricted freq bands in the restrict.dat.inst
#ifndef CLR_BANDS_MAX
#define CLR_BANDS_MAX 6
#endif

#define RESTRICT_SHM_SIZE (RESTRICT_NUM * 2 * sizeof(int))          // 2 = start and end freqs
#define SAMPLES_SHM_SIZE (ANTENNAS_NUM * SAMPLES_NUM * sizeof(int)) // TODO: Finish 2x2500 test
#define CLR_BANDS_SHM_SIZE (CLR_BANDS_MAX * 4 * 3)         // TODO: Round to convert freqs to int again 
// #define SAMPLES_SHM_SIZE (2 * SAMPLES_NUM * ANTENNAS_NUM * sizeof(int))  // Size for a 2x14x2500 array of integers

// Shared Memory and Semaphore Names 
#define SAMPLES_SHM_NAME    "/samples"
#define RESTRICT_SHM_NAME   "/restricted_freq"
#define CLRFREQ_SHM_NAME    "/clear_freq"
#define PARAM_NUM 3

#define SEM_CLIENT  "/sem_client"               // For Sync and reserving client and server roles during data transfer
#define SEM_SERVER  "/sem_server"               
#define SEM_DATA    "/sem_data"                 // For multiple data transfers on single instance 
#define SEM_SAMPLES "/sem_samples"              // For Data locking b/w write/reads
#define SEM_INIT    "/sem_init"                 // init = initialization
#define SEM_CLRFREQ "/sem_clrfreq"
#define SEM_NUM 6

typedef struct shm_obj{
    const char* name;
    int* shm_ptr;
    int shm_fd;
    size_t size;
} shm_obj;

typedef struct semaphore {
    const char* name;
    sem_t* sem;
} semaphore;



// Build with the following flags:
// -lrt -pthread -lfftw3 -lm

// TODO: Read/Write of specific data types (only write samples and restrict vs read clr freq)
// TODO: Calc all beam directions at recieve of samples
// TODO: Num of Radar {has table of all beam dirc {has Clr freq in beam direction}}} while sharing restricting assigned freq 
// TODO: Rewrite of usrp sample send/clr freq request timing logic 

struct semaphore s_client   = {SEM_CLIENT,  NULL};
struct semaphore s_server   = {SEM_SERVER,  NULL};
struct semaphore s_data     = {SEM_DATA,    NULL};
struct semaphore s_samples  = {SEM_SAMPLES, NULL};
struct semaphore s_init     = {SEM_INIT,    NULL};
struct semaphore s_clrfreq  = {SEM_CLRFREQ, NULL};
struct semaphore *semaphores[SEM_NUM] = {
    &s_client,
    &s_server,    
    &s_data,
    &s_samples,
    &s_init,
    &s_clrfreq,
};

shm_obj samples_obj = {SAMPLES_SHM_NAME, NULL, -1, SAMPLES_SHM_SIZE};
shm_obj restrict_obj = {RESTRICT_SHM_NAME, NULL, -1, RESTRICT_SHM_SIZE};
shm_obj clrfreq_obj = {CLRFREQ_SHM_NAME, NULL, -1, CLR_BANDS_SHM_SIZE};
struct shm_obj *objects[PARAM_NUM] = {
    &samples_obj, 
    &restrict_obj, 
    &clrfreq_obj,
};





// TODO: Read shm_obj ptr into temp_data
// Note: obj and data must be the same data type
void read_restrict_shm(freq_band *restricted_freq, int *restrict_shm_ptr) {
    for (int i = 0; i < RESTRICT_NUM; i++)
    {
        // Store Restricted Freq
        restricted_freq[i].f_start = restrict_shm_ptr[i * 2];
        restricted_freq[i].f_end = restrict_shm_ptr[i * 2 + 1];
        // if (i < 3) printf("restricted_freq[%d]: |%d -- %d| ...\n", i, restricted_freq[i].f_start, restricted_freq[i].f_end);
    }
}

void read_sample_shm(fftw_complex **temp_samples, int *samples_shm_ptr) {
    // Store sample data into complex form
    for (int i = 0; i < ANTENNAS_NUM; i++)
    {
        for (int j = 0; j < SAMPLES_NUM; j += 2)
        {
            temp_samples[i][j] = samples_shm_ptr[i * SAMPLES_NUM + j] + I * samples_shm_ptr[i * SAMPLES_NUM + j + 1];

            // Debug: Print 20 complex of each antenna batch
            // if (j < 20) {
            //     printf("shm[%d]      =   %d + i%d\n", i * SAMPLES_NUM + j, (int)samples_shm_ptr[i * SAMPLES_NUM + j], (int)samples_shm_ptr[i * SAMPLES_NUM + j + 1]);
            //     printf("vs\n");
            //     printf("temp_samples[%d][%d] =  %f + i%f\n\n", i, j, creal(temp_samples[i][j]), cimag(temp_samples[i][j]));
            // }
        }
    }
}

/**
 * @brief  Reads in Clear Frequency Bands from its shared memory pointer.
 * @note   Used for debugging
 * @param  *clr_bands: Clear Frequency Bands
 * @param  *ptr: Shared Memory Pointer for Clear Frequency Bands
 * @retval None
 */
void read_clrfreq_shm(freq_band *clr_bands, int *ptr) {
    int elements_per_band = 3;
    // Store into clr_bands
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        clr_bands[i].f_start= ptr[i * elements_per_band];
        clr_bands[i].f_end  = ptr[i * elements_per_band + 1];
        clr_bands[i].noise  = ptr[i * elements_per_band + 2];
    }
}

/**
 * @brief  Writes Clear Freq Bands to its shared memory pointer.
 * @note   
 * @param  *clr_bands: Clear Frequency Bands
 * @param  *ptr: Shared Memory Pointer for Clear Frequency Bands
 * @retval None
 */
void write_clrfreq_shm(freq_band *clr_bands, int *ptr) {
    int elements_per_band = 3;
    for (int i = 0; i < CLR_BANDS_MAX; i++) {
        ptr[i * elements_per_band]      = clr_bands[i].f_start;
        ptr[i * elements_per_band + 1]  = clr_bands[i].noise;
        ptr[i * elements_per_band + 2]  = clr_bands[i].f_end;
    }
}

/**
 * @brief  Unlinks/deallocates Shared Memory Object mapping (ptr), file 
 * * descriptor (fd), and semaphore name (name).
 * @note   
 * @param  shm_obj: Shared Memory Object struct.
 * @retval None
 */
void clean(shm_obj shm_obj) {
    if (shm_obj.shm_ptr) munmap(shm_obj.shm_ptr, shm_obj.size);
    if (shm_obj.shm_fd >= 0) close(shm_obj.shm_fd);
    sem_unlink(shm_obj.name);
}

void clean_sem(semaphore sem) {
    if ((sem.sem)) sem_close( (sem.sem));
}

/**
 * @brief  Deallocates all service semaphores and SHM pointers.
 * @note   
 * @retval None
 */
void cleanup() {
    for (int i = 0; i < SEM_NUM; i++) clean_sem(*semaphores[i]);
    // if (sem_server) sem_close(sem_server);
    // if (sem_client) sem_close(sem_client);
    for (int i = 0; i < PARAM_NUM; i++) clean(*objects[i]);

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
    for (int i = 0; i < PARAM_NUM; i++){
        objects[i]->shm_fd = shm_open(objects[i]->name, O_CREAT | O_RDWR, 0666);
        if (objects[i]->shm_fd == -1) {
            printf("[Frequency Server] shm_open failed for %s\n", objects[i]->name);
            exit(EXIT_FAILURE);
        }
    }
    printf("[Frequency Server] Created Shared Memory Objects...\n");
    
    // Set Size of Shared Memory Object
    for (int i = 0; i < PARAM_NUM; i++){
        if (ftruncate(objects[i]->shm_fd, objects[i]->size) == -1) {
            perror("[Frequency Server] ftruncate failed\n");
            exit(EXIT_FAILURE);
        }
    }
    // Request Block of Memory
    printf("[Frequency Server] Requesting Shared Memory Cache...\n");    
    for (int i = 0; i < PARAM_NUM; i++) {
        objects[i]->shm_ptr = mmap(0, objects[i]->size, PROT_WRITE | PROT_READ, MAP_SHARED, objects[i]->shm_fd, 0);
        if (objects[i]->shm_ptr == MAP_FAILED) {
            printf("[Frequency Server] Memory Mapping failed for %s\n", objects[i]->name);
            exit(EXIT_FAILURE);
        }
    }
    printf("[Frequency Server] Memory successfully cached...\n");


    // Open Semaphores for synchronization     
    printf("[Frequency Server] Opening Communication Semaphores...\n");    
    for (int i = 0; i < SEM_NUM; i++) {
        semaphores[i]->sem = sem_open(semaphores[i]->name, O_CREAT, 0666, 0);
        if (semaphores[i]->sem == SEM_FAILED) {
            printf("[Frequency Server] \"%s\" sem_open failed.\n", semaphores[i]->name);
            exit(EXIT_FAILURE);    
        } 
        // else {
        //     printf("[Frequency Server] \"%s\" sem_open success.\n", semaphores[i]->name);
        // }
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
        printf("[Frequency Server] Requesting new client to respond...\n\n");
        sem_post(s_client.sem); 
        
        // TODO: Client writes to param shm
        printf("[Frequency Server] Awaiting client response...\n");
        sem_wait(s_server.sem);   

        // TODO: Add Param Buffer for meta data and clear freq range
        printf("[Frequency Server] Processing client frequency data...\n");
        // Debug: Print Shared Memory
        // for (int i = 0; i < 6; i += 2) printf("    %d + j%d\n", samples_obj.shm_ptr[i], samples_obj.shm_ptr[i+1]);
        // for (int i = 0; i < RESTRICT_NUM; i++) printf("    restrict[%d]: %d\n", i, restrict_obj.shm_ptr[i]);

        // Read data in
        read_sample_shm(temp_samples, samples_obj.shm_ptr);
        read_restrict_shm(restricted_freq, restrict_obj.shm_ptr);
        // TODO: Add the rest of meta data storing (clear_freq_range, etc.)
        printf("[Frequency Server; from Client] \n");

        // Process data for clear freqs and store result
        clear_freq_search(temp_samples, clr_bands, restricted_freq, RESTRICT_NUM);
        // XXX: Store result into radar table here
        printf("[Frequency Server] Clear Freq Bands updated...\n");
        
        // Write data to Shared Memory
        printf("[Frequency Server] Writing clear frequency data to Shared Memory...\n");
        write_clrfreq_shm(clr_bands, clrfreq_obj.shm_ptr);
        if (msync(clrfreq_obj.shm_ptr, CLR_BANDS_SHM_SIZE, MS_SYNC) == -1) {    // Synchronize data writes with program counter
            perror("msync failed");
        }
        sem_post(s_data.sem);
        printf("[Frequency Server] Processed Client successfully...\n");
    }

    return 0;
}
