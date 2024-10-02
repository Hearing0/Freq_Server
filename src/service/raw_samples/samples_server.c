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
#include <time.h>
#include "../../clear_freq_search.c"

// #define INT_SIZE = 4
// #define DOUBLE_SIZE = 8


#define SAMPLES_NUM 2500
#define ANTENNAS_NUM 14
#define BEAM_NUM_ELEM 1
#define RESTRICT_NUM 15 //16             // Number of restricted freq bands in the restrict.dat.inst
#ifndef CLR_BANDS_MAX
#define CLR_BANDS_MAX 6
#endif
#define CLR_STORAGE_NUM 5
#define CLR_STORE_FILEPATH "../../../utils/csv_dump/clr_band_storage/"

#define SAMPLES_SHM_SIZE    (ANTENNAS_NUM * SAMPLES_NUM * sizeof(int)) // TODO: Finish 2x2500 test
#define CLR_RANGE_SHM_SIZE  (2 * sizeof(double))
#define FCENTER_SHM_SIZE    (1 * sizeof(int))
#define BEAM_NUM_SHM_SIZE   (1 * sizeof(int))
#define SAMPLE_SEP_SHM_SIZE (1 * sizeof(double))
#define RESTRICT_SHM_SIZE   (RESTRICT_NUM * 2 * sizeof(int))          // 2 = start and end freqs
#define META_DATA_SHM_SIZE  (5 * sizeof(double))
#define CLR_BANDS_SHM_SIZE  (CLR_BANDS_MAX * 4 * 3)     // TODO: Round to convert freqs to int again 

// Shared Memory and Semaphore Names 
#define SAMPLES_SHM_NAME        "/samples"
#define CLR_RANGE_SHM_NAME      "/clear_freq_range"
#define FCENTER_SHM_NAME        "/fcenter"
#define BEAM_NUM_SHM_NAME       "/beam_num"
#define SAMPLE_SEP_SHM_NAME     "/sample_sep"
#define RESTRICT_SHM_NAME       "/restricted_freq"
#define META_DATA_SHM_NAME      "/meta_data"
#define CLRFREQ_SHM_NAME        "/clear_freq"
#define SAMPLE_PARAM_NUM 5
#define RESTRICT_PARAM_NUM 2
#define PARAM_NUM 8

#define SEM_F_CLIENT    "/sf_client"               // For Sync and reserving client and server roles during data transfer
#define SEM_F_SERVER    "/sf_server"    
#define SEM_F_SAMPLES   "/sf_samples"
#define SEM_F_INIT      "/sf_init"           
#define SEM_F_CLRFREQ   "/sf_clrfreq"              // For multiple data transfers on single instance 
#define SEM_L_SAMPLES   "/sl_samples"              // For Data locking b/w write/reads
#define SEM_L_INIT      "/sl_init"                 // init = initialization
#define SEM_L_CLRFREQ   "/sl_clrfreq"
#define SL_NUM 3
#define SEM_NUM 8

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




// TODO: Read of specific data types 
// TODO: Calc all beam directions at recieve of samples
// TODO: Num of Radar {has table of all beam dirc {has Clr freq in beam direction}}} while sharing restricting assigned freq 
// TODO: Rewrite of usrp sample send/clr freq request timing logic 

// Semaphores Locks prevent race conditions
// Semaphore Flags allow client and server to signal specific data transfers
struct semaphore sf_client   = {SEM_F_CLIENT,  NULL};
struct semaphore sf_server   = {SEM_F_SERVER,  NULL};
struct semaphore sf_samples  = {SEM_F_SAMPLES, NULL};
struct semaphore sf_init     = {SEM_F_INIT,    NULL};
struct semaphore sf_clrfreq  = {SEM_F_CLRFREQ, NULL};
struct semaphore sl_samples  = {SEM_L_SAMPLES, NULL};
struct semaphore sl_init     = {SEM_L_INIT,    NULL};
struct semaphore sl_clrfreq  = {SEM_L_CLRFREQ, NULL};
struct semaphore *semaphores[SEM_NUM] = {
    &sf_client,
    &sf_server,    
    &sf_samples,
    &sf_init,
    &sf_clrfreq,
    &sl_samples,
    &sl_init,
    &sl_clrfreq,
};

shm_obj samples_obj = {SAMPLES_SHM_NAME, NULL, -1, SAMPLES_SHM_SIZE};
shm_obj clr_range_obj = {CLR_RANGE_SHM_NAME, NULL, -1, CLR_RANGE_SHM_SIZE};
shm_obj fcenter_obj = {FCENTER_SHM_NAME, NULL, -1, FCENTER_SHM_SIZE};
shm_obj beam_num_obj = {BEAM_NUM_SHM_NAME, NULL, -1, BEAM_NUM_SHM_SIZE};
shm_obj sample_sep_obj = {SAMPLE_SEP_SHM_NAME, NULL, -1, SAMPLE_SEP_SHM_SIZE};
shm_obj restrict_obj = {RESTRICT_SHM_NAME, NULL, -1, RESTRICT_SHM_SIZE};
shm_obj meta_obj = {META_DATA_SHM_NAME, NULL, -1, META_DATA_SHM_SIZE};
shm_obj clrfreq_obj = {CLRFREQ_SHM_NAME, NULL, -1, CLR_BANDS_SHM_SIZE};
struct shm_obj *objects[PARAM_NUM] = {
    &samples_obj,
    &clr_range_obj,
    &fcenter_obj,
    &beam_num_obj,
    &sample_sep_obj,
    &restrict_obj,
    &meta_obj,
    &clrfreq_obj,
};




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
 * @brief  Reads shm integer data into the result ptr.
 * @note   
 * @param  *result: Destination ptr where data will be stored.
 * @param  *shm_ptr: Shared Memory Pointer where data will read-in from.
 * @param  elem_num: Number of elements to read in from shm_ptr. 
 * @retval None
 */
void read_int(void *result, int *shm_ptr, int elem_num) {
    if (elem_num > 1) {
        int *result_ptr = (int *) result;
        for (int i = 0; i < elem_num; i++) {
            result_ptr[i] = shm_ptr[i];
            printf("    read_int: %d\n", i);
        }
    } else {
        printf("Error: Use read_int_single() for singleton variables\n");
    }
}

void read_int_single(int result, int *shm_ptr) {
    result = shm_ptr[0];
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
void clean_obj(shm_obj shm_obj) {
    // if (shm_obj.shm_ptr) 
    munmap(shm_obj.shm_ptr, shm_obj.size);
    // if (shm_obj.shm_fd >= 0) 
    close(shm_obj.shm_fd);
    sem_unlink(shm_obj.name);
}

void clean_sem(semaphore sem) {
    // if ((sem.sem)) 
    sem_close(sem.sem);
    sem_unlink(sem.name);
}

/**
 * @brief  Deallocates all service semaphores and SHM pointers.
 * @note   
 * @retval None
 */
void cleanup() {
    for (int i = 0; i < SEM_NUM; i++) clean_sem(*semaphores[i]);
    for (int i = 0; i < PARAM_NUM; i++) clean_obj(*objects[i]);
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

void write_clr_log_csv(freq_band **clr_storage, int clr_num) {
    // Timestamp Variables
    time_t raw_time;
    struct tm *time_info;
    int buffer_size = 100;
    char timestamp[buffer_size];
    char name[buffer_size]; 

    // Generate timestamp
    time(&raw_time);
    time_info = localtime(&raw_time);
    strftime(timestamp, buffer_size, "%Y.%m.%d_%H:%M:%S", time_info);
    snprintf(name, sizeof(name), "utils/csv_dump/clr_log/clrlog_%s.csv", timestamp);

    // Generate clear log file
    FILE *file = fopen(name, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Start Frequency,End Frequency,Noise,Clear Freq Start,Clear Freq End\n");
    for (int clr_batch_idx = 0; clr_batch_idx < clr_num; clr_batch_idx++) {
        freq_band *clr_bands = clr_storage[clr_batch_idx];

        // Find Start and End of Clear Freq Range
        int clr_start = RAND_MAX;
        int clr_end = 0;
        for (int i = 0; i < CLR_BANDS_MAX; i++) {
            if (clr_bands[i].f_start < clr_start && clr_bands[i].noise < RAND_MAX) clr_start = clr_bands[i].f_start;
            if (clr_bands[i].f_end > clr_end && clr_bands[i].noise < RAND_MAX) clr_end = clr_bands[i].f_end;
        }    

        // Record each Clear Freq
        for (int i = 0; i < CLR_BANDS_MAX; i++) {
            // Special: Print Clear Freq Range on Line 0
            if (i == 0) fprintf(file, "%d,%d,%f,%d,%d\n", clr_bands[i].f_start, clr_bands[i].f_end, clr_bands[i].noise,clr_start,clr_end);
            else fprintf(file, "%d,%d,%f\n", clr_bands[i].f_start, clr_bands[i].f_end, clr_bands[i].noise);
        }
    }

    fclose(file);
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
        if (i < (SEM_NUM - SL_NUM)) semaphores[i]->sem = sem_open(semaphores[i]->name, O_CREAT, 0666, 0);
        else semaphores[i]->sem = sem_open(semaphores[i]->name, O_CREAT, 0644, 1);
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
    freq_band **clr_bands_storage = NULL;
    clr_bands_storage = (freq_band **)malloc(CLR_BANDS_MAX * sizeof(freq_band *));
    if (clr_bands_storage == NULL) {
        perror("Error allocating memory for clr_bands_storage pointers");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < CLR_STORAGE_NUM; i++) {
        clr_bands_storage[i] = (freq_band *)malloc(CLR_BANDS_MAX * sizeof(freq_band));
        if (clr_bands_storage[i] == NULL) {
            perror("Error allocating memory for clr_bands_storage elements");
            exit(EXIT_FAILURE);
        }
    }
    int clr_storage_i = 0;
    int clr_range[2] = {0};
    int fcenter = 0;
    int beam_num = 0;
    double sample_sep = 0;
    sample_meta_data meta_data = {0};
    // freq_band restricted_freq[RESTRICT_NUM] = {0};
            
    // Read-in Restricted Frequencies
    char *restrict_file = "/home/df/Desktop/PSU-SuperDARN/Freq_Server/utils/misc_param/restrict.dat.inst";
    read_restrict(restrict_file, restricted_freq, RESTRICT_NUM);

    for (int i = 0; i < RESTRICT_NUM; i++) {
        printf("restrict[%d]: %d -- %d\n", i, restricted_freq[i].f_start, restricted_freq[i].f_end);
    }

    // Continuously process clients via shared memory
    while (1) {
        printf("[Frequency Server] Requesting new client to respond...\n\n");
        sem_post(sf_client.sem); 
        
        printf("[Frequency Server] Awaiting client response...\n");
        sem_wait(sf_server.sem);   

        // If initialization data flagged, read and store data
        if (sem_trywait(sf_init.sem)){
            printf("[Frequency Server] Awaiting initialization data unlock...\n");
            sem_wait(sl_init.sem);
            printf("[Frequency Server] Initialization data read...\n");
            // if (restrict_obj.shm_ptr[0] != 0) read_restrict_shm(restricted_freq, restrict_obj.shm_ptr);
            // read(meta_data)
            sem_post(sl_init.sem);
            printf("[Frequency Server] Initialization data read; processing...\n");
            // storeInRadarTable(restrict_freq, meta_data)
            printf("[Frequency Server] Initialization data processed...\n");
        }
        // If samples flagged (& intialized), process clear frequency
        if (sem_trywait(sf_samples.sem) && restricted_freq != NULL){
            // Wait to read-in data
            printf("[Frequency Server] Awaiting sample data unlock...\n");
            sem_wait(sl_samples.sem);

            // Process Sample relevant data
            printf("[Frequency Server] Processing client sample data...\n");
            read_sample_shm(temp_samples, samples_obj.shm_ptr);
            printf("[Frequency Server] Samples done...\n");
            
            if (clr_range_obj.shm_ptr[0] != 0) read_int(clr_range, clr_range_obj.shm_ptr, 2);
            printf("clr_range: %d -- %d\n", clr_range[0], clr_range[1]);
            printf("[Frequency Server] Clear Range done...\n");

            if (fcenter_obj.shm_ptr[0] != 0) {
                // int *fcenter_ptr;
                read_int_single(fcenter, fcenter_obj.shm_ptr);
                // fcenter = *fcenter_ptr;
                printf("fcenter: %d\n", fcenter);
            }
            printf("[Frequency Server] Freq Center done...\n");


            if (beam_num_obj.shm_ptr[0] != 0) {
                read_int_single(beam_num, beam_num_obj.shm_ptr);
                printf("beam_num: %d\n", beam_num);
            }
            printf("[Frequency Server] Beam Number done...\n");

            sem_post(sl_samples.sem);

            // store sample data for debugging
            // debug function here

            for (int i = 0; i < RESTRICT_NUM; i++) {
                printf("restrict[%d]: %d -- %d\n", i, restricted_freq[i].f_start, restricted_freq[i].f_end);
            }

            // Process Clear Freq
            clear_freq_search(
                temp_samples, 
                clr_range,
                restricted_freq, 
                RESTRICT_NUM,
                clr_bands                
            );
            // update_clr_table(clr_bands);
            
            // Write Clear Freq Data
            printf("[Frequency Server] Writing clear frequency data to Shared Memory...\n");
            sem_wait(sl_clrfreq.sem);
            write_clrfreq_shm(clr_bands, clrfreq_obj.shm_ptr);
            if (msync(clrfreq_obj.shm_ptr, CLR_BANDS_SHM_SIZE, MS_SYNC) == -1) {    // Synchronize data writes with program counter
                perror("msync failed");
            }
            printf("[Frequency Server] clrfreq_shm written...\n");
            sem_post(sl_clrfreq.sem);
            sem_post(sf_clrfreq.sem);
            printf("[Frequency Server] Processed Clear Freq Request successfully...\n");

            // Debug: Store prior clear freq band sets
            clr_bands_storage[clr_storage_i] = clr_bands;
            clr_storage_i++;
            if (clr_storage_i >= CLR_STORAGE_NUM) {
                write_clr_log_csv(clr_bands_storage, clr_storage_i);
                clr_storage_i = 0;
            }
            printf("[Frequency Server] Clr Freq Log Batch: %d/%d\n", clr_storage_i, CLR_STORAGE_NUM);
        }
        else if (restricted_freq == NULL) {
            printf("[Frequency Server] ERROR: Called for Clear Freq without prior Initialization\n");
        }
        
        printf("[Frequency Server] Processed Client successfully...\n");

    }

    return 0;
}
