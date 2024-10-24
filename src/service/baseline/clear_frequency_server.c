#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>

#define SMH_NAME "/shared_memory"
#define SHM_SIZE 1024
#define SEM_SERVER "/sem_server"
#define SEM_CLIENT "/sem_client"


int main() {

    // Open Shared Memory Object
    printf("[Frequency Server] Initializing Shared Memory Object...\n");
    int shm_fd = shm_open(SMH_NAME, O_CREAT | O_RDWR, 0666);
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
    char *shm_ptr = mmap(0, SHM_SIZE, PROT_WRITE | PROT_READ, MAP_SHARED, shm_fd, 0);
    if(shm_ptr == MAP_FAILED) {
        perror("[Frequency Server] Memory Mapping failed.\n");
        exit(EXIT_FAILURE);
    } else {
        printf("[Frequency Server] Memory successfully cached...\n");
    }

    // Open Semaphores for synchronization
    printf("[Frequency Server] Openning Communication Semaphores...\n");
    sem_t *sem_server = sem_open(SEM_SERVER, O_CREAT, 0666, 0);
    sem_t *sem_client = sem_open(SEM_CLIENT, O_CREAT, 0666, 0);
    if ( sem_server == SEM_FAILED || sem_client == SEM_FAILED ) {
        perror("[Frequency Server] One or more sem_open failed.\n");
        exit(EXIT_FAILURE);
    } else {
        printf("[Frequency Server] All semaphores have been established...\n");
    }

    printf("[Frequency Server] Done Initializing...\n\n");

    
    // Continuously Send and Receive messages via Shared Memory
    while (1) {
        // Write Message to Shared Memory
        printf("[Frequency Server] Writing frequency data to Shared Memory...\n");
        const char *message = "Frequency data from Frequency Server";
        strcpy(shm_ptr, message); 
        printf("[Frequency Server] Requesting new client to respond...\n\n");
        sem_post(sem_client); // Request Client
        
        // Await Server Request
        printf("[Frequency Server] Awaiting client response...\n");
        sem_wait(sem_server);   

        // Process Client Message
        printf("[Frequency Server] Processing client frequency data...\n");
        printf("[Frequency Server; from Client] %s\n", shm_ptr);
        printf("[Frequency Server] Processed Client response successfully...\n");
    }


    // Clean up
    munmap(shm_ptr, SHM_SIZE);
    close(shm_fd);
    sem_close(sem_server);
    sem_close(sem_client);
    sem_unlink(SMH_NAME);
    sem_unlink(SEM_SERVER);
    sem_unlink(SEM_CLIENT);

    // Prompt exit to terminal  
    printf("[Frequency Server] Main processes and communication terminated.\n"
            "Goodbye.\n");

    return 0;
}