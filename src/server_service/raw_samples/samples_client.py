import mmap
import os
import struct
import time
import posix_ipc

# Shared Memory Object and Semaphores
SHM_NAME = "/shared_memory"                 # For Data Transmission
SHM_SIZE = 2 * 2500 * 4  # 2x2500 array of integers (each integer is 4 bytes)
SEM_SERVER = "/sem_server"                  # For Synchronization 
SEM_CLIENT = "/sem_client"
ACTIVE_CLIENTS_SHM_NAME = "/active_clients" # For closing on exit of last client

RETRY_ATTEMPTS = 5
RETRY_DELAY = 2  # seconds

def initialize_shared_memory():
    """ Initialize Shared Memory Object for data transmission between Server 
        and Clients. Attempts to check for already initialized object (from 
        server).

    Returns:
        Integer: On success, returns file descriptor of shared memory object.
    """
    attempts = 0
    while attempts < RETRY_ATTEMPTS:
        try:
            print(f"[Frequency Client] Attempting to initialize Shared Memory Object (Attempt {attempts + 1}/{RETRY_ATTEMPTS})...")
            shm_fd = os.open(f"/dev/shm{SHM_NAME}", os.O_RDWR)
            print("[Frequency Client] Created Shared Memory Object...")
            return shm_fd
        except FileNotFoundError:
            print("[Frequency Client] Shared Memory Object not found. Retrying...")
            attempts += 1
            time.sleep(RETRY_DELAY)
    print("[Frequency Client] Failed to initialize Shared Memory Object after multiple attempts. Exiting.")
    exit(1)

def initialize_semaphores():
    """ Initializes Synchronization Semaphores. Attempts to check for already 
        initialized object (from server).

    Returns:
        Tuple: On success, returns tuple of semaphores (sem_server, sem_client).
    """
    attempts = 0
    while attempts < RETRY_ATTEMPTS:
        try:
            print(f"[Frequency Client] Attempting to initialize Semaphores (Attempt {attempts + 1}/{RETRY_ATTEMPTS})...")
            sem_server = posix_ipc.Semaphore(SEM_SERVER)
            sem_client = posix_ipc.Semaphore(SEM_CLIENT)
            print("[Frequency Client] Shared Memory and Semaphores ready...")
            return sem_server, sem_client
        except posix_ipc.ExistentialError:
            print("[Frequency Client] Semaphore(s) not found. Retrying...")
        attempts += 1
        time.sleep(RETRY_DELAY)
    print("[Frequency Client] Failed to initialize Semaphores after multiple attempts. Exiting.")
    exit(1)

def initialize_active_clients_counter():
    attempts = 0
    while attempts < RETRY_ATTEMPTS:
        try:
            print(f"[Frequency Client] Attempting to initialize Active Clients Counter (Attempt {attempts + 1}/{RETRY_ATTEMPTS})...")
            active_clients_fd = os.open(f"/dev/shm{ACTIVE_CLIENTS_SHM_NAME}", os.O_RDWR | os.O_CREAT, 0o666)
            os.ftruncate(active_clients_fd, struct.calcsize('i'))  # Ensure the size of the shared memory object is large enough for an integer
            # Initialize counter to 0 if it's the first time
            with mmap.mmap(active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
                m.seek(0)
                current_value = struct.unpack('i', m.read(struct.calcsize('i')))[0]
                if current_value < 0 or current_value > 1000:  # Arbitrary threshold to detect uninitialized state
                    m.seek(0)
                    m.write(struct.pack('i', 0))
            print("[Frequency Client] Created Active Clients Counter...")
            return active_clients_fd
        except FileNotFoundError:
            print("[Frequency Client] Active Clients Counter not found. Retrying...")
        except PermissionError:
            print("[Frequency Client] Permission error while accessing Active Clients Counter. Retrying...")
        except OSError as e:
            print(f"[Frequency Client] OS error while accessing Active Clients Counter: {e}. Retrying...")
        attempts += 1
        time.sleep(RETRY_DELAY)
    print("[Frequency Client] Failed to initialize Active Clients Counter after multiple attempts. Exiting.")
    exit(1)

def increment_active_clients(active_clients_fd):
    with mmap.mmap(active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
        m.seek(0)
        active_clients = struct.unpack('i', m.read(struct.calcsize('i')))[0]
        active_clients += 1
        m.seek(0)
        m.write(struct.pack('i', active_clients))
        return active_clients

def decrement_active_clients(active_clients_fd):
    with mmap.mmap(active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
        m.seek(0)
        active_clients = struct.unpack('i', m.read(struct.calcsize('i')))[0]
        active_clients -= 1
        m.seek(0)
        m.write(struct.pack('i', active_clients))
        return active_clients

def main():
    """ Waits for client requests, then processes server data, writes client 
        data, and requests server to process new data. When process is 
        terminated, the try/finally block cleans up.
    """
    
    try:
        
        shm_fd = initialize_shared_memory()
        active_clients_fd = initialize_active_clients_counter()
        sem_server, sem_client = initialize_semaphores()
        print("[Frequency Client] Done Initializing...\n\n")
        
        active_clients = increment_active_clients(active_clients_fd)
        print(f"[Frequency Client] Active clients count: {active_clients}\n")
        
        # Map shared memory object as m
        with mmap.mmap(shm_fd, SHM_SIZE, mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
            while True:    
                # Await for a Client Request
                print("[Frequency Client] Awaiting Client Request...\n")
                sem_client.acquire()
                print("[Frequency Client] Acquired Client Request...")
                
                # Read data from shared memory object
                print("[Frequency Client] Reading data from Shared Memory...")
                m.seek(0)
                data = struct.unpack('i' * (2 * 2500), m.read(SHM_SIZE))
                print("[Frequency Client] Data read from Shared Memory:", data[:10], "...")  # Print first 10 integers for brevity
                print("[Frequency Client] Done reading data from Shared Memory...")
                    
                # Write new data to Shared Memory Object
                time.sleep(2)  # Simulate processing time
                m.seek(0)
                new_data = [i * 2 for i in range(2 * 2500)]  # Example: double each value
                m.write(struct.pack('i' * (2 * 2500), *new_data))
                
                # Request Server 
                print("[Frequency Client] Requesting Server Response...")
                sem_server.release()
    except KeyboardInterrupt:
        print("\n[Frequency Client] Keyboard interrupt received. Exiting...")
    finally:
        # Clean up
        active_clients = decrement_active_clients(active_clients_fd)
        print(f"[Frequency Client] Active clients count after decrement: {active_clients}")

        os.close(shm_fd)
        os.close(active_clients_fd)
        sem_server.close()
        sem_client.close()

        # Special: No active clients remaining; unlink shared memory and semaphores 
        if active_clients == 0:
            print("[Frequency Client] No active clients remaining, cleaning up shared resources.")
                        
            resources = [
                (SHM_NAME, posix_ipc.unlink_shared_memory),
                (ACTIVE_CLIENTS_SHM_NAME, posix_ipc.unlink_shared_memory),
                (SEM_SERVER, posix_ipc.unlink_semaphore),
                (SEM_CLIENT, posix_ipc.unlink_semaphore)
            ]
            
            for name, unlink_function in resources:
                try:
                    unlink_function(name)
                    print(f"    {name} unlinked successfully.")
                except posix_ipc.ExistentialError:
                    print(f"    {name} has already been unlinked.")

if __name__ == "__main__":
    main()
