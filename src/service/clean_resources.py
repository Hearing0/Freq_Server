import posix_ipc

# Shared Memory Object and Semaphores
SHM_NAMES = {"/samples", "/restricted_freq", "/clear_freq"}
ACTIVE_CLIENTS_SHM_NAME = "/active_clients"
SEM_SERVER = "/sem_server"
SEM_CLIENT = "/sem_client"

def cleanup_shared_memory_and_semaphores():
    for SHM_NAME in SHM_NAMES:
        try:
            posix_ipc.unlink_shared_memory(SHM_NAME)
            print(f"Unlinked shared memory {SHM_NAME}")
        except posix_ipc.ExistentialError:
            print(f"Shared memory {SHM_NAME} does not exist")

    try:
        posix_ipc.unlink_shared_memory(ACTIVE_CLIENTS_SHM_NAME)
        print(f"Unlinked shared memory {ACTIVE_CLIENTS_SHM_NAME}")
    except posix_ipc.ExistentialError:
        print(f"Shared memory {ACTIVE_CLIENTS_SHM_NAME} does not exist")

    try:
        posix_ipc.unlink_semaphore(SEM_SERVER)
        print(f"Unlinked semaphore {SEM_SERVER}")
    except posix_ipc.ExistentialError:
        print(f"Semaphore {SEM_SERVER} does not exist")

    try:
        posix_ipc.unlink_semaphore(SEM_CLIENT)
        print(f"Unlinked semaphore {SEM_CLIENT}")
    except posix_ipc.ExistentialError:
        print(f"Semaphore {SEM_CLIENT} does not exist")

if __name__ == "__main__":
    cleanup_shared_memory_and_semaphores()
