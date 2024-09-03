import mmap
import os
import struct
import time
import posix_ipc
import pickle       # To read in pickle test samples

RESTRICT_FILE = "/home/df/Desktop/PSU-SuperDARN/Freq_Server/utils/misc_param/restrict.dat.inst"

def read_restrict_file(restrict_file):
    restricted_frequencies = []
    with open(restrict_file, 'r') as f:
        for line in f:
            if line[0] == '#' or line[0] == 'd' or len(line) < 8:
                continue
            line = line.split(' ')
            restrict_start = int(line[0]) * 1e3 # convert kHz units in restrict to Hz
            restrict_end = int(line[1]) * 1e3 # convert kHz units in restrict to Hz
            restricted_frequencies.append([restrict_start, restrict_end])

    return restricted_frequencies; 




class ClearFrequencyService():
    # Shared Memory Object and Semaphores
    SAMPLES_NUM = 2500
    ANTENNAS_NUM = 2
    RESTRICT_NUM = 15
    RESTRICT_ELEM_NUM = RESTRICT_NUM * 2                # 2 = start & stop of restrict freq band
    CLR_BAND_MAX = 6
    RESTRICT_SHM_SIZE = RESTRICT_NUM * 4 * 2            # 4 * 2 = int size * (start & stop of restrict freq band)

    SAMPLES_SHM_SIZE = ANTENNAS_NUM * SAMPLES_NUM * 4    # 2x2500 array of integers (each integer is 4 bytes)
    
    RETRY_ATTEMPTS = 5
    RETRY_DELAY = 2  # seconds
    
    def __init__(self):
        # Shared Memory Object and Semaphores
        self.RESTRICT_SHM_SIZE = self.RESTRICT_NUM * 4 * 2
        self.SAMPLES_SHM_SIZE = self.ANTENNAS_NUM * self.SAMPLES_NUM * 4  # 14x2500 array of integers (each integer is 4 bytes)
        self.SAMPLES_SHM_NAME = "/samples"             # For Data Transmission
        self.RESTRICT_SHM_NAME = "/restricted_freq"
        self.SEM_SERVER = "/sem_server"                              # For Synchronization 
        self.SEM_CLIENT = "/sem_client"
        self.ACTIVE_CLIENTS_SHM_NAME = "/active_clients" # For Debugging
        
        self.samples_shm_fd = self.initialize_shared_memory(self.SAMPLES_SHM_NAME)
        self.restrict_shm_fd = self.initialize_shared_memory(self.RESTRICT_SHM_NAME)
        self.active_clients_fd = None 
        self.initialize_active_clients_counter()
        self.sem_server, self.sem_client = self.initialize_semaphores()
        print("[clearFrequencyService] Done Initializing...\n\n")

    def initialize_shared_memory(self, shm_name):
        """ Initialize Shared Memory Object for data transmission between Server 
            and Clients. Attempts to check for already initialized object (from 
            server).

        Returns:
            Integer: On success, returns file descriptor of shared memory object.
        """
        attempts = 0
        while attempts < self.RETRY_ATTEMPTS:
            try:
                print(f"[clearFrequencyService] Attempting to initialize {shm_name} Shared Memory Object (Attempt {attempts + 1}/{self.RETRY_ATTEMPTS})...")
                shm_fd = os.open(f"/dev/shm{shm_name}", os.O_RDWR)
                print("[clearFrequencyService] Created {shm_name} Shared Memory Object...")
                return shm_fd
            except FileNotFoundError:
                print("[clearFrequencyService] Shared Memory Object not found. Retrying...")
                attempts += 1
                time.sleep(self.RETRY_DELAY)
        print("[clearFrequencyService] Failed to initialize {shm_name} Shared Memory Object after multiple attempts. Exiting.")
        exit(1)

    def initialize_semaphores(self):
        """ Initializes Synchronization Semaphores. Attempts to check for already 
            initialized object (from server).

        Returns:
            void: On success, returns tuple of semaphores (sem_server, sem_client).
        """
        attempts = 0
        while attempts < self.RETRY_ATTEMPTS:
            try:
                print(f"[clearFrequencyService] Attempting to initialize Semaphores (Attempt {attempts + 1}/{self.RETRY_ATTEMPTS})...")
                self.sem_server = posix_ipc.Semaphore(self.SEM_SERVER)
                self.sem_client = posix_ipc.Semaphore(self.SEM_CLIENT)
                print("[clearFrequencyService] Shared Memory and Semaphores ready...")
                return self.sem_server, self.sem_client
            except posix_ipc.ExistentialError:
                print("[clearFrequencyService] Semaphore(s) not found. Retrying...")
            attempts += 1
            time.sleep(self.RETRY_DELAY)
        print("[clearFrequencyService] Failed to initialize Semaphores after multiple attempts. Exiting.")
        exit(1)

    def initialize_active_clients_counter(self):
        attempts = 0
        while attempts < self.RETRY_ATTEMPTS:
            try:
                print(f"[clearFrequencyService] Attempting to initialize Active Clients Counter (Attempt {attempts + 1}/{self.RETRY_ATTEMPTS})...")
                self.active_clients_fd = os.open(f"/dev/shm{self.ACTIVE_CLIENTS_SHM_NAME}", os.O_RDWR | os.O_CREAT, 0o666)
                os.ftruncate(self.active_clients_fd, struct.calcsize('i'))  # Ensure the size of the shared memory object is large enough for an integer
                # Initialize counter to 0 if it's the first time
                with mmap.mmap(self.active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
                    m.seek(0)
                    current_value = struct.unpack('i', m.read(struct.calcsize('i')))[0]
                    if current_value < 0 or current_value > 1000:  # Arbitrary threshold to detect uninitialized state
                        m.seek(0)
                        m.write(struct.pack('i', 0))
                print("[clearFrequencyService] Created Active Clients Counter...")
                return 
            except FileNotFoundError:
                print("[clearFrequencyService] Active Clients Counter not found. Retrying...")
            except PermissionError:
                print("[clearFrequencyService] Permission error while accessing Active Clients Counter. Retrying...")
            except OSError as e:
                print(f"[clearFrequencyService] OS error while accessing Active Clients Counter: {e}. Retrying...")
            attempts += 1
            time.sleep(self.RETRY_DELAY)
        print("[clearFrequencyService] Failed to initialize Active Clients Counter after multiple attempts. Exiting.")
        exit(1)

    def increment_active_clients(self):
        with mmap.mmap(self.active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
            m.seek(0)
            active_clients = struct.unpack('i', m.read(struct.calcsize('i')))[0]
            active_clients += 1
            m.seek(0)
            m.write(struct.pack('i', active_clients))
            print(f"[clearFrequencyService] Incremented Active Clients Counter: {active_clients}")
            return active_clients

    def decrement_active_clients(self):
        with mmap.mmap(self.active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
            m.seek(0)
            active_clients = struct.unpack('i', m.read(struct.calcsize('i')))[0]
            active_clients -= 1
            m.seek(0)
            m.write(struct.pack('i', active_clients))
            print(f"[clearFrequencyService] Decremented Active Clients Counter: {active_clients}")
            return active_clients
                
    def sendSamples(self, raw_samples, restrict_data):
        """ Waits for client requests, then processes server data, writes client 
            data, and requests server to process new data. When process is 
            terminated, the try/finally block cleans up.
        """
        
        # Special: Find missing data
        if restrict_data is None:
            restrict_data = read_restrict_file(RESTRICT_FILE)
        
        # Get in Queue
        active_clients = self.increment_active_clients()
        print(f"[clearFrequencyService] Active clients count: {active_clients}\n")
        
        try:
            # Map shared memory object as m
            # TODO: Why does nested map break!!!!!!!!!!!! 
            with (
                mmap.mmap(self.samples_shm_fd, self.SAMPLES_SHM_SIZE, mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m_samples,
                mmap.mmap(self.restrict_shm_fd, self.RESTRICT_SHM_SIZE, mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m_restrict
                ):  
                # Await for a Client Request
                print("[clearFrequencyService] Awaiting Client Request...\n")
                self.sem_client.acquire()
                print("[clearFrequencyService] Acquired Client Request...")
                
                # Read data from shared memory object
                print("[clearFrequencyService] Reading data from Shared Memory...")
                m_samples.seek(0)
                # m_restrict.seek(0)
                samples_data = struct.unpack('i' * (self.ANTENNAS_NUM * self.SAMPLES_NUM), m_samples.read(self.SAMPLES_SHM_SIZE))
                print("[clearFrequencyService] Data read from Samples:", samples_data[:5], "...")  # Print first 10 integers for brevity
                # TODO: Fix Fail check
                # if restrict_data is not None:
                #     print("[clearFrequencyService] Data read from Restricted Freqs:", restrict_data[:5], "...")  # Print first 10 integers for brevity
                print("[clearFrequencyService] Done reading data from Shared Memory...")
                   
                # Repack read data into its orignal format 
                print("[Frequency Client] Repacking data to Shared Memory...")    
                trimmed_samples = raw_samples[:1]           #HACK: writes only first antenna's samples; write all antenna samples
                print("[Frequency Client] trimmed len: ", len(trimmed_samples))
                new_sample_data = []
                for antenna_sample in trimmed_samples:  
                    print("[Frequency Client] sample_arr len: ", len(antenna_sample))
                    for sample in antenna_sample:
                        new_sample_data.append(int(sample.real))
                        new_sample_data.append(int(sample.imag))
                        
                        # Debug: Display samples
                        # print("[Frequency Client] Flattening samples: ", sample)
                        # print("[Frequency Client]                   : ", new_data[-2])
                        # print("[Frequency Client]                   : ", new_data[-1])
                new_restrict_data = []
                print("[Frequency Client] restrict_arr len: ", len(restrict_data))
                for restricted_freq in restrict_data:   
                    new_restrict_data.append(int(restricted_freq[0]))
                    new_restrict_data.append(int(restricted_freq[1]))
                # print(new_restrict_data)
                    
                # Write new data to Shared Memory Object
                print("[Frequency Client] Writng data to Shared Memory...")    
                print("[Frequency Client] new_data len of samples: ", len(new_sample_data))
                print("[Frequency Client] Writing sample data:\n", new_sample_data[:10], "...")
                m_samples.seek(0)
                m_samples.write(struct.pack('i' * (self.ANTENNAS_NUM * self.SAMPLES_NUM), *new_sample_data))
                print("[Frequency Client] Writing restricted freq data:\n", restrict_data[:2], "...")
                
                m_restrict.seek(0)
                m_restrict.write(struct.pack('i' * (self.RESTRICT_NUM * 2), *new_restrict_data))
                print("[Frequency Client] Done writing data to Shared Memory...")
                
                # Request Server 
                print("[clearFrequencyService] Requesting Server Response...")
                self.sem_server.release()
        except KeyboardInterrupt:
            print("[clearFrequencyService] Keyboard interrupt received. Exiting...")
        finally:
            # Clean up
            active_clients = self.decrement_active_clients()
            print(f"[clearFrequencyService] Active clients count after decrement: {active_clients}")

            os.close(self.samples_shm_fd)
            os.close(self.restrict_shm_fd)
            os.close(self.active_clients_fd)
            self.sem_server.close()
            self.sem_client.close()

            if active_clients == 0:
                print("[clearFrequencyService] No active clients remaining, but not cleaning up shared resources to keep service idle.")
                

def read_sample_pickle(pickle_file):
    """ Reads in raw sample data and sample meta data from pickle for a Python 
    script. 

    Args:
        pickle_file (string): read filepath for pickle file

    Returns:
        raw_samples: antenna_num by sample_num by complex (3D) array
    """
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
        
    raw_samples = data['raw_samples']
    sample_meta = data['sample_data']
    
    return raw_samples, sample_meta



# if __name__ == "__main__":
    # main()



raw_samples, meta_data = read_sample_pickle("/home/df/Desktop/PSU-SuperDARN/Freq_Server/utils/pickle_input/clrfreq_dump.1.pickle")

# Send Empty sample info
CFS = ClearFrequencyService()
CFS.sendSamples(raw_samples, read_restrict_file(RESTRICT_FILE))