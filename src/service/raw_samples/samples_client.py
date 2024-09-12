import mmap
import os
import struct
import time
import posix_ipc
import pickle       # To read in pickle test samples
import numpy as np

RESTRICT_FILE = "/home/df/Desktop/PSU-SuperDARN/Freq_Server/utils/misc_param/restrict.dat.inst"

# Used to mimick USRP functionality
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
    ANTENNAS_NUM = 2                                    # 14 for real
    RESTRICT_NUM = 15
    RESTRICT_ELEM_NUM = RESTRICT_NUM * 2                # 2 = start & stop of restrict freq band
    CLR_BAND_MAX = 6
    CLR_BAND_ELEM_NUM = CLR_BAND_MAX * 3                # 3 = start & stop freqs and noise
    
    RESTRICT_SHM_SIZE = RESTRICT_NUM * 4 * 2            # 4 * 2 = int size * (start & stop of restrict freq band)
    SAMPLES_SHM_SIZE = ANTENNAS_NUM * SAMPLES_NUM * 4    # 2x2500 array of integers (each integer is 4 bytes)
    CLR_FREQ_SHM_SIZE = CLR_BAND_ELEM_NUM* 4 
    
    RETRY_ATTEMPTS = 5
    RETRY_DELAY = 2  # seconds
    
    # Shared Memory Object and Semaphores Names
    SAMPLES_SHM_NAME  = "/samples"             # For Data Transmission
    RESTRICT_SHM_NAME = "/restricted_freq"
    CLR_FREQ_SHM_NAME = "/clear_freq"
    
    SEM_SERVER     = "/sem_server"              # For Sync and reserving client and server roles during data transfer 
    SEM_CLIENT     = "/sem_client"
    SEM_DATA       = "/sem_data"                # For multiple data transfers on single instance
    SEM_SAMPLES    = "/sem_samples"             # For Data locking b/w write/reads
    SEM_INIT       = "/sem_init"                # init = initialization
    SEM_CLRFREQ    = "/sem_clrfreq"
    ACTIVE_CLIENTS_SHM_NAME = "/active_clients" # For Debugging
    
    
    def __init__(self):
        # Shared Memory Object and Semaphores
        self.shm_objects = [
            self.create_shm_obj(self.SAMPLES_SHM_NAME ,self.SAMPLES_SHM_SIZE), 
            self.create_shm_obj(self.RESTRICT_SHM_NAME,self.RESTRICT_SHM_SIZE), 
            self.create_shm_obj(self.CLR_FREQ_SHM_NAME,self.CLR_FREQ_SHM_SIZE)
        ]
        self.semaphores = [
            self.create_semaphore(self.SEM_CLIENT),
            self.create_semaphore(self.SEM_SERVER),
            self.create_semaphore(self.SEM_SERVER),
            self.create_semaphore(self.SEM_SAMPLES),
            self.create_semaphore(self.SEM_INIT),
            self.create_semaphore(self.SEM_CLRFREQ),
        ]

        for obj in self.shm_objects:
            obj['shm_fd'] = self.initialize_shared_memory(obj['name'])
        self.sem_server, self.sem_client = self.semaphores[1]['sem'], self.semaphores[0]['sem']
        # self.sem_server, self.sem_client = self.initialize_semaphores()
                
        self.active_clients_fd = None 
        self.initialize_active_clients_counter()
        print("[clearFrequencyService] Done Initializing...\n\n")
        
    def create_shm_obj(self, name, size):
        """ Returns a dictionary containing pre-filled fields for shared memory (SHM) object data.
            
        Args:
            name (string): name of the shared memory object
            size (integer): size of the file/SHM object

        Returns:
            dictionary: Contains commonly referenced info of a shared memory object
                Contains:
                    - 'name' 
                    - 'shm_ptr' or pointer 
                    - 'shm_fd' or file descriptor 
                    - 'size' of the file/SHM object
        """
        return {
            'name': name,
            'shm_ptr': None,
            'shm_fd': -1,
            'size': size
        }
        
    def create_semaphore(self, name):
        return {
            'name': name,
            'sem':  self.initialize_semaphores(name)
        }

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
                print(f"[clearFrequencyService] Created {shm_name}. Shared Memory Object...")
                return shm_fd
            except FileNotFoundError:
                print("[clearFrequencyService] Shared Memory Object not found. Retrying...")
                attempts += 1
                time.sleep(self.RETRY_DELAY)
        print(f"[clearFrequencyService] Failed to initialize {shm_name}. Shared Memory Object after multiple attempts. Exiting.")
        exit(1)

    def initialize_semaphores(self, name):
        """ Initializes Synchronization Semaphores. Attempts to check for already 
            initialized object (from server).

        Returns:
            void: On success, returns tuple of semaphores (sem_server, sem_client).
        """
        attempts = 0
        while attempts < self.RETRY_ATTEMPTS:
            try:
                print(f"[clearFrequencyService] Attempting to initialize Semaphore {name} (Attempt {attempts + 1}/{self.RETRY_ATTEMPTS})...")
                semaphore = posix_ipc.Semaphore(name)
                print(f"[clearFrequencyService] Semaphore {name} ready...")
                return semaphore        # self.sem_server, self.sem_client
            except posix_ipc.ExistentialError:
                print("[clearFrequencyService] Semaphore not found. Retrying...")
            attempts += 1
            time.sleep(self.RETRY_DELAY)
        print(f"[clearFrequencyService] Failed to initialize Semaphore {name} after multiple attempts. Exiting.")
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
        
    # def write_data():
        # If ___ flags are on, trigger its function
        
    def write_m_data(data_shm_fd, data_shm_size, array_data, data_size, complex=False):
        with mmap.mmap(data_shm_fd, data_shm_size, mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m_data:
            if complex:
                # If complex, flatten and separate real and imaginary parts
                flattened_data = np.column_stack((array_data.real.flatten(), array_data.imag.flatten())).flatten()
            else:
                # Otherwise, just flatten
                flattened_data = array_data.flatten()
            
            print("[Frequency Client] new_data len of: ", len(flattened_data))
            print("[Frequency Client] Writing sample data:\n", flattened_data[:10], "...")
            m_data.seek(0)
            m_data.write(struct.pack('i' * data_size, *flattened_data))   
    
                
    def read_m_data(data_shm_fd, data_size, data_shm_size, data_sub_size = 1):
        with mmap.mmap(data_shm_fd, data_shm_size, mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m_data:
            m_data.seek(0)
            read_data = struct.unpack('i' * (data_size * data_sub_size), m_data.read(data_shm_size))
            print("[clearFrequencyService] Data read from Shm: ", read_data[:5], "...")  # Print first 10 integers for brevity
        return read_data
    
    
    def repack_data(read_data, clr_freq = False, data_size = 1, data_sub_size = 1):
        """Repacks read data from SHM into its proper format. Currently repacks
        the following:
        - Clear Frequency

        Args:
            read_data (_type_): _description_
            data_size (_type_): _description_
            data_sub_size (int, optional): _description_. Defaults to 1.
            clr_freq (bool, optional): Interpret as Clear Frequency (returns 
                centerFreq, Noise). Defaults to False.
        """
        packed_data = []
        if clr_freq:
            noise_data = []
            for start_freq, noise, end_freq in zip(read_data[::3], read_data[1::3], read_data[2::3]):
                # Return Center Freq and Noise
                packed_data.append((start_freq + end_freq) / 2)
                noise_data.append(noise)
            return packed_data, noise_data
            
                
    def sendSamples(self, raw_samples, restrict_data):
        """ Waits for client requests, then processes server data, writes client 
            data, and requests server to process new data. When process is 
            terminated, the try/finally block cleans up.
        """
        # Special: Re-initialize ClearFreqService
        if self.active_clients_fd == None and self.clear_shm_fd == None:
            # TODO: Test 
            self.__init__()
        
        # Special: Find missing data
        if restrict_data is None:
            restrict_data = read_restrict_file(RESTRICT_FILE)
        
        # Get in Queue
        active_clients = self.increment_active_clients()
        print(f"[clearFrequencyService] Active clients count: {active_clients}\n")
        
        try:
            # Map shared memory object as m
            # with mmap.mmap(self.shm_objects[0]['shm_fd'], self.shm_objects[0]['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m_samples:
            #     with mmap.mmap(self.shm_objects[1]['shm_fd'], self.shm_objects[1]['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m_restrict:  
                    
            for obj in self.shm_objects:
                obj['shm_ptr'] = mmap.mmap(self.shm_objects[0]['shm_fd'], self.shm_objects[0]['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE)
            
            # Await for a Client Request
            print("[clearFrequencyService] Awaiting Client Request...\n")
            self.sem_client.acquire()
            print("[clearFrequencyService] Acquired Client Request...")
            
            # Read data from shared memory object
            print("[clearFrequencyService] Reading data from Shared Memory...")
            self.shm_objects[0]['shm_ptr'].seek(0)
            # m_restrict.seek(0)
            samples_data = struct.unpack('i' * (self.ANTENNAS_NUM * self.SAMPLES_NUM), self.shm_objects[0]['shm_ptr'].read(self.shm_objects[0]['size']))
            print("[clearFrequencyService] Data read from Samples: ", samples_data[:5], "...")  # Print first 10 integers for brevity
            print("[clearFrequencyService] Done reading data from Shared Memory...")
            
            # Flatten data arrays for SHM writing 
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

            # self.write_m_data()
                
            # Write new data to Shared Memory Object
            print("[Frequency Client] Writng data to Shared Memory...")    
            print("[Frequency Client] new_data len of samples: ", len(new_sample_data))
            print("[Frequency Client] Writing sample data:\n", new_sample_data[:10], "...")
            self.shm_objects[0]['shm_ptr'].seek(0)
            self.shm_objects[0]['shm_ptr'].write(struct.pack('i' * (self.ANTENNAS_NUM * self.SAMPLES_NUM), *new_sample_data))
            
            print("[Frequency Client] Writing restricted freq data:\n", restrict_data[:2], "...")
            self.shm_objects[1]['shm_ptr'].seek(0)
            self.shm_objects[1]['shm_ptr'].write(struct.pack('i' * (self.RESTRICT_NUM * 2), *new_restrict_data))
            print("[Frequency Client] Done writing data to Shared Memory...")
            
            # Request Server 
            print("[clearFrequencyService] Requesting Server Response...")
            # self.sem_server.release()
            self.semaphores[1]['sem'].release()
        except KeyboardInterrupt:
            print("[clearFrequencyService] Keyboard interrupt received. Exiting...")
        finally:
            # Clean up
            active_clients = self.decrement_active_clients()
            print(f"[clearFrequencyService] Active clients count after decrement: {active_clients}")

            # os.close(self.samples_shm_fd)
            for obj in self.shm_objects: 
                os.close(obj['shm_fd'])
            os.close(self.active_clients_fd)
            # self.sem_server.close()
            # self.sem_client.close()
            for sem in self.semaphores:
                sem['sem'].close()

            if active_clients == 0:
                # print("[clearFrequencyService] No active clients remaining; cleaning up shared resources.")
                # self.cleanup_shm()
                print("[clearFrequencyService] No active clients remaining, but not cleaning up shared resources to keep service idle.")
    
    def cleanup_shm(self):
        for obj in self.shm_objects:
            try:
                posix_ipc.unlink_shared_memory(obj['name'])
                print(f"Unlinked shared memory {obj['name']}")
            except posix_ipc.ExistentialError:
                print(f"Shared memory {obj['name']} does not exist")

        try:
            posix_ipc.unlink_shared_memory(self.ACTIVE_CLIENTS_SHM_NAME)
            print(f"Unlinked shared memory {self.ACTIVE_CLIENTS_SHM_NAME}")
        except posix_ipc.ExistentialError:
            print(f"Shared memory {self.ACTIVE_CLIENTS_SHM_NAME} does not exist")

        for sem in self.semaphores:
            try:
                posix_ipc.unlink_semaphore(sem['name'])
                print(f"Unlinked semaphore {sem['name']}")
            except posix_ipc.ExistentialError:
                print(f"Semaphore {sem['name']} does not exist")

        # try:
        #     posix_ipc.unlink_semaphore(self.SEM_SERVER)
        #     print(f"Unlinked semaphore {self.SEM_SERVER}")
        # except posix_ipc.ExistentialError:
        #     print(f"Semaphore {self.SEM_SERVER} does not exist")

        # try:
        #     posix_ipc.unlink_semaphore(self.SEM_CLIENT)
        #     print(f"Unlinked semaphore {self.SEM_CLIENT}")
        # except posix_ipc.ExistentialError:
        #     print(f"Semaphore {self.SEM_CLIENT} does not exist")
                

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