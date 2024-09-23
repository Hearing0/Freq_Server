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
    # TODO: Look into loading Constants by .ini or .env
    # from dotenv import load_dotenv
    # load_dotenv(".env")

    # Debugging Flags
    CLEAN_ON_INACTIVE = False           # Cleans all semaphores and shared memory objects when there are no Active Clients
    
    # Static Constants
    INT_SIZE = 4
    DOUBLE_SIZE = 8
    
    # Shared Memory Object and Semaphores Constants
    SAMPLES_NUM = 2500
    ANTENNAS_NUM = 2                                    # 14 for real
    RESTRICT_NUM = 15
    # RESTRICT_ELEM_NUM = RESTRICT_NUM * 2                # 2 = start & stop of restrict freq band
    CLR_BAND_MAX = 6
    
    SAMPLES_ELEM_NUM =  ANTENNAS_NUM * SAMPLES_NUM
    CLR_RANGE_ELEM_NUM= 2
    RESTRICT_ELEM_NUM = RESTRICT_NUM * 2
    META_ELEM_NUM =     5
    CLR_BANDS_ELEM_NUM = CLR_BAND_MAX * 3                # 2 = start & stop freqs and noise
    
    SAMPLES_SHM_SIZE =    (ANTENNAS_NUM * SAMPLES_NUM * INT_SIZE) 
    CLR_RANGE_SHM_SIZE =  (2 * INT_SIZE)
    FCENTER_SHM_SIZE =    (1 * INT_SIZE)
    BEAM_NUM_SHM_SIZE =   (1 * INT_SIZE)
    SAMPLE_SEP_SHM_SIZE = (1 * DOUBLE_SIZE)
    RESTRICT_SHM_SIZE =   (RESTRICT_NUM * 2 * INT_SIZE)          # 2 = start and end freqs
    META_DATA_SHM_SIZE =  (5 * DOUBLE_SIZE)
    CLR_BANDS_SHM_SIZE =  (CLR_BAND_MAX * INT_SIZE * 3)     # TODO: Round to convert freqs to int again 
    
    RETRY_ATTEMPTS = 5
    RETRY_DELAY = 2  # seconds
    
    # Shared Memory Object and Semaphores Names
    SAMPLES_SHM_NAME =          "/samples"          # For Data Transmission
    CLR_RANGE_SHM_NAME =        "/clear_freq_range"
    FCENTER_SHM_NAME =          "/fcenter"
    BEAM_NUM_SHM_NAME =         "/beam_num"
    SAMPLE_SEP_SHM_NAME =       "/sample_sep"
    RESTRICT_SHM_NAME =         "/restricted_freq"
    META_DATA_SHM_NAME =        "/meta_data"
    CLRFREQ_SHM_NAME =          "/clear_freq"
    ACTIVE_CLIENTS_SHM_NAME =   "/active_clients"   # For Debugging
    
    SAMPLE_PARAM_NUM =      5
    RESTRICT_PARAM_NUM =    2
    PARAM_NUM =             8
    
    SEM_F_CLIENT =  "/sf_client"               # For reserving client and server roles during data transfer
    SEM_F_SERVER =  "/sf_server"               # And for signalling specific data transfers 
    SEM_F_SAMPLES = "/sf_samples"
    SEM_F_INIT =    "/sf_init"           
    SEM_F_CLRFREQ = "/sf_clrfreq"              
    SEM_L_SAMPLES = "/sl_samples"              # For Data locking b/w write/reads
    SEM_L_INIT =    "/sl_init"                 # init = initialization
    SEM_L_CLRFREQ = "/sl_clrfreq"
    
    SEM_NUM =       8
    SL_NUM =        3
    
    # Initialization Check
    isRestricted = False
    
    semaphores = []
    shm_objects = []
    
    def __init__(self):
        # Shared Memory Object and Semaphores
        self.sf_client  = self.create_semaphore(self.SEM_F_CLIENT)
        self.sf_server  = self.create_semaphore(self.SEM_F_SERVER)
        self.sf_samples = self.create_semaphore(self.SEM_F_SAMPLES)
        self.sf_init    = self.create_semaphore(self.SEM_F_INIT)
        self.sf_clrfreq = self.create_semaphore(self.SEM_F_CLRFREQ)
        self.sl_samples = self.create_semaphore(self.SEM_L_SAMPLES)
        self.sl_init    = self.create_semaphore(self.SEM_L_INIT)
        self.sl_clrfreq = self.create_semaphore(self.SEM_L_CLRFREQ)
        self.semaphores = [
            self.sf_client,
            self.sf_server,
            self.sf_samples,
            self.sf_init,
            self.sf_clrfreq,
            self.sl_samples,
            self.sl_init,
            self.sl_clrfreq,
        ]
        self.shm_objects = [
            self.create_shm_obj(self.SAMPLES_SHM_NAME ,     self.SAMPLES_SHM_SIZE   , self.SAMPLES_ELEM_NUM), 
            self.create_shm_obj(self.CLR_RANGE_SHM_NAME,    self.CLR_RANGE_SHM_SIZE , self.CLR_RANGE_ELEM_NUM), 
            self.create_shm_obj(self.FCENTER_SHM_NAME,      self.FCENTER_SHM_SIZE   , ),
            self.create_shm_obj(self.BEAM_NUM_SHM_NAME,     self.BEAM_NUM_SHM_SIZE  , ), 
            self.create_shm_obj(self.SAMPLE_SEP_SHM_NAME,   self.SAMPLE_SEP_SHM_SIZE, ),
            self.create_shm_obj(self.RESTRICT_SHM_NAME,     self.RESTRICT_SHM_SIZE  , self.RESTRICT_ELEM_NUM), 
            self.create_shm_obj(self.META_DATA_SHM_NAME,    self.META_DATA_SHM_SIZE , self.META_ELEM_NUM),
            self.create_shm_obj(self.CLRFREQ_SHM_NAME,      self.CLR_BANDS_SHM_SIZE , self.CLR_BANDS_ELEM_NUM), 
        ]

        for obj in self.shm_objects:
            obj['shm_fd'] = self.initialize_shared_memory(obj['name'])
        
        # self.sem_server, self.sem_client = self.semaphores[1]['sem'], self.semaphores[0]['sem']
        # self.sem_server, self.sem_client = self.initialize_semaphores()
                
        self.active_clients_fd = None 
        self.initialize_active_clients_counter()
        print("[clearFrequencyService] Done Initializing...\n\n")
        
    def create_shm_obj(self, name, size, elem_num=1):
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
            'size': size,
            'elem_num': elem_num
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
                print(f"[clearFrequencyService] Created {shm_name} Shared Memory Object...")
                return shm_fd
            except FileNotFoundError:
                print("[clearFrequencyService] Shared Memory Object not found. Retrying...")
                attempts += 1
                time.sleep(self.RETRY_DELAY)
        print(f"[clearFrequencyService] Failed to initialize {shm_name} Shared Memory Object after multiple attempts. Exiting.")
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
        
           
    def write_m_data(self, obj, array_data, complex=False):
        try:    
            if complex:
                # If complex, flatten and separate real and imaginary parts
                # flattened_data = np.column_stack((array_data.real.flatten(), array_data.imag.flatten())).flatten()
                for antenna_sample in array_data:  
                    # print("[Frequency Client] sample_arr len: ", len(antenna_sample))
                    
                    # TODO: Time trial w/ zip(sample.real, sample.imag)
                    for sample in antenna_sample:
                        flattened_data.append(int(sample.real))
                        flattened_data.append(int(sample.imag))
            else:
                # Otherwise, just flatten
                flattened_data = []
                if type(array_data[0]) is list:  
                    for row in array_data:
                        flattened_data += row
                else: 
                    flattened_data = array_data
                
            print("[Frequency Client] new_data len of: ", len(flattened_data))
            print("[Frequency Client] Writing sample data:\n", flattened_data[:10], "...")
            obj['shm_ptr'].seek(0)
            obj['shm_ptr'].write(struct.pack('i' * obj['elem_size'], *flattened_data))   
            
        except AttributeError:
            print("[Frequency Client] ERROR: Element Size is incorrect. send()'s parameters were likely not assigned properly. Please verify...")
            
            # print("[Frequency Client] new_data len of: ", len(array_data))
            # print("[Frequency Client] Writing sample data:\n", array_data[:10], "...")
            # obj['shm_ptr'].seek(0)
            # obj['shm_ptr'].write(struct.pack('i' * len(flattened_data), *array_data))   
  
        # print("[Frequency Client] new_data len of: ", len(flattened_data))
        # print("[Frequency Client] Writing sample data:\n", flattened_data[:10], "...")
        # obj['shm_ptr'].seek(0)
        # obj['shm_ptr'].write(struct.pack('i' * obj['elem_num'], *flattened_data))   
                    
    def read_m_data(self, obj):
        """Reads in data from the shared memory file descriptor.

        Args:
            obj (dict): Data object containing file descriptor and shared memory size.
            data_size (int): Shared Memory Size.
            data_sub_size (int, optional): Shared Memory sub element size; used for 2D arrays. Defaults to 1.

        Returns:
            list: Contains 1D list of Shared Memory data.
        """
        obj['shm_ptr'].seek(0)
        read_data = struct.unpack('i' * obj['elem_num'], obj['shm_ptr'].read(obj['size']))
        
        # Debug: Verify format of data object's raw data
        # print("[clearFrequencyService] Data read from Shm: ", read_data[:5], "...")  # Print first 10 integers for brevity
        
        return read_data
    
    
    def repack_data(self, read_data, clr_freq = False, data_size = 1, data_sub_size = 1):
        """Repacks read data from SHM into its specified format. Currently repacks
        the following:
        - Clear Frequency

        Args:
            read_data (list): 1D list of Shared Memory data.
            data_size (int): Shared Memory Size.
            data_sub_size (int, optional): Shared Memory sub element size; used for 2D arrays. Defaults to 1.
            clr_freq (bool, optional): Interpret as Clear Frequency Flag (returns 
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
            
                
    def sendSamples(self, raw_samples, clr_range=None, fcenter=None, beam_num=None, sample_sep=None, restrict_data=None, meta_data=None):
        """ Waits for client requests, then processes server data, writes client 
            data, and requests server to process new data. When process is 
            terminated, the try/finally block cleans up.
        """
        
        input_data = [
            raw_samples, 
            clr_range, 
            fcenter, 
            beam_num,
            sample_sep, 
            restrict_data, 
            meta_data,
        ]
        
        # Special: Re-initialize ClearFreqService
        # if self.active_clients_fd == None:
        #     # TODO: Test 
        #     self.__init__()
        
        # Special: Find missing data
        if self.isRestricted is False and restrict_data is None:
            restrict_data = read_restrict_file(RESTRICT_FILE)
        
        # Get in Queue
        active_clients = self.increment_active_clients()
        print(f"[clearFrequencyService] Active clients count: {active_clients}\n")
        
        try:
            # Map shared memory object pointers
            for obj in self.shm_objects:
                obj['shm_ptr'] = mmap.mmap(obj['shm_fd'], obj['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE)
            
            # Await for a Client Request
            print("[clearFrequencyService] Awaiting Client Request...\n")
            self.sf_client['sem'].acquire()
            print("[clearFrequencyService] Acquired Client Request...")
            
            if restrict_data is not None:
                self.sl_init['sem'].acquire()
                
                if restrict_data is not None:    
                    new_restrict_data = []
                    # print("[Frequency Client] restrict_arr len: ", len(restrict_data))
                    for restricted_freq in restrict_data:   
                        new_restrict_data.append(int(restricted_freq[0]))
                        new_restrict_data.append(int(restricted_freq[1]))
                        
                    print("[Frequency Client] Writing restricted freq data:\n", restrict_data[:2], "...")
                    self.shm_objects[5]['shm_ptr'].seek(0)
                    self.shm_objects[5]['shm_ptr'].write(struct.pack('i' * (self.RESTRICT_NUM * 2), *new_restrict_data))
                    
                self.sl_init['sem'].release()
                self.sf_init['sem'].release()
                                
            
            if raw_samples is not None:
                print("[clearFrequencyService] Awaiting Sample Semphore Lock...")
                self.sl_samples['sem'].acquire()
                
                # self.write_m_data(self.shm_objects[0]['shm_ptr'], raw_samples[:1], (self.ANTENNAS_NUM * self.SAMPLES_NUM), complex = True)
                
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
                        
                
                # Write Samples to Shared Memory Object
                print("[Frequency Client] Writng data to Shared Memory...")    
                print("[Frequency Client] new_data len of samples: ", len(new_sample_data))
                print("[Frequency Client] Writing sample data:\n", new_sample_data[:10], "...")
                self.shm_objects[0]['shm_ptr'].seek(0)
                self.shm_objects[0]['shm_ptr'].write(struct.pack('i' * (self.ANTENNAS_NUM * self.SAMPLES_NUM), *new_sample_data))

                # If given, write that data
                for i in range(self.SAMPLE_PARAM_NUM):
                    print(f"[Frequency Client] Data Write Progress: {i + 1}/{self.SAMPLE_PARAM_NUM} {self.shm_objects[i]['name']}")
                    if i != 0: 
                        if input_data[i] is not None:
                            self.write_m_data(self.shm_objects[i], input_data[i])
                
                                    
                # print("[Frequency Client] Writing restricted freq data:\n", restrict_data[:2], "...")
                # self.shm_objects[5]['shm_ptr'].seek(0)
                # self.shm_objects[5]['shm_ptr'].write(struct.pack('i' * (self.RESTRICT_NUM * 2), *new_restrict_data))
                
                self.sl_samples['sem'].release()
                print("[Frequency Client] Done writing data to Shared Memory...")
                
                # Request Server 
                print("[clearFrequencyService] Requesting Server Response...")
                self.sf_server['sem'].release()
                # self.semaphores[1]['sem'].release()
                
                # Read-in Clear Freq data
                print("[clearFrequencyService] Awaiting Server Response...")
                self.sf_clrfreq['sem'].acquire()
                self.sl_clrfreq['sem'].acquire()
                print("[clearFrequencyService] Recieved Server Response. Reading Clear Freq data...")
                new_noise_data = []
                new_clrfreq_data = self.read_m_data(self.shm_objects[7])
                new_clrfreq_data, new_noise_data = self.repack_data(new_clrfreq_data, True)
                for clr_freq in zip(new_clrfreq_data, new_noise_data):
                    print(f"[clearFrequencyService] Clear Freq Band: | {clr_freq[0]} (Hz), {clr_freq[1]} (N/A) |")
                    
                self.sl_clrfreq['sem'].release()
                        
        except KeyboardInterrupt:
            print("[clearFrequencyService] Keyboard interrupt received. Exiting...")
        finally:
            # Clean up
            active_clients = self.decrement_active_clients()
            print(f"[clearFrequencyService] Active clients count after decrement: {active_clients}")

            # for obj in self.shm_objects: 
            #     os.close(obj['shm_fd'])
            # os.close(self.active_clients_fd)
            # for sem in self.semaphores:
            #     sem['sem'].close()
            

            if active_clients == 0:
                print("[clearFrequencyService] No active clients remaining; cleaning up shared resources.")
                # self.cleanup_shm()
                # print("[clearFrequencyService] No active clients remaining, but not cleaning up shared resources to keep service idle.")
    
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
CFS.sendSamples(raw_samples, restrict_data=read_restrict_file(RESTRICT_FILE))
# CFS.sendSamples(raw_samples)