import mmap
import os
import struct
import time
import timeit
import posix_ipc
import pickle       # To read in pickle test samples
import numpy as np




class ClearFrequencyService():
    # TODO: Look into loading Constants by .ini or .env
    # from dotenv import load_dotenv
    # load_dotenv(".env")

    # Program Flags
    CLEAN_ON_INACTIVE   = False           # Cleans all semaphores and shared memory objects when there are no Active Clients
    soft_kill = False
    
    # Static Constants
    CHAR_SIZE = 1
    INT_SIZE = 4
    DOUBLE_SIZE = 8
    
    # Shared Memory Object and Semaphores Constants
    SAMPLES_NUM  = 2500
    ANTENNA_NUM = 16
    RESTRICT_NUM = 20
    META_ELEM    = 3                                    # 3 = 4 - 1 (fcenter has unique obj)
    CLR_BAND_MAX = 6
    
    SAMPLES_ELEM_NUM    = ANTENNA_NUM * SAMPLES_NUM * 2
    CLR_RANGE_ELEM_NUM  = 2
    RESTRICT_ELEM_NUM   = RESTRICT_NUM * 2
    META_ELEM_NUM       = META_ELEM + ANTENNA_NUM
    CLR_BANDS_ELEM_NUM  = 1 * 3                         # 3     = start & stop freqs and noise
    SITE_ID_ELEM_NUM    = 1 * 3                         # 1 * 3 = one instance of a 3 letter identifier
    
    SAMPLES_SHM_SIZE        = (ANTENNA_NUM * SAMPLES_NUM * 2 * INT_SIZE) 
    CLR_RANGE_SHM_SIZE      = (2 * INT_SIZE)
    FCENTER_SHM_SIZE        = (1 * INT_SIZE)
    BEAM_NUM_SHM_SIZE       = (1 * INT_SIZE)
    SAMPLE_SEP_SHM_SIZE     = (1 * INT_SIZE)
    RESTRICT_SHM_SIZE       = (RESTRICT_NUM * 2 * INT_SIZE)          # 2 = start and end freqs
    META_DATA_SHM_SIZE      = ((META_ELEM + ANTENNA_NUM) * DOUBLE_SIZE)
    ANTENNA_SHM_SIZE        = (1 * INT_SIZE)
    CLR_BANDS_SHM_SIZE      = (1 * INT_SIZE * 3)     # TODO: Round to convert freqs to int again 
    SITE_ID_SHM_SIZE        = (3 * CHAR_SIZE)
    RADAR_ID_SHM_SIZE       = (1 * INT_SIZE)
    ACTIVE_CLIENTS_SHM_SIZE = (1 * INT_SIZE)

    RETRY_ATTEMPTS = 3
    RETRY_DELAY = 2  # seconds
    
    # Shared Memory Object and Semaphores Names
    SAMPLES_SHM_NAME =          "/samples"          # For Data Transmission
    CLR_RANGE_SHM_NAME =        "/clear_freq_range"
    FCENTER_SHM_NAME =          "/fcenter"
    BEAM_NUM_SHM_NAME =         "/beam_num"
    SAMPLE_SEP_SHM_NAME =       "/sample_sep"
    RESTRICT_SHM_NAME =         "/restricted_freq"
    META_DATA_SHM_NAME =        "/meta_data"
    ANTENNA_SHM_NAME =          "/antenna_num"
    CLRFREQ_SHM_NAME =          "/clear_freq"
    SITE_ID_SHM_NAME =          "/site_id"
    RADAR_ID_SHM_NAME =         "/radar_id"
    ACTIVE_CLIENTS_SHM_NAME =   "/active_clients"   # For Debugging

    # Semaphore Constants
    SAMPLE_PARAM_NUM =      2
    RESTRICT_PARAM_NUM =    2
    PARAM_NUM =             10
    
    SEM_F_CLIENT =      "/sf_client"               # For reserving client and server roles during data transfer
    SEM_F_SERVER =      "/sf_server"               # And for signalling specific data transfers 
    SEM_F_SAMPLES =     "/sf_samples"
    SEM_F_INIT =        "/sf_init"           
    SEM_F_CLRFREQ =     "/sf_clrfreq"              
    SEM_F_PROCESSED =   "/sf_processed"            # For processed data transfer
    SEM_L_SAMPLES =     "/sl_samples"              # For Data locking b/w write/reads
    SEM_L_INIT =        "/sl_init"                 # init = initialization
    SEM_L_CLRFREQ =     "/sl_clrfreq"
    
    SEM_NUM =       9
    SL_NUM =        3
    
    # Service Variables
    semaphores = []
    shm_objects = []
    cur_antenna_num = ANTENNA_NUM
    old_meta_data = [[], 0, 0.0, 0.0]
    old_smsep = 0
    
    
    def __init__(self, sid = 'lab'):
        # Process Site ID during Sample Send 
        ClearFrequencyService.sid = sid
        
        try:
            # Skip Initialization if SHMs exists
            if (len(ClearFrequencyService.semaphores) > 0 and len(ClearFrequencyService.shm_objects) > 0):
                print("[clearFrequencyService] Existing Shared Memory Objects and Semaphores found. Skipping Initialization...")
                return
            
            # Shared Memory Object and Semaphores
            ClearFrequencyService.sf_client      = self.create_semaphore(self.SEM_F_CLIENT)
            ClearFrequencyService.sf_server      = self.create_semaphore(self.SEM_F_SERVER)
            ClearFrequencyService.sf_samples     = self.create_semaphore(self.SEM_F_SAMPLES)
            ClearFrequencyService.sf_init        = self.create_semaphore(self.SEM_F_INIT)
            ClearFrequencyService.sf_clrfreq     = self.create_semaphore(self.SEM_F_CLRFREQ)
            ClearFrequencyService.sf_processed   = self.create_semaphore(self.SEM_F_PROCESSED)
            ClearFrequencyService.sl_samples     = self.create_semaphore(self.SEM_L_SAMPLES)
            ClearFrequencyService.sl_init        = self.create_semaphore(self.SEM_L_INIT)
            ClearFrequencyService.sl_clrfreq     = self.create_semaphore(self.SEM_L_CLRFREQ)
            ClearFrequencyService.semaphores = [
                ClearFrequencyService.sf_client,
                ClearFrequencyService.sf_server,
                ClearFrequencyService.sf_samples,
                ClearFrequencyService.sf_init,
                ClearFrequencyService.sf_clrfreq,
                ClearFrequencyService.sf_processed,
                ClearFrequencyService.sl_samples,
                ClearFrequencyService.sl_init,
                ClearFrequencyService.sl_clrfreq,
            ]
            ClearFrequencyService.shm_objects = [
                self.create_shm_obj(self.SAMPLES_SHM_NAME ,         self.SAMPLES_SHM_SIZE       , self.SAMPLES_ELEM_NUM), 
                self.create_shm_obj(self.FCENTER_SHM_NAME,          self.FCENTER_SHM_SIZE       , ),
                self.create_shm_obj(self.CLR_RANGE_SHM_NAME,        self.CLR_RANGE_SHM_SIZE     , self.CLR_RANGE_ELEM_NUM), 
                self.create_shm_obj(self.BEAM_NUM_SHM_NAME,         self.BEAM_NUM_SHM_SIZE      , ), 
                self.create_shm_obj(self.SAMPLE_SEP_SHM_NAME,       self.SAMPLE_SEP_SHM_SIZE    , ),
                self.create_shm_obj(self.RESTRICT_SHM_NAME,         self.RESTRICT_SHM_SIZE      , self.RESTRICT_ELEM_NUM), 
                self.create_shm_obj(self.META_DATA_SHM_NAME,        self.META_DATA_SHM_SIZE     , self.META_ELEM_NUM),
                self.create_shm_obj(self.ANTENNA_SHM_NAME,          self.ANTENNA_SHM_SIZE       , ),
                self.create_shm_obj(self.CLRFREQ_SHM_NAME,          self.CLR_BANDS_SHM_SIZE     , self.CLR_BANDS_ELEM_NUM), 
                self.create_shm_obj(self.SITE_ID_SHM_NAME,          self.SITE_ID_SHM_SIZE       , self.SITE_ID_ELEM_NUM),
                self.create_shm_obj(self.RADAR_ID_SHM_NAME,         self.RADAR_ID_SHM_SIZE      , ),
                self.create_shm_obj(self.ACTIVE_CLIENTS_SHM_NAME,   self.ACTIVE_CLIENTS_SHM_SIZE, )
            ]

            for obj in ClearFrequencyService.shm_objects:
                obj['shm_fd'] = self.initialize_shared_memory(obj['name'])
                            
            ClearFrequencyService.active_clients_fd = None 
            self.initialize_active_clients_counter()
            print("[clearFrequencyService] Done Initializing...\n\n")

        except ValueError:
            print("[ClearFrequencyService] Initialization Failed. Cleaning up SHM Objects and Semaphores...")
            ClearFrequencyService.soft_kill = True
            self.cleanup_shm()
        except KeyboardInterrupt:
            print("[CFS] Keyboard Interupt triggered during Initialization... Canceling and cleaning up...")
            ClearFrequencyService.soft_kill = True
            self.cleanup_shm()
            
        
    def create_shm_obj(self, name: str, size: int, elem_num= 1):
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
                    - 'elem_num' or number of elements the object contains
        """
        return {
            'name': name,
            'shm_ptr': None,
            'shm_fd': -1,
            'size': size,
            'elem_num': elem_num
        }
        
    def create_semaphore(self, name: str):
        return {
            'name': name,
            'sem':  self.initialize_semaphore(name)
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

    def initialize_semaphore(self, name):
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
                return semaphore
            except posix_ipc.ExistentialError:
                print("[clearFrequencyService] Semaphore not found. Retrying...")
            attempts += 1
            time.sleep(self.RETRY_DELAY)
        print(f"[clearFrequencyService] Failed to initialize Semaphore {name} after multiple attempts. Exiting.")
        exit(1)

    @classmethod
    def initialize_active_clients_counter(self):
        attempts = 0
        while attempts < self.RETRY_ATTEMPTS:
            try:
                # Init counter
                print(f"[clearFrequencyService] Attempting to initialize Active Clients Counter (Attempt {attempts + 1}/{self.RETRY_ATTEMPTS})...")
                self.active_clients_fd = os.open(f"/dev/shm", os.O_RDWR | os.O_TMPFILE, 0o666)
                os.ftruncate(self.active_clients_fd, struct.calcsize('i'))  # Ensure the size of the shared memory object is large enough for an integer
                # If abnormal num of clients, set to 0
                with mmap.mmap(self.active_clients_fd, struct.calcsize('i'), mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) as m:
                    m.seek(0)
                    current_value = struct.unpack('i', m.read(struct.calcsize('i')))[0]
                    if current_value < 0:  # Arbitrary threshold to detect abnormal num of clients
                        m.seek(0)
                        m.write(struct.pack('i', 0))
                    m.seek(0)
                    current_value = struct.unpack('i', m.read(struct.calcsize('i')))[0]
                print("[clearFrequencyService] Created Active Clients Counter... @ ", current_value)
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
            print(f"[clearFrequencyService] Decremented Active Clients Counter: {active_clients}\n")
            return active_clients
      
              
    def detect_dtype(self, var):
        """Recursively detects whether the variable is an int (i) or float (d). 
        Note that float is considered 'd' so that it can be used for struct.pack()

        Args:
            var (any): Variable to be detected

        Raises:
            ValueError: If the variable doesn't contain either ints or floats.

        Returns:
            _type_: Either 'd' or 'i' for float or int respectively.
        """
        # Continuation
        if isinstance(var, (list, set)):
            return self.detect_dtype(var[0]) 
        
        # Break Condition and Break
        elif isinstance(var, float):
            return 'd'
        elif isinstance(var, int):
            return 'i'
        else:
            raise ValueError(f"var ({var}) is contains neither float nor integers.")
                   
    def find_list_of_lists(self, var):
        """Recursively finds the point in the variable where the following 
        is true:
        (list_of_list -> list -> elem)
        This point is hereby called list_of_lists (LoL) for simplicity and is used 
        to flatten these arrays where applicable. 
        
        Note this function only works for variables that contain either ints or floats!

        Args:
            var (any): Variable to be parsed for LoL.

        Raises:
            ValueError: If the variable doesn't contain either ints or floats.

        Returns:
            _type_: Can be a list (the list_of_lists), 'd' or 'i' if var is a singular var.
        """
        # Continuation
        if isinstance(var, (list, set)):
            result = self.find_list_of_lists(var[0])
            
            # Break Conditions (return list_of_list->list->elem)
            if result  == 'elem':
                return 'list'
            elif result == 'list':
                return var
        
        # Break
        elif isinstance(var, (int, float)):
            return 'elem'
        else:
            raise ValueError(f"var ({var}) is contains neither float nor integers.")
    
    def write_data(self, obj, array_data, atype=''):
        """Writes data from array_data onto the object's shared memory pointer. 
        "Sends data from array_data across the obj's channel"

        Args:
            obj (dict): Object Dictionary contain Shared Memory data for the object. 
            array_data (list): List of data points 
            complex (bool, optional): Flag to write and unpack array_data 
                from its complex notation. Defaults to False.

        Raises:
            ValueError: If the variable doesn't contain either ints or floats.
        """
        try:    
            # Debug: Record start time
            start_time = time.time()
            
            flattened_data = []
            if atype == 'complex':
                # Convert to np array
                array_data_np = np.array(array_data, dtype=np.complex64)
                
                # Flatten and interleave real and imaginary parts as integers
                interleaved_data = np.empty(array_data_np.size * 2, dtype=np.int32)
                interleaved_data[0::2] = array_data_np.real.astype(np.int32).ravel()
                interleaved_data[1::2] = array_data_np.imag.astype(np.int32).ravel()
                
                # Print set per 2500 elemnents (till 5 set) in interleaved_data to verify
                for i in range(0, interleaved_data.size // 5000):
                    print(f"[Frequency Client] interleaved_data: ", interleaved_data[i * 5000:(i + 1) * 5000], "...")
                
                # Write directly to shared memory
                obj['shm_ptr'].seek(0)
                obj['shm_ptr'].write(interleaved_data.tobytes())
                
                return
            elif atype == 'meta':
                for i in range (1, len(array_data)):
                    flattened_data.append(array_data[i])
                # Place antenna list last
                flattened_data += array_data[0]
            elif atype == "sid":
                for letter in array_data:
                    flattened_data.append(bytes(letter, 'ascii'))
            else:
                # Otherwise, just flatten                
                list_of_lists = self.find_list_of_lists(array_data)
                
                # Element/1D List Found
                if type(list_of_lists) is str:
                    flattened_data = array_data
                # Greater-than-1D list Found
                elif type(list_of_lists) is list: 
                    print("2D list detected! Flattening...")
                    for row in list_of_lists:
                        flattened_data += row
                # Fail: Unexpected value found
                else: 
                    raise ValueError(f"An unexpected value occured: {list_of_lists}")
            
            
            # Determine dtype for Packing
            print(f"flattened array type: {type(flattened_data)}")
            dtype = 'i'
            if atype == 'meta':
                dtype = 'd'
            elif atype == "sid":
                dtype = b'c'
            else: 
                dtype = self.detect_dtype(flattened_data)
            print(f"dtype: {dtype}, elem_num: {obj['elem_num']}, ")
                
            print(f"flattened array type: {type(flattened_data)}")
                
            # Pack and write data
            if type(flattened_data) is list or type(flattened_data) is str:  
                print("[Frequency Client] new_data len of: ", len(flattened_data))
                if atype == 'complex':
                    print("[Frequency Client] Writing data:\n", flattened_data[:1], "...")
                    
                else:
                    print("[Frequency Client] Writing data:\n", flattened_data)                
                    
                obj['shm_ptr'].seek(0)
                if atype == 'sid':
                    print(f"ascii bytes: {flattened_data}")
                    print(f"dtype argument: {dtype * obj['elem_num']}")
                #     obj['shm_ptr'].write(struct.pack(dtype * obj['elem_num'], bytes(flattened_data, 'ascii'))) 
                # else: 
                obj['shm_ptr'].write(struct.pack(dtype * obj['elem_num'], *flattened_data)) 
            else:
                print("[Frequency Client] new_data len of: ", 1)
                print("[Frequency Client] Writing data:\n", flattened_data)
                
                obj['shm_ptr'].seek(0)
                obj['shm_ptr'].write(struct.pack(dtype * 1, flattened_data))
        
        
        except AttributeError as e:
            # Display error if element size is incorrect
            print("[Frequency Client] ERROR: Element Size is incorrect. send()'s parameters were likely not assigned properly. Please verify...")
            print(f"AttributeError: {e}")
            print(f"Object: {obj}, Attributes: {dir(obj)}")
            raise
        
        finally:
            # Debug: Print time to write
            end_time = time.time()
            elapsed_time = end_time - start_time
            
            # print(f"[Frequency Client] Time to write {obj['elem_num']} elements: {elapsed_time:.6f} seconds")
                    
    def read_m_data(self, obj):
        """Reads in data from the shared memory file descriptor.

        Args:
            obj (dict): Data object containing file descriptor, shared memory size, and number of elements expected.

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
                packed_data.append(int(((start_freq + end_freq) / 2) / 1000))
                noise_data.append(noise)
            return packed_data, noise_data 
        
    def premap_shm(self, meta_data=None):
        """Premaps all shared memory objects' pointers to their memory addresses.  
        """
        # If no SHM mapping and meta_data exist, map all SHM objects
        if self.shm_objects[0]['shm_ptr'] == None and meta_data != None:
            
            ## Check for Premapped antenna num
            # Map shared memory object pointer for antenna num
            print(f"Mapping {self.shm_objects[7]['name']}")
            self.shm_objects[7]['shm_ptr'] = mmap.mmap(self.shm_objects[7]['shm_fd'], self.shm_objects[7]['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE)
            
            # Check if Antenna Num changed, update corresponding values before they're mapped
            shm_ant_num = self.read_m_data(self.shm_objects[7])[0]
            print("SHM Antenna_num:  ", shm_ant_num)
            print("Meta Antenna num: ", len(meta_data['antenna_list']))
            if shm_ant_num != self.cur_antenna_num or self.cur_antenna_num != len(meta_data['antenna_list']) or self.shm_objects[0]['elem_num'] != (len(meta_data['antenna_list']) * int(meta_data['number_of_samples']) * 2):
                print("Antenna_num has been changed, updating SHM values before further SHM mapping...")
                self.cur_antenna_num = len(meta_data['antenna_list'])
                
                # Update meta SHM values
                meta_obj = self.shm_objects[6]
                meta_obj['elem_num'] = len(meta_data['antenna_list']) + self.META_ELEM
                meta_obj['size'] = meta_obj['elem_num'] * self.DOUBLE_SIZE
                os.ftruncate(meta_obj['shm_fd'], meta_obj['size'])
                
                # Update samples SHM values
                samples_obj = self.shm_objects[0]
                samples_obj['elem_num'] = len(meta_data['antenna_list']) * int(meta_data['number_of_samples']) * 2
                samples_obj['size'] = samples_obj['elem_num'] * self.INT_SIZE
                os.ftruncate(samples_obj['shm_fd'], samples_obj['size'])
                
            # Map shared memory object pointers
            print(f"Mapping Shared Memory for Objects...\n")
            for obj in self.shm_objects:
                # Special: Skip Antenna_Num mapping
                if obj['name'] == '/antenna_num':
                    continue
                print(f"Mapping {obj['name']}")
                obj['shm_ptr'] = mmap.mmap(obj['shm_fd'], obj['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE) 
    
    
    def send_samples(self, raw_samples, radar_id, fcenter=None, meta_data=None):
        """ Waits for client requests, then processes server data, writes client 
            data, and requests server to process new data. When process is 
            terminated, the try/finally block cleans up.
            
            Note: fcenter and meta_data can be None after being passed as arguments on the first send_samples() method call.
        """
        input_data = [
            raw_samples, 
            fcenter, 
        #     clr_range,
        ]
        
        meta_data_list = [
                        meta_data['antenna_list'],
                        meta_data['number_of_samples'],
                        meta_data['x_spacing'],
                        meta_data['usrp_rf_rate'],
                    ]
        
        # Special: Halt all future ClearFreqService
        if self.soft_kill is True:
            return
                
        # Fail: If no antennas, skip to recover on next Clear Search Cycle 
        if meta_data is None or len(meta_data['antenna_list']) == 0:
            print("[clearFrequencyService] ERROR: No antennas found. Skipping...")
            return
                
        # Get in Queue
        active_clients = self.increment_active_clients()
        print(f"[clearFrequencyService] Active clients count: {active_clients}")
        
        try:
            self.premap_shm(meta_data)
                                        
            # Await for a Client Request
            print("[clearFrequencyService] Awaiting Client Request...\n")
            self.sf_client['sem'].acquire()
            print("[clearFrequencyService] Acquired Client Request...")
            
            # Check & Send Initialization Data
            if meta_data is not None:
                print("[clearFrequencyService] Requesting Initialization Semaphore...")
                self.sl_init['sem'].acquire()
                print("[clearFrequencyService] Initialization Semaphore Acquired...")

                # If meta_data has changed
                if self.old_meta_data != meta_data_list:
                    self.old_meta_data = meta_data_list
                    shm_ant_num = self.read_m_data(self.shm_objects[7])
                    
                    # If antenna length or sample_num has changed, send, set, and sync with server
                    if self.cur_antenna_num != len(meta_data['antenna_list']) or self.shm_objects[0]['elem_num'] != (len(meta_data['antenna_list']) * int(meta_data['number_of_samples']) * 2):
                        print(f"[Frequency Client] Antenna_num changed. Reallocating memory")
                        self.cur_antenna_num = len(meta_data['antenna_list'])
                        
                        # Send Antenna Num
                        print(f"[Frequency Client] Data Write Progress: {self.shm_objects[7]['name']}")
                        self.write_data(self.shm_objects[7], len(meta_data['antenna_list']))
                                                
                        # Reallocate meta SHM
                        meta_obj = self.shm_objects[6]
                        meta_obj['elem_num'] = len(meta_data['antenna_list']) + self.META_ELEM
                        meta_obj['size'] = meta_obj['elem_num'] * self.DOUBLE_SIZE
                        os.ftruncate(meta_obj['shm_fd'], meta_obj['size'])
                        meta_obj['shm_ptr'] = mmap.mmap(meta_obj['shm_fd'], meta_obj['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE)
                        
                        # Reallocate samples SHM
                        samples_obj = self.shm_objects[0]
                        samples_obj['elem_num'] = len(meta_data['antenna_list']) * int(meta_data['number_of_samples']) * 2
                        samples_obj['size'] = samples_obj['elem_num'] * self.INT_SIZE
                        os.ftruncate(samples_obj['shm_fd'], samples_obj['size'])
                        samples_obj['shm_ptr'] = mmap.mmap(samples_obj['shm_fd'], samples_obj['size'], mmap.MAP_SHARED, mmap.PROT_READ | mmap.PROT_WRITE)
                    
                    # If server's antenna num is outdated, update it
                    elif shm_ant_num != self.cur_antenna_num:
                        # Send
                        print(f"[Frequency Client] Data Write Progress: {self.shm_objects[7]['name']}")
                        self.write_data(self.shm_objects[7], len(meta_data['antenna_list']))
                        
                    print(f"[Frequency Client] Data Write Progress: {self.shm_objects[6]['name']}")
                    
                    # Rearrange meta_data ordering
                    self.write_data(self.shm_objects[6], meta_data_list, 'meta')
                    
                # Write Site ID (SID)
                print(f"[Frequency Client] Data Write Progress: {self.shm_objects[9]['name']}")
                print(f"    len of objects list is {len(self.shm_objects)}")
                self.write_data(self.shm_objects[9], self.sid, 'sid')
    
                self.sl_init['sem'].release()
                self.sf_init['sem'].release()
                print("[clearFrequencyService] Initialization Semaphore Released ...")
                print("[clearFrequencyService] Server Initialization Flag raised ...")

                                
            if raw_samples is not None:
                print("[clearFrequencyService] Awaiting Sample Semphore Lock...")
                self.sl_samples['sem'].acquire()

                # Write Sample data
                self.write_data(self.shm_objects[0], raw_samples, 'complex')
                                
                # If Sample-relevant Data given, write it                
                for i in range(1, self.SAMPLE_PARAM_NUM):
                    print(f"[Frequency Client] Data Write Progress: {i}/{self.SAMPLE_PARAM_NUM} {self.shm_objects[i]['name']}")
                    
                    # General: Write updated input data 
                    if input_data[i] is not None:
                        self.write_data(self.shm_objects[i], input_data[i])
                        
                # Write Radar ID
                self.write_data(self.shm_objects[10], radar_id)
                
                self.sl_samples['sem'].release()
                self.sf_samples['sem'].release()
                print("[Frequency Client] Done writing data to Shared Memory...")
                
                # Request Server 
                print("[clearFrequencyService] Requesting Server Response...")
                self.sf_server['sem'].release()
                
                        
        except KeyboardInterrupt:
            print("[clearFrequencyService] Keyboard interrupt received. Exiting...")
        except posix_ipc.ExistentialError or ValueError or AttributeError:
                print("[clearFrequencyService] Shared memory has been delinked. Exiting...")
        finally:
            active_clients = self.decrement_active_clients()
                
        return 
                
    def request_clr_freq(self, radar_id, beam_num=None, sample_sep=None, clr_range=None, ):
        """ Waits for client requests, then processes server data, writes client 
            data, and requests server to process new data. When process is 
            terminated, the try/finally block cleans up.\
                
            Note that sample_sep is not expected to change for each clr freq request.
        """
        
        input_data = [
            clr_range,
            beam_num, 
            sample_sep,
        ]
        
        # Special: Halt all future ClearFreqService
        if self.soft_kill is True:
            return
                
        # Get in Queue
        active_clients = self.increment_active_clients()
        print(f"[clearFrequencyService] Active clients count: {active_clients}")
        
        try:
            self.premap_shm()
            
            # Await for a Client Request
            print("[clearFrequencyService] Awaiting Client Request...\n")
            self.sf_client['sem'].acquire()
            print("[clearFrequencyService] Acquired Client Request...")
            
            # Write Input Data present
            print("[clearFrequencyService] Requesting ClrFreq Semaphore...")
            self.sl_clrfreq['sem'].acquire()
            print("[clearFrequencyService] ClrFreq Semaphore Acquired...")
            
            for i in range(self.SAMPLE_PARAM_NUM, self.SAMPLE_PARAM_NUM + 3):
                
                # If sample separation present, send and update, else skip
                if input_data[i - self.SAMPLE_PARAM_NUM] is not None and self.old_smsep != input_data[i - self.SAMPLE_PARAM_NUM]:
                    self.old_smsep = input_data[i - self.SAMPLE_PARAM_NUM]
                else: continue
                
                # Write present data
                if input_data[i - self.SAMPLE_PARAM_NUM] is not None:
                    print(f"[Frequency Client] Data Write: {self.shm_objects[i]['name']}") 
                    self.write_data(self.shm_objects[i], input_data[i - self.SAMPLE_PARAM_NUM])
                
            # Write Radar ID
            print(f"[Frequency Client] Data Write: {self.shm_objects[10]['name']}")
            self.write_data(self.shm_objects[10], radar_id)
                
            self.sl_clrfreq['sem'].release()
            print("[clearFrequencyService] ClrFreq Semaphore Released ...")
                                            
                                            
            # Send Clear Frequency and Server Request 
            print("[clearFrequencyService] Requesting Clear Freq...")
            self.sf_clrfreq['sem'].release()
            print("[clearFrequencyService] Requesting Server Response...")
            self.sf_server['sem'].release()
            
            
            # Read-in Clear Freq data
            print("[clearFrequencyService] Awaiting Server Response...")
            self.sf_processed['sem'].acquire()
            print("[clearFrequencyService] Recieved Server Response. Reading Clear Freq data...")
            self.sl_clrfreq['sem'].acquire()
            new_noise_data = []
            new_clrfreq_data = self.read_m_data(self.shm_objects[8])
            new_clrfreq_data, new_noise_data = self.repack_data(new_clrfreq_data, True)
            for clr_freq_and_noise in zip(new_clrfreq_data, new_noise_data):
                print(f"[clearFrequencyService] Clear Freq Band: | {clr_freq_and_noise[0]} (kHz), {clr_freq_and_noise[1]} (N/A) |")
            clr_freq, noise = new_clrfreq_data[0], new_noise_data[0]
            
            self.sl_clrfreq['sem'].release()
                    
        except KeyboardInterrupt:
            print("[clearFrequencyService] Keyboard interrupt received. Exiting...")
        except posix_ipc.ExistentialError or ValueError or AttributeError:
                print("[clearFrequencyService] Shared memory has been delinked. Exiting...")
        finally:
            active_clients = self.decrement_active_clients()
                
        return clr_freq, noise
    
    @classmethod
    def cleanup_shm(self):
        if self.soft_kill is False or self.CLEAN_ON_INACTIVE is False:
            print("[clearFrequencyService] No active clients remaining, but not cleaning up shared resources to keep service idle.")
            try:
                posix_ipc.unlink_shared_memory(self.ACTIVE_CLIENTS_SHM_NAME)
                print(f"Unlinked shared memory {self.ACTIVE_CLIENTS_SHM_NAME}")
            except posix_ipc.ExistentialError or ValueError or AttributeError:
                print(f"Shared memory {self.ACTIVE_CLIENTS_SHM_NAME} does not exist")

        else:
            print("[clearFrequencyService] No active clients remaining; cleaning up shared resources.")
            for obj in self.shm_objects:
                try:
                    posix_ipc.unlink_shared_memory(obj['name'])
                    print(f"Unlinked shared memory {obj['name']}")
                except posix_ipc.ExistentialError or ValueError or AttributeError:
                    print(f"Shared memory {obj['name']} does not exist")
            try:
                posix_ipc.unlink_shared_memory(self.ACTIVE_CLIENTS_SHM_NAME)
                print(f"Unlinked shared memory {self.ACTIVE_CLIENTS_SHM_NAME}")
            except posix_ipc.ExistentialError or ValueError or AttributeError:
                print(f"Shared memory {self.ACTIVE_CLIENTS_SHM_NAME} does not exist")

            for sem in self.semaphores:
                try:
                    posix_ipc.unlink_semaphore(sem['name'])
                    print(f"Unlinked semaphore {sem['name']}")
                except posix_ipc.ExistentialError:
                    print(f"Semaphore {sem['name']} does not exist")
                    
    def flag_debug(self, t1 = 0, t2 = 0, t3 = 0):
        
        # Await for a Client Request
        print("[clearFrequencyService] Awaiting Client Request...\n")
        self.sf_client['sem'].acquire()
        print("[clearFrequencyService] Acquired Client Request...")
                                        
        if t1 == 1: 
            self.sl_init['sem'].acquire()
            self.sl_init['sem'].release()
            
            self.sf_init['sem'].release()
            print("[clearFrequencyService] Processed init flags...\n")
            
        if t2 == 1:                             
            print("[clearFrequencyService] Awaiting Sample Semphore Lock...")
            self.sl_samples['sem'].acquire()
            self.sl_samples['sem'].release()
            
            self.sf_samples['sem'].release()
            print("[Frequency Client] Done writing data to Shared Memory...")
            
            # Request Server 
            print("[clearFrequencyService] Requesting Server Response...\n\n")
            self.sf_server['sem'].release()
        elif t3 == 1: 
            print("[clearFrequencyService] Requesting Sample Semaphore for beam num...")
            self.sl_samples['sem'].acquire()
            time.sleep(1)
            print("[clearFrequencyService] Sample Semaphore Acquired...")
            self.sl_samples['sem'].release()
            print("[clearFrequencyService] Sample Semaphore Released ...")
                                            
                                            
            # Send Clear Frequency Request
            print("[clearFrequencyService] Requesting Clear Freq...")
            self.sf_clrfreq['sem'].release()
            
            # Request Server 
            print("[clearFrequencyService] Requesting Server Response...")
            self.sf_server['sem'].release()
            
            
            # Read-in Clear Freq data
            print("[clearFrequencyService] Awaiting Server Response...")
            self.sf_clrfreq['sem'].acquire()
            self.sl_clrfreq['sem'].acquire()
            print("[clearFrequencyService] Recieved Server Response. Reading Clear Freq data...\n\n")
            
            self.sl_clrfreq['sem'].release()

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

# Used to mimick USRP functionality
def read_restrict_file(restrict_file):
    print(f"Opening file {restrict_file} ...")
    
    restricted_frequencies = []
    with open(restrict_file, 'r') as f:
        for line in f:
            print(line)
    
            if line[0] == '#' or line[0] == 'd' or len(line) < 8:
                continue
            line = line.split(' ')
            restrict_start = int(line[0]) * 1e3 # convert kHz units in restrict to Hz
            restrict_end = int(line[1]) * 1e3 # convert kHz units in restrict to Hz
            restricted_frequencies.append([restrict_start, restrict_end])

    return restricted_frequencies; 

# if __name__ == "__main__":
    # main()




RESTRICT_FILE = "/home/df/Desktop/PSU-SuperDARN/Freq_Server/utils/misc_param/restrict.dat.inst"
CFS = ClearFrequencyService(sid='lab')

# raw_samples, meta_data = read_sample_pickle("/data/repos/Freq_Server/utils/pickle_input/clrfreq_dump.1.pickle")
raw_samples, meta_data = read_sample_pickle("/home/df/Desktop/PSU-SuperDARN/Freq_Server/utils/pickle_input/clrfreq_dump.1.pickle")
clear_freq_range = [ int(12 * pow(10,6)), int(12.5 * pow(10,6)) ]
# restrict_data=read_restrict_file(RESTRICT_FILE)

# Test: only first two antennas and 2000 samples
meta_ant_full = meta_data['antenna_list']
meta_ant_partial = [0,2] 
trimmed_samples = []
# for i in meta_ant_partial:
#     trimmed_samples.append(raw_samples[i][:2000])
trimmed_samples = raw_samples[:2]
 
# trimmed_samples = raw_samples[:2][:2000]

# print(f"samples: {raw_samples[:2][:10]}")

# print(f"raw_samples size: {len(raw_samples) * len(raw_samples[0])}")
# print(f"raw_samples shape: {len(raw_samples)} x {len(raw_samples[0])} x {2} (antenna_num x sample_num x complex)")
# print(f"trimmed_s size: {len(trimmed_samples) * len(trimmed_samples[0])}")
print(f"trimmed_samples shape: {len(trimmed_samples)} x {len(trimmed_samples[0])} x {2} (antenna_num x sample_num x complex)")

# CFS.flag_debug(trimmed_samples, 
#                 clr_range=clear_freq_range, 
#                 fcenter=12000,
#                 beam_num=1,
#                 sample_sep=340,
#                 meta_data=meta_data
#                 )

# Note: this is used to test the flag_debug function
# CFS.flag_debug(t1=1, t2=1)
# while (True):    CFS.flag_debug(t3=1)

def flatten_raw_into_int_bytes(arr):
    """Stack and flatten a 3D array into a 2D array.
    
    Args:
        arr (np.ndarray): Input 3D array.
        
    Returns:
        np.ndarray: Flattened 2D array.
    """
    return np.stack((arr.real, arr.imag), axis=-1).reshape(-1, 2)

# # Test stack flattening for raw samples
# stack_flat_samples_1 = np.stack((raw_samples.real, raw_samples.imag), axis=-1).reshape(-1, 2)
# # print(f"stack_flat_samples_1 shape: {stack_flat_samples_1.shape}")
# print(f"stack_flat_samples_1: {stack_flat_samples_1[:10]}")



# Test flatten speed
# timeit.timeit()

i = 0
while (i < 20):
    clear_freq_range = [int(12 * pow(10,6)), int(12.5 * pow(10,6))]
   
    meta_data['antenna_list'] = meta_ant_full
    meta_data['number_of_samples'] = 2500
    CFS.send_samples(
        raw_samples, 
        radar_id=0,
        fcenter=12000,
        meta_data=meta_data
    )

    # Test dynamic SHM reallocation due to antenna resizing
    if trimmed_samples is not None:
        meta_data['antenna_list'] = meta_ant_partial
        # meta_data['number_of_samples'] = 2000
        CFS.send_samples(
            trimmed_samples, 
            radar_id=1,
            fcenter=12000,
            meta_data=meta_data
        )
    
    i += 1

    for j in range(0, 3):
        clear_freq_range = [int(12 * pow(10,6)), int(12.5 * pow(10,6))]
        CFS.request_clr_freq(
            radar_id=0,
            beam_num=j,
            clr_range=clear_freq_range, 
            sample_sep=340,     # only necesary on first request or if changing
        )
        
        # Test: Ensure radars have separate clear frequency ranges
        clear_freq_range = [int(12 * pow(10,6)), int(12.25 * pow(10,6))]
        CFS.request_clr_freq(
            radar_id=1,
            beam_num=j,
            clr_range=clear_freq_range, 
            sample_sep=340,     # only necesary on first request or if changing
        )
        

    # break

        
    #     # break

    #     CFS.request_clr_freq(
    #         beam_num=1,
    #     )
