import numpy as np

def parse_array(arr, complex_data=False, parse_in_pairs=False):
    if complex_data:
        # If complex, flatten and separate real and imaginary parts
        flattened_array = np.column_stack((arr.real.flatten(), arr.imag.flatten())).flatten()
    else:
        # Otherwise, just flatten
        flattened_array = arr.flatten()
    
    if parse_in_pairs:
        # Parse in pairs
        return [(flattened_array[i], flattened_array[i+1]) for i in range(0, len(flattened_array), 2)]
    else:
        return flattened_array
    
def repack_data(read_data, clr_freq=False):
        """Repacks read data from SHM into its proper format. Currently repacks
        the following:
        - Clear Frequency

        Args:
            read_data (_type_): _description_
            data_size (_type_): _description_
            data_sub_size (int, optional): _description_. Defaults to 1.
            complex (bool, optional): _description_. Defaults to False.
        """
        packed_data = []
        if clr_freq:
            for start_freq, noise, end_freq in zip(read_data[::3], read_data[1::3], read_data[2::3]):
                # Return Center Freq and Noise
                packed_data.append([(start_freq + end_freq) / 2, noise])
                print(start_freq)
                print(end_freq)
                print(noise)

# Example usage
arr = np.array([[1000, 30000], [50000, 7000000], [900000 , 1100000]])
print(arr)
result = parse_array(arr, complex_data=False, parse_in_pairs=False)
print(result)
repacked_result = repack_data(result, clr_freq=True)
print(repacked_result)
