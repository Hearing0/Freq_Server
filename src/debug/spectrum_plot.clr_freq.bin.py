import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import struct

# Config Flags
FULL_SPECTRUM = False

# Directory Constants
FFT_DIRECTORY_PATH = 'log/fft_spectrum/*.bin'
CLR_DIRECTORY_PATH = 'log/clr_freq/*.bin'

INT_SIZE = 4
DOUBLE_SIZE = 8

# Create list of bin files in the specified directory matching the pattern
import glob
fft_files       = glob.glob(FFT_DIRECTORY_PATH)
clr_files  = glob.glob(CLR_DIRECTORY_PATH)

# Sort the files by modification time to get the latest one
fft_files.sort(reverse=True)
clr_files.sort(reverse=True)
if not fft_files:
    raise FileNotFoundError(f"No fft .bin files found in the specified directory: {FFT_DIRECTORY_PATH}")
if not clr_files:
    raise FileNotFoundError(f"No clr freq .bin files found in the specified directory: {CLR_DIRECTORY_PATH}")
latest_fft_file = fft_files[0]
latest_clr_file = clr_files[0]


# Load the latest spectral magnitude data from .bin file
with open(latest_fft_file, 'rb') as file:
    # Print the name of the file being processed
    print(f"Processing file: {latest_fft_file}")
    
    # Read number of samples
    num_samples = int.from_bytes(file.read(INT_SIZE), 'little')

    print(num_samples)
    
    # Read frequency vector
    freq_vector = struct.unpack('d' * num_samples, file.read(DOUBLE_SIZE * num_samples))

    print(freq_vector[10])
    
    # Read spectrum power data
    power_data = struct.unpack('d' * num_samples, file.read(DOUBLE_SIZE * num_samples))
    
    # Create a DataFrame from the frequency and power data
    data = pd.DataFrame({
        'Frequency': freq_vector,
        'Power': power_data
    })
    
# Load the latest clear freq data from .bin file
clr_search_range = []
clr_freq_data = []
with open(latest_clr_file, 'rb') as clr_file:
    # Print the name of the file being processed
    print(f"Processing file: {latest_clr_file}")

    # Read Clear Search Range
    clr_search_range = struct.unpack('i' * 2, clr_file.read(INT_SIZE * 2))

    # Read freq data
    f_start_data  = []
    noise_data    = []
    f_end_data    = []
    for idx in range(0, 6):
        f_start = struct.unpack('i', clr_file.read(INT_SIZE))[0]
        f_start_data.append(f_start)
        noise = struct.unpack('d', clr_file.read(DOUBLE_SIZE))[0]
        noise_data.append(noise)
        f_end = struct.unpack('i', clr_file.read(INT_SIZE))[0]
        f_end_data.append(f_end)

    # freq_data = []
    # freq_data = struct.unpack('i' * num_samples,)

    # Parse data into freq data
    for start_freq, noise, end_freq in zip(f_start_data, noise_data, f_end_data):
        clr_freq_data.append([start_freq, noise, end_freq])
        print(f"| {start_freq} -- Noise: {noise} -- {end_freq} |")

# Ensure the expected columns exist
if 'Frequency' not in data.columns or 'Power' not in data.columns:
    print(f"Available columns: {data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the file.")

# Plot the spectrum data (Freq in MHz)
Power = data['Power'].to_numpy()
Freq = data['Frequency'].to_numpy() / 1e6

# Plot spectrum data
plt.figure(figsize=(16, 12))  # Increase the figure size
if FULL_SPECTRUM:
    plt.xlim(Freq[0], Freq[-1])
    # plt.xticks(np.arange(Freq[0], Freq[-1], .02))
else:
    margin = 0.05  # MHz, adjust as needed
    x_min = clr_search_range[0] / 1e6 - margin
    x_max = clr_search_range[1] / 1e6 + margin
    plt.xlim(x_min, x_max)
    plt.xticks(np.arange(x_min, x_max, .02))
plt.plot(Freq, Power, label='Spectrum')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Power')
plt.title('Spectrum Analysis')
plt.grid(True)
plt.ylim(0, 10000)
plt.yticks(np.arange(0, 25e3, 2000))


# Create Color Palette 
colors = sns.color_palette("colorblind", len(clr_freq_data))

# Plot the clear frequency bands
idx = 0
for clr_freq in clr_freq_data:
    start_freq = clr_freq[0] / 1e6  # Convert to MHz
    noise = clr_freq[1]
    end_freq = clr_freq[2] / 1e6  # Convert to MHz
    
    # Ignore Dummy Bands 
    if noise < 100000:
        # Mark the frequency band
        plt.axvspan(start_freq, end_freq, color=colors[idx], alpha=.5)        
        mid_freq = (start_freq + end_freq) / 2
        plt.text(mid_freq, 20000 + 1500 * (idx % 2), f"{round(start_freq, 2)}-{round(end_freq, 2)} MHz\nNoise: {round(noise, 2)}", 
                horizontalalignment='center', verticalalignment='bottom',fontsize=9, bbox=dict(facecolor='white', alpha=0.8))
        
    idx += 1


# Display plot
plt.legend()
plt.savefig("plots/debug/spectrum_plot.clrbands.bin.png")
