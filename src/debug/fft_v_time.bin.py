import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import struct

# Note that this ignores non-TCS files as they do not accurately depict the spectrum over time  

# Config Vars
FULL_SPECTRUM = False
N_FFT_PLOTTED = 12

# Directory Constants
FFT_DIRECTORY_PATH = 'log/fft_spectrum/*.bin'
CLR_DIRECTORY_PATH = 'log/clr_freq/*.bin'

INT_SIZE = 4
DOUBLE_SIZE = 8

# Create list of bin files in the specified directory matching the pattern
import glob
import re
fft_files = glob.glob(FFT_DIRECTORY_PATH)
clr_files = glob.glob(CLR_DIRECTORY_PATH)

# Sort the files by modification time from latest to oldest
fft_files.sort(reverse=True)
clr_files.sort(reverse=True)
if not fft_files:
    raise FileNotFoundError(f"No fft .bin files found in the specified directory: {FFT_DIRECTORY_PATH}")
if not clr_files:
    raise FileNotFoundError(f"No clr freq .bin files found in the specified directory: {CLR_DIRECTORY_PATH}")


# Parse N-hours of data to plot
fft_data = []
clr_freq_data = []
for latest_fft_file, latest_clr_file in fft_files, clr_files:

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
        
        # Pack data
        fft_data.append(pd.DataFrame({
            'Frequency': freq_vector,
            'Power': power_data
        }))
        

# Ensure the expected columns exist
if 'Frequency' not in fft_data.columns or 'Power' not in fft_data.columns:
    print(f"Available columns: {fft_data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the file.")

# Plot the spectrum data (Freq in MHz)
Power = fft_data['Power'].to_numpy()
Freq = fft_data['Frequency'].to_numpy() / 1e6

# Plot spectrum data
plt.figure(figsize=(16, 12))  # Increase the figure size
if FULL_SPECTRUM:
    plt.xlim(Freq[0], Freq[-1])
    # plt.xticks(np.arange(Freq[0], Freq[-1], .02))
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
