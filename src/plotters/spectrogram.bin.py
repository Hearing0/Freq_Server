import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import struct
from datetime import datetime, timezone, timedelta

# Note that this ignores non-TCS files as they do not accurately depict the spectrum over time  

# Config Vars
FULL_SPECTRUM = False
N_FFT_PLOTTED = 12

# Directory Constants
FFT_DIRECTORY_PATH = 'log/fft_spectrum/*.tcs.bin'
CLR_DIRECTORY_PATH = 'log/clr_freq/*.bin'
SAVE_PLOT_AS       = 'plots/debug/spectrogram.kod.png'

INT_SIZE = 4
DOUBLE_SIZE = 8
TIME_SIZE = 8

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

print(fft_files[0])


# Parse N-hours of data to plot
num_samples = 5000
set_count = 0

clr_freq_data = []
data = [[] for _ in range(3)]

print(data)

# Load the latest spectral magnitude data from .bin file
with open(fft_files[0], 'rb') as file:
    # num_samples_bytes = file.read(INT_SIZE)
    # if len(num_samples_bytes) < INT_SIZE:
    #     raise ValueError(f"File {fft_files[0]} is too short to contain num_samples int header.")
    # num_samples = struct.unpack('i', num_samples_bytes)[0]
    # print(num_samples)

    while True:
        # Store Time
        time_bytes = file.read(TIME_SIZE)
        if not time_bytes or len(time_bytes) < TIME_SIZE:
            break    
        timestamp = datetime.fromtimestamp(struct.unpack("<Q", time_bytes)[0], tz = timezone.utc)
        data[0].append(timestamp)
        time_start = timestamp

        # Store Frequency
        freq_bytes = file.read(DOUBLE_SIZE * num_samples)
        if not freq_bytes or len(freq_bytes) < DOUBLE_SIZE * num_samples:
            break
        freq_vector = struct.unpack('d' * num_samples, freq_bytes)

        # Store Power 
        power_bytes = file.read(DOUBLE_SIZE * num_samples)
        if not power_bytes or len(power_bytes) < DOUBLE_SIZE * num_samples:
            break
        power_data = struct.unpack('d' * num_samples, power_bytes)
        data[2].append(np.array(power_data))

        # Print starting data
        if set_count == 0: 
            print(f"Timestamp: {timestamp}")
            print(freq_vector[0])
            print(power_data[0])

        set_count += 1

    # Convert Freq from Hz to MHz
    data[1].append(np.array(freq_vector) / 1e6)

print(f"# of sets: {set_count}")


# Reformat data for simple plotting
print("Reformatting data for plotting...")
power   = np.arange(set_count * num_samples).reshape(num_samples, set_count)
for i in range(0,num_samples):
    for j in range(0, set_count):
        power[i][j] = data[2][j][i]

time    = np.array(data[0])
freq    = np.array(data[1][0])

print("shape: ", time.shape)
print("shape: ", freq.shape)
print("shape: ", power.shape)
print("Checking for NaN in Power...")
print("NaNs in power: ", np.isnan(power).any())

if len(time) == power.shape[0]: print("time axis alligned")
if len(freq) == power.shape[1]: print("freq axis alligned")


# Plot spectrum data
plt.figure(figsize=(16, 14))  # Increase the figure size
# if FULL_SPECTRUM:
#     plt.xlim(data[1][0], data[1][-1])
#     plt.xticks(np.arange(Freq[0], Freq[-1], .02))
print("Creating Temporal Spectrum Analysis Plot...")
plt.pcolormesh(
    time,
    freq,
    power,
    shading='auto',
    cmap='viridis',
    vmin=np.percentile(power, 5),
    vmax=np.percentile(power, 95)
)
cbar = plt.colorbar(label='Power (N/A)')
# cbar.ax.yaxis.set_major_formatter(plt.ticker.LogFormatter(labelOnlyBase=False))

# Format UTC time x-axis
ax = plt.gca()
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%H:%M:%S"))
plt.xticks(rotation=45)
plt.xlabel('Time (UTC)')

# Format Freq y-axis
plt.ylabel('Frequency (MHz)')
plt.suptitle(f'{time[0].date()}', y=.023, fontsize=10)
plt.title('Radar Spectrogram', y=1.01)
plt.grid(True)

# Display plot
plt.tight_layout()
plt.savefig(SAVE_PLOT_AS)


print("Plotting Finished")
print("Saved as:", SAVE_PLOT_AS)