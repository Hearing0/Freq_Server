
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import os

""" Creates a spectrum plot (Freqency vs Power) with frequency in MHz
    trimmed to only display the Clear Frequency Bands. 
    Uses .csv files for both the spectrum power and frequency steps for each 
    power sample (data), and for the clear frequency band and its min and max
    frequencies (clr_freq_data).

Raises:
    ValueError: "The required columns 'Frequency' and 'Power' are not present in the CSV file."
    ValueError: "The required columns 'Start Frequency', 'End Frequency', and 'Noise' are not present in the 'clr_freq_data' CSV file."
"""

fft_files = glob.glob('log/fft_spectrum/*.csv')
clr_log_files = glob.glob('log/clr_freq/clr_freq.*.csv')

# Sort the files by modification time to get the latest one
fft_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
clr_log_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
if not fft_files and clr_log_files:
    raise FileNotFoundError("No .bin files found in the specified directory.")

# debug: print the sorted file lists
print("Sorted FFT files:", fft_files)
print("Sorted Clear Frequency Band files:", clr_log_files)

latest_fft_file = fft_files[5]
latest_clr_log_file = clr_log_files[5] 

# Load the latest spectral magnitude data from csv file
print(f"Processing file: {latest_fft_file}")
spectra_data = pd.read_csv (latest_fft_file)

# Load the latest clear frequency band data from file
print(f"Processing file: {latest_clr_log_file}")
clr_freq_data = pd.read_csv(latest_clr_log_file)

    
    
# Print the DataFrame to check its content
print(spectra_data.head())
print(clr_freq_data.head(6))

# Ensure the expected columns exist
if 'Frequency' not in spectra_data.columns or 'Power' not in spectra_data.columns:
    print(f"Available columns: {spectra_data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the 'data' CSV file.")

if 'Start Frequency' not in clr_freq_data.columns or 'End Frequency' not in clr_freq_data.columns or 'Noise' not in clr_freq_data.columns:
    print(f"Available columns: {clr_freq_data.columns}")
    raise ValueError("The required columns 'Start Frequency', 'End Frequency', and 'Noise' are not present in the 'clr_freq_data' CSV file.")

if 'Clear Freq Start' not in clr_freq_data.columns or 'Clear Freq End' not in clr_freq_data.columns:
    print(f"Optional Columns 'Clear Freq Start' and 'Clear Freq End' not present.")
    print(f"Clear Band Labels WILL overlap!!!")

# Convert data to numpy arrays
Power = spectra_data['Power'].to_numpy() 
Freq = spectra_data['Frequency'].to_numpy() / 1e6  # Convert Frequency to MHz
# if 'Clear Freq Start' in clr_freq_data.columns and 'Clear Freq End' in clr_freq_data.columns:
#     clear_start = clr_freq_data['Clear Freq Start']
#     clear_end = clr_freq_data['Clear Freq End']
#     print("Clear |{clear_start} -- {clear_end}|")

# Plot spectrum data
plt.figure(figsize=(16, 12))  # Increase the figure size
if 'Clear Freq Start' in clr_freq_data.columns and 'Clear Freq End' in clr_freq_data.columns:
    plt.xlim(clr_freq_data['Clear Freq Start'][0] / 1e6, clr_freq_data['Clear Freq End'][0] / 1e6)
    plt.xticks(np.arange(clr_freq_data['Clear Freq Start'][0] / 1e6, clr_freq_data['Clear Freq End'][0] / 1e6, .02))
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
count = 0
for index, row in clr_freq_data.iterrows():
    start_freq = row['Start Frequency'] / 1e6  # Convert to MHz
    end_freq = row['End Frequency'] / 1e6  # Convert to MHz
    noise = row['Noise']
    
    # Ignore Dummy Bands 
    if noise < 100e3:
        # Mark the frequency band
        plt.axvspan(start_freq, end_freq, color=colors[index], alpha=.5)        
        mid_freq = (start_freq + end_freq) / 2
        plt.text(mid_freq, 20000 + 1500 * (index % 2), f"{round(start_freq, 2)}-{round(end_freq, 2)} MHz\nNoise: {round(noise, 2)}", 
                horizontalalignment='center', verticalalignment='bottom',fontsize=9, bbox=dict(facecolor='white', alpha=0.8))
        

# Display plot
plt.legend()
plt.savefig("plots/debug/spectrum_plot.clrbands.png")
