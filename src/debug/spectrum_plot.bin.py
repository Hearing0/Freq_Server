import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Create list of bin files in the specified directory matching the pattern
import glob
bin_files = glob.glob('log/fft_spectrum/*.bin')

# Sort the files by modification time to get the latest one
bin_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
if not bin_files:
    raise FileNotFoundError("No .bin files found in the specified directory.")
latest_bin_file = bin_files[0]

# Load the latest spectral magnitude data from .bin file
with open(latest_bin_file, 'rb') as file:
    # Print the name of the file being processed
    print(f"Processing file: {latest_bin_file}")
    
    # Read number of samples
    num_samples = int.from_bytes(file.read(4), 'little')
    
    # Read frequency vector
    freq_vector = np.fromfile(file, dtype=np.float32, count=num_samples)
    
    # Read spectrum power data
    power_data = np.fromfile(file, dtype=np.float32, count=num_samples)
    
    # Create a DataFrame from the frequency and power data
    data = pd.DataFrame({
        'Frequency': freq_vector,
        'Power': power_data
    })
    

# Ensure the expected columns exist
if 'Frequency' not in data.columns or 'Power' not in data.columns:
    print(f"Available columns: {data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the file.")

# Plot the spectrum data (Freq in MHz)
Power = data['Power'].to_numpy() 
Freq = data['Frequency'].to_numpy() / 1e6


plt.plot(Freq, Power)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Power')
plt.title('Spectrum Analysis')
plt.grid(True)
plt.ylim(0,10000)
plt.yticks(np.arange(0, 11e3, 2000))

# Print completion message
print("Plotting completed successfully.")

# Save the plot
plt.savefig("plots/debug/avg_sample_plot.bin.png")
