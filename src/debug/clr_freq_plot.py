import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configuration flags
USE_CSV = True  # Switch to False when implementing .bin support
FILE_EXTENSION = '.csv' if USE_CSV else '.bin'

# Path constants
BASE_PATH = pathlib.Path('/logs/fft_spectrum')

CSV_DIR = BASE_PATH / 'csv_files'
BIN_DIR = BASE_PATH / 'bin_files'
PLOT_DIR = BASE_PATH / ''

# Ensure directories exist
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# Data storage
freqs_list = []
noise_list = []


def retrieve_data(file_path, mode='csv'):
    """Retrieve data from file based on format"""
    
    if mode == 'csv':
        df = pd.read_csv(file_path)
        return df['Start Frequency'].to_numpy(), df['Noise'].to_numpy()
    elif mode == 'bin':
        # Implement .bin parsing later
        raise NotImplementedError("BIN format not implemented yet")
    else:
        raise ValueError(f"Unsupported mode: {mode}")
    
    
# Get appropriate directory based on file extension
source_dir = CSV_DIR if USE_CSV else BIN_DIR

# Process files
for file_path in source_dir.glob(f'*{FILE_EXTENSION}'):
    if file_path.is_file():
        try:
            freqs, noise = retrieve_data(file_path, mode=FILE_EXTENSION)
            freqs_list.extend(freqs)
            noise_list.extend(noise)
        except Exception as e:
            print(f"Error processing {file_path.name}: {str(e)}")

# Convert to numpy arrays
freqs_array = np.array(freqs_list)
noise_array = np.array(noise_list)

# Create plot
plt.figure(figsize=(10, 6))
plt.scatter(freqs_array, noise_array, s=10, alpha=0.5)
plt.title('Frequency vs Noise')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Noise')
plt.grid(True)

# Save plot
plot_path = PLOT_DIR / 'frequency_vs_noise.png'
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Plot saved to {plot_path}")