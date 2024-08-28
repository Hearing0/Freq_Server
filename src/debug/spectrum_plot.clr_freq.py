
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

""" Creates a spectrum plot (Freqency vs Power) with frequency in MHz
    trimmed to only display the Clear Frequency Bands. 
    Uses .csv files for both the spectrum power and frequency steps for each 
    power sample (data), and for the clear frequency band and its min and max
    frequencies (clr_freq_data).

Raises:
    ValueError: "The required columns 'Frequency' and 'Power' are not present in the CSV file."
    ValueError: "The required columns 'Start Frequency', 'End Frequency', and 'Noise' are not present in the 'clr_freq_data' CSV file."
"""

# Load the spectrum data from CSV files
data = pd.read_csv('utils/csv_dump/spectrum_output.csv')
clr_freq_data = pd.read_csv('utils/csv_dump/clr_freq_output.csv')

# Print the DataFrame to check its content
print(data.head())
print(clr_freq_data.head(6))

# Ensure the expected columns exist
if 'Frequency' not in data.columns or 'Power' not in data.columns:
    print(f"Available columns: {data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the 'data' CSV file.")

if 'Start Frequency' not in clr_freq_data.columns or 'End Frequency' not in clr_freq_data.columns or 'Noise' not in clr_freq_data.columns:
    print(f"Available columns: {clr_freq_data.columns}")
    raise ValueError("The required columns 'Start Frequency', 'End Frequency', and 'Noise' are not present in the 'clr_freq_data' CSV file.")

if 'Clear Freq Start' not in clr_freq_data.columns or 'Clear Freq End' not in clr_freq_data.columns:
    print(f"Optional Columns 'Clear Freq Start' and 'Clear Freq End' not present.")
    print(f"Clear Band Labels WILL overlap!!!")

# Convert data to numpy arrays
Power = data['Power'].to_numpy() 
Freq = data['Frequency'].to_numpy() / 1e6  # Convert Frequency to MHz
if 'Clear Freq Start' in clr_freq_data.columns and 'Clear Freq End' in clr_freq_data.columns:
    Clear_Start = clr_freq_data['Clear Freq Start']
    Clear_End = clr_freq_data['Clear Freq End']

# Plot spectrum data
plt.figure(figsize=(16, 12))  # Increase the figure size
# if 'Clear Freq Start' in clr_freq_data.columns and 'Clear Freq End' in clr_freq_data.columns:
#     plt.xlim(clr_freq_data['Clear Freq Start'][0], clr_freq_data['Clear Freq End'][0])
plt.plot(Freq, Power, label='Spectrum')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Power')
plt.title('Spectrum Analysis')
plt.grid(True)
plt.ylim(0, 10000)
plt.yticks(np.arange(0, 11e3, 2000))

# Create Color Palette 
colors = sns.color_palette("colorblind", len(clr_freq_data))

# Plot the clear frequency bands
for index, row in clr_freq_data.iterrows():
    start_freq = row['Start Frequency'] / 1e6  # Convert to MHz
    end_freq = row['End Frequency'] / 1e6  # Convert to MHz
    noise = row['Noise']
    
    # Ignore Dummy Bands 
    if noise < 100e3:
        # Mark the frequency band
        plt.axvspan(start_freq, end_freq, color=colors[index], alpha=.5)

        # Annotate start frequency, end frequency, and noise
        if index % 2:
            v_alignment = 'bottom'
        else:
            v_alignment = 'top'
        mid_freq = (start_freq + end_freq) / 2
        plt.text(mid_freq, 10500, f"{round(start_freq, 2)}-{round(end_freq, 2)} MHz\nNoise: {round(noise, 2)}", 
                horizontalalignment='center', verticalalignment=v_alignment,fontsize=9, bbox=dict(facecolor='white', alpha=0.8))


# Display plot
plt.legend()
plt.savefig("plot_with_bands.png")
