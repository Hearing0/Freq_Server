import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the spectrum data from CSV file
data = pd.read_csv('utils/csv_dump/spectrum_output.csv')
# data = pd.read_csv('../Freq_Server/utils/csv_dump/avg_sample.csv')
# data = pd.read_csv('utils/csv_dump/spectrum_output.avg1.csv')
# data = pd.read_csv('utils/csv_dump/spectrum_output.NoAvg.csv')



# Print the DataFrame to check its content
print(data.head())

# Ensure the expected columns exist
if 'Frequency' not in data.columns or 'Power' not in data.columns:
    print(f"Available columns: {data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the CSV file.")

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

# Display the plot
# plt.savefig("plots/debug/avg_sample_plot.png")
plt.savefig("plots/debug/avg_sample_plot.Sequential.1.png")
# plt.savefig("plots/debug/spectrum_plot.NoAvg.png")
# plt.savefig("plots/debug/avg_spectrum_plot.2.png")
