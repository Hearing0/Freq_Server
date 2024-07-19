import matplotlib
# matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the spectrum data from CSV file
data = pd.read_csv('spectrum_output.csv')

# Print the DataFrame to check its content
print(data.head())

# # Ensure the expected columns exist
if 'Frequency' not in data.columns or 'Power' not in data.columns:
    print(f"Available columns: {data.columns}")
    raise ValueError("The required columns 'Frequency' and 'Power' are not present in the CSV file.")

# Plot the spectrum data (Freq in MHz)
# data['Power'] = data['Power'] 
Power = data['Power'].to_numpy() 
Freq = data['Frequency'].to_numpy() / 1e6


plt.plot(Freq, Power)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Power')
plt.title('Spectrum Analysis')
plt.grid(True)
plt.ylim(0,10000)
# plt.xticks(np.arange(min(Freq), max(Freq) + 1e3, 1e6))
plt.yticks(np.arange(0, 11e3, 2000))
# TODO: scale freq to 1e6 instead


# Display the plot
plt.savefig("plot.png")