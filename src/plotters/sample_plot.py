import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Run script
# python3 src/debug/sample_plot.py 



base_filepath = 'utils/csv_dump/samples/'
files = ['avg_sample.csv']
# 'sample_im_output.csv', 'sample_re_output.csv', 

def retrieve_samples(filename):
    data = pd.read_csv(filename)

    # Print the DataFrame to check its content
    print(data.head())

    # Ensure the expected columns exist
    if 'Samples' not in data.columns or 'Power' not in data.columns:
        print(f"Available columns: {data.columns}")
        raise ValueError("The required columns 'Samples' and 'Power' are not present in the CSV file.")

    # Plot the spectrum data 
    return data

# Initialize arrays for Power and Samples
Power = []
Samples = []

# Load the spectrum data from CSV file
for i in len(2):
    file_name = base_filepath + files[i]
    data = retrieve_samples(file_name)
    
    Power.append(data['Power'].to_numpy())
    
    # Ensure Samples are the same for both files
    if i == 0:
        Samples = data['Samples'].to_numpy()

# Convert Power to a 2D array where N is the number of samples
Power = np.array(Power)

print(Power)

plt.plot(Samples, Power[0], color='r', label='Imaginary')
plt.plot(Samples, Power[1], color='g', label='Real')
plt.xlabel('Samples')
plt.ylabel('Power')
plt.title('Spectrum Analysis')
plt.grid(True)
plt.ylim(0,200)
plt.yticks(np.arange(0, 400, 50))
plt.legend()

# Display the plot
plt.savefig("plots/test_Plot.png")
