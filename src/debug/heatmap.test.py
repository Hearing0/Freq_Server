import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# sample data
set_count = 20
num_samples = 625

time = pd.date_range("2023-01-01", periods=set_count, freq="3s", tz="UTC")

freq = np.linspace(9.5,14.5, num_samples)

power = np.random.rand(num_samples, set_count)


plt.figure(figsize=(12,8))
plt.scatter(4, freq[np.argmax(power) // num_samples], c='red')
# plt.pcolormesh(time, freq, power, shading='auto', cmap='viridis')
# plt.colorbar

# Format UTC time axis
ax = plt.gca()
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%H:%M:%S"))
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig("plots/debug/spectrum_plot.heattest.bin.png")

