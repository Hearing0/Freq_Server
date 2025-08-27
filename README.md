# Clear Frequency Server (CFS)

This is an isolated testbed for the Clear Frequency Service routine, which was designed for the US SuperDARN Group. It specializes in processing time-domain sample set from a Python script, identifying the lowest noise frequencies in the sample set, then returning the lowest noise frequency (clear freq/band) to the Python Client. It can process multiple sample sends and clear frequency requests for multi-radar and -channel sites. 

Additional features include:
- Average Antenna Power Tracker (good for determining problematic antennas)
- Temporal Clear Search (increases accuracy after ~1 min of runtime, i.e. 20x 3-sec sample sets)
- Spectrogram Plotter for TCS FFT Spectra data
- CFS Logger
- FFT Spectrum logs
- Clear Freq logs
