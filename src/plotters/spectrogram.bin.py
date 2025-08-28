import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import struct
from datetime import datetime, timezone, timedelta
import argparse
import glob
import re

### Note that this plotter ignores non-TCS files as they do not accurately depict the spectrum over time  
# Change station ID here or via command line argument
DEFAULT_STID = 'kod'

# Config Vars
FULL_SPECTRUM = False
N_FFT_PLOTTED = 12

# Directory Constants
FFT_DIRECTORY_PATH  = 'log/fft_spectrum/*.tcs.bin'
PLOT_DIRECTORY_PATH = 'plots/debug/'

INT_SIZE = 4
DOUBLE_SIZE = 8
TIME_SIZE = 8

# num_samples = 5000
channels = ['a', 'b', 'c', 'd']

# Parse the date string
def parse_date_argument(date_string):
    if date_string == 'yesterday' or date_string == 'y':
        return datetime.now(timezone.utc) - timedelta(days=1)
    try:
        # Example format: YYYY-MM-DD
        return datetime.strptime(date_string, '%Y-%m-%d')
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid date format: '{date_string}'. Expected format: YYYY-MM-DD"
        )
def parse_time_argument(time_string):
    try:
        # Example format: YYYY-MM-DD HH:MM:SS
        return datetime.strptime(time_string, '%H:%M:%S').time()
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid time format: '{time_string}'. Expected format: HH:MM:SS"
        )

parser = argparse.ArgumentParser(description='Generate spectrogram from FFT TSC binary files.')
parser.add_argument(
    '-stid',
    type=str,
    default=DEFAULT_STID,
    help='Specify the station ID. Default is "kod".'
)
parser.add_argument(
    # If no channel is specified, plot all channels found in the files
    '-ch',
    '--channel',
    type=str,
    help='Specify the channel to plot.'
)
parser.add_argument(
    '-d',
    '--date',
    type=parse_date_argument,
    help="Specify a UTC date in 'YYYY-MM-DD' format. Defaults to current UTC date."
)
parser.add_argument(
    '-ts',
    '--time_start',
    type=parse_time_argument,
    default=datetime.strptime('00:00:00', '%H:%M:%S').time(),
    help="Specify start time in 'HH:MM:SS' format. Default is 00:00:00."
)
parser.add_argument(
    '-te',
    '--time_end',
    type=parse_time_argument,
    default=datetime.strptime('23:59:59', '%H:%M:%S').time(),
    help="Specify end time in 'HH:MM:SS' format. Default is 23:59:59."
)
args = parser.parse_args()

# Default to current UTC date if no date argument is provided
if args.date:
    # Use the provided date argument
    date_str = args.date.strftime('%Y-%m-%d')
else:
    date_str = datetime.now(timezone.utc).strftime('%Y-%m-%d') 

print(f"Generating {args.stid} spectrogram for {date_str} from {args.time_start} to {args.time_end} UTC")

all_fft_files = glob.glob(FFT_DIRECTORY_PATH)

# If a channel is specified, only filter that channel
if args.channel:
    channels = [args.channel]
    print(f"Channel specified: {args.channel}")

# If no channel is specified, filter files and plot per channel
else:
    for ch in channels:
        print(f"Processing channel: {ch}")
        # Filter files based on the specified date
        date_pattern = re.compile(          # Pattern must match 20250828.14.kod.d.fft.tcs.bin where kod is the station ID and d is the current channel to plot
                fr'(\d{{4}})(\d{{2}})(\d{{2}})\.(\d{{2}})\.{args.stid}\.{ch}\.fft\.tcs\.bin$'
        )
        def file_matches_date(f, date_str):
            m = date_pattern.search(os.path.basename(f))
            if not m:
                return False
            file_date = f"{m.group(1)}-{m.group(2)}-{m.group(3)}"
            return file_date == date_str

        fft_files = [f for f in all_fft_files if file_matches_date(f, date_str)]

        # Also include the latest file from the previous day to capture early UTC data
        previous_date = (args.date - timedelta(days=1)) if args.date else (datetime.now(timezone.utc) - timedelta(days=1))
        previous_date_str = previous_date.strftime('%Y-%m-%d')
        print("Also checking previous day:", previous_date_str)
        previous_day_files = [f for f in all_fft_files if file_matches_date(f, previous_date_str)]
        print(f"Found {len(previous_day_files)} files from previous day")
        if previous_day_files:
            latest_previous_file = max(previous_day_files, key=os.path.getmtime)
            fft_files.append(latest_previous_file)

        # Sort files by name (which includes timestamp)
        fft_files.sort()
        print(f"Found {len(fft_files)} FFT files for date {date_str}")
        if len(fft_files) == 0:
            print("No FFT files found for the specified date. Checking next channel...")
            continue



        # Load the latest spectral magnitude data from .bin file
        set_count = 0
        found_time_before_start_bound = False
        found_time_end_bound = False
        clr_freq_data = []
        data = [[] for _ in range(3)]
        for fft_file in fft_files:
            print("Loading file:", os.path.basename(fft_file))

            # Check if file exists
            if not os.path.isfile(fft_file):
                print("File not found:", os.path.basename(fft_file))
                continue

            with open(fft_file, 'rb') as file:
                while True:
                    found_time_before_start_bound = False

                    # Read number of samples
                    num_samples_bytes = file.read(INT_SIZE)
                    if not num_samples_bytes or len(num_samples_bytes) < INT_SIZE:
                        break
                    num_samples = struct.unpack("<I", num_samples_bytes)[0]
                    print("Num Samples:", num_samples)


                    # Read beam number (not used in this plotter)
                    beam_bytes = file.read(INT_SIZE)
                    if not beam_bytes or len(beam_bytes) < INT_SIZE:
                        break
                    beam_number = struct.unpack("<I", beam_bytes)[0]
                    print("Beam Number:", beam_number)


                    # Read Time
                    time_bytes = file.read(TIME_SIZE)
                    if not time_bytes or len(time_bytes) < TIME_SIZE:
                        break
                    timestamp = datetime.fromtimestamp(struct.unpack("<Q", time_bytes)[0], tz=timezone.utc)
                    print("Timestamp:", timestamp)

                    # If date argument is provided, filter data by time range
                    if args.date:
                        # Skip data before start time
                        if timestamp.time() < args.time_start and timestamp.date() == args.date.date() or timestamp.date() < args.date.date():
                            found_time_before_start_bound = True
                        
                        # Stop processing if end time reached
                        elif timestamp.time() > args.time_end and timestamp.date() >= args.date.date():
                            print("FFT file out of time range: ", timestamp)
                            found_time_end_bound = True
                            break


                    # Read Frequency
                    # print("Reading Freq Vector...")
                    freq_bytes = file.read(DOUBLE_SIZE * num_samples)
                    if not freq_bytes or len(freq_bytes) < DOUBLE_SIZE * num_samples:
                        break


                    # Read Power 
                    power_bytes = file.read(DOUBLE_SIZE * num_samples)
                    if not power_bytes or len(power_bytes) < DOUBLE_SIZE * num_samples:
                        break
                    power_data = struct.unpack('d' * num_samples, power_bytes)

                    # Only store data within the specified time range            
                    if found_time_before_start_bound is False:
                        # Store valid data
                        data[0].append(timestamp)
                        f_time_start = timestamp
                        freq_vector = struct.unpack('d' * num_samples, freq_bytes)
                        data[2].append(np.array(power_data))
                        set_count += 1

                        # Print processed data
                        if set_count % 2500 == 0: 
                            print(f"Processed: {timestamp}")
                            # print(freq_vector[0])
                            # print(power_data[0])

                # Store Freq vector only once
                if set_count > 0:
                    data[1].append(np.array(freq_vector) / 1e6)

                if found_time_end_bound:
                    break

        print(f"# of sets: {set_count}")

        # Check if we have data to plot
        if set_count == 0:
            print("No valid data found for plotting.")
            exit(1)

        # Reformat data for simple plotting
        print("Reformatting data for plotting...")
        power   = np.arange(set_count * num_samples).reshape(num_samples, set_count)
        for i in range(0,num_samples):
            for j in range(0, set_count):
                power[i][j] = data[2][j][i]

        time    = np.array(data[0])
        freq    = np.array(data[1][0])

        # Debug: Print data shapes to verify if valid plotting data
        # print("shape: ", time.shape)
        # print("shape: ", freq.shape)
        # print("shape: ", power.shape)
        # print("Checking for NaN in Power...")
        # print("NaNs in power: ", np.isnan(power).any())

        # Plot spectrum data
        plt.figure(figsize=(16, 14))
        print("Creating Temporal Spectrum Analysis Plot...")
        plt.pcolormesh(
            time,
            freq,
            power,
            shading='auto',
            cmap='viridis',
            vmin=0,
            vmax=50e3
        )
        cbar = plt.colorbar(label='Power (N/A)')

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
        plot_name = f"spectrogram.{args.stid}.{ch}.{date_str}"
        plot_name += f"_{args.time_start.strftime('%H%M%S')}-{args.time_end.strftime('%H%M%S')}"
        plot_name += ".png"
        save_plot_as = os.path.join(PLOT_DIRECTORY_PATH, plot_name)
        plt.savefig(save_plot_as)


        print("Plotting Finished")
        print("Saved as:", save_plot_as)