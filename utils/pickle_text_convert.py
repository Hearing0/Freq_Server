import os
import pickle
import numpy as np

def save_to_text(pickle_file, text_file):
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    raw_samples = data['raw_samples']
    sample_meta_data = data['sample_data']

    with open(text_file, 'w') as f:
        # Save metadata
        f.write(f"number_of_samples: {sample_meta_data['number_of_samples']}\n")
        f.write(f"beam_angle: {data['beam_angle']}\n")
        f.write(f"usrp_rf_rate: {sample_meta_data['usrp_rf_rate']}\n")
        f.write(f"usrp_fcenter: {sample_meta_data['usrp_fcenter']}\n")
        f.write(f"x_spacing: {sample_meta_data['x_spacing']}\n")
        f.write("antenna_list: " + ",".join(map(str, sample_meta_data['antenna_list'])) + "\n")

        # Save raw samples
        f.write("raw_samples:\n")
        for sample_array in raw_samples:
            for sample in sample_array:
                f.write(f"{sample.real},{sample.imag}\n")
# Example usage

# TODO auto-convert all files in input folder
pickle_dir = "pickle_input"  # Directory where the pickle file is stored
output_dir = "clear_freq_input"  # Directory where the output text file will be saved

pickle_file = os.path.join(pickle_dir, 'clrfreq_dump.1.pickle')
text_file = os.path.join(output_dir, 'clrfreq_dump.1.txt')

save_to_text(pickle_file, text_file)