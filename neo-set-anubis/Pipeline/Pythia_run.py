import argparse
import os
import pythia_sim
from datetime import datetime

def create_pythia_generator(config_file, lhe_output, hepmc_output, num_events):
    return pythia_sim.create_pythia_generator(config_file, lhe_output, hepmc_output, "", num_events)

def process_file(config_file, output_lhe_dir, output_hepmc_dir, num_events, suffix, include_time):
    base_name = os.path.splitext(os.path.basename(config_file))[0]
    if include_time:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        base_name = f"{timestamp}_{base_name}"
    lhe_output = os.path.join(output_lhe_dir, f"{base_name}_{suffix}.lhe")
    hepmc_output = os.path.join(output_hepmc_dir, f"{base_name}_{suffix}.hepmc")
    generator = create_pythia_generator(config_file, lhe_output, hepmc_output, num_events)
    generator.generate_events()

def ensure_directories(base_dir, sub_dirs):
    for sub_dir in sub_dirs:
        full_path = os.path.join(base_dir, sub_dir)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
    return [os.path.join(base_dir, sub_dir) for sub_dir in sub_dirs]

def main():
    parser = argparse.ArgumentParser(description="Pythia Event Generator")
    parser.add_argument("-n", "--num_events", type=int, default=2000, help="Number of events to generate")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input .cmnd file, list of files, or directory containing .cmnd files")
    parser.add_argument("-d", "--input_dir", type=str, default="db/Temp/Pythia/cmnd", help="Input directory containing .cmnd files")
    parser.add_argument("-o", "--output_dir", type=str, default="db/Temp/Pythia", help="Base output directory for generated files")
    parser.add_argument("-s", "--suffix", type=str, default="output", help="Suffix for the output files")
    parser.add_argument("-t", "--time", action='store_true', help="Include timestamp in output filenames")

    args = parser.parse_args()

    output_lhe_dir, output_hepmc_dir = ensure_directories(args.output_dir, ['lhe', 'hepmc'])

    if os.path.isdir(args.input):
        config_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".cmnd")]
    elif os.path.isfile(args.input) or '/' in args.input:
        config_files = [args.input]
    else:
        config_files = [os.path.join(args.input_dir, f) for f in args.input.split(',') if f.endswith(".cmnd")]

    for config_file in config_files:
        process_file(config_file, output_lhe_dir, output_hepmc_dir, args.num_events, args.suffix, args.time)

if __name__ == "__main__":
    main()

