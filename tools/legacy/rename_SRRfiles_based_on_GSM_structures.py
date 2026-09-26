import argparse
import os
import re
import csv
import subprocess
from collections import defaultdict

def generate_and_execute_rename_commands(base_dir):
    """
    Automatically detect GSM directories under the base directory, generate commands
    to rename files in those directories by replacing 'SRR' with GSM directory name
    and assigning different labels based on SRR IDs. The rename details are saved
    to a CSV file named based on the last part of the base directory.

    :param base_dir: The base directory containing GSM directories
    :return: None
    """
    rename_details = []  # List to store details for CSV output

    # Automatically find GSM directories under the base directory
    gsm_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

    for gsm_dir in gsm_dirs:
        full_path = os.path.join(base_dir, gsm_dir)
        files = os.listdir(full_path)

        # Extract and sort SRR IDs
        srr_files = defaultdict(list)
        for file in files:
            match = re.search(r"(SRR\d+)_", file)
            if match:
                srr_id = match.group(1)
                srr_files[srr_id].append(file)

        sorted_srr_ids = sorted(srr_files.keys())

        # Generate new file names with different labels based on SRR IDs
        lane_counter = 1
        for srr_id in sorted_srr_ids:
            for file in srr_files[srr_id]:
                new_file = re.sub(r"SRR\d+", gsm_dir, file)
                new_file = re.sub(r"_L\d{3}_", f"_L{str(lane_counter).zfill(3)}_", new_file)
                original_file_path = os.path.join(full_path, file)
                new_file_path = os.path.join(full_path, new_file)

                # Execute the rename command
                subprocess.run(["mv", original_file_path, new_file_path], check=True)

                # Store the details for CSV output
                rename_details.append([gsm_dir, file, new_file])
            lane_counter += 1

    # Extract the last part of the base directory for the CSV file name
    base_dir_name = os.path.basename(os.path.normpath(base_dir))
    csv_file = os.path.join(base_dir, f"rename_details_{base_dir_name}.csv")

    # Write rename details to CSV
    with open(csv_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["GSM Directory", "Original File Name", "New File Name"])
        writer.writerows(rename_details)


def main():
    parser = argparse.ArgumentParser(description='Rename files in GSM directories based on directory names.')
    parser.add_argument('--dir', type=str, required=True, help='Base directory containing GSM directories.')
    args = parser.parse_args()

    generate_and_execute_rename_commands(args.dir)

if __name__ == "__main__":
    main()
