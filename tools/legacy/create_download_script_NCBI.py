import os
import pandas as pd
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate a wget download script from a TSV file')
parser.add_argument('--filereport_read_run_tsv', required=True, help='Path to the TSV file')
parser.add_argument('--download_script_outputdir', required=True, help='Path to the output directory')
args = parser.parse_args()

# Get file name and extension
file_name, _ = os.path.splitext(os.path.basename(args.filereport_read_run_tsv))

# Set the file name for the download script
if not os.path.exists(args.download_script_outputdir):
    os.mkdir(args.download_script_outputdir)

download_script = os.path.join(args.download_script_outputdir, f"{file_name}_download_srr.sh")

# If the download script already exists, delete it
if os.path.exists(download_script):
    os.remove(download_script)

# Load TSV file
df = pd.read_csv(args.filereport_read_run_tsv, sep='\t')

# Extract URLs from the 'sra_ftp' column
urls = df['run_accession'].dropna().unique()

# Options for wget command.
# Use non-verbose output so project logs are not flooded by wget progress lines.
wget_options = "-c -t 20 --waitretry=100 --retry-connrefused --timeout=1000 --read-timeout=1000 --no-verbose"

# Generate the download script
with open(download_script, 'w') as file:
    for url in urls:
        file.write(f'wget {wget_options} "https://sra-pub-run-odp.s3.amazonaws.com/sra/{url}/{url}"\n')


print(f'Download script has been generated: {download_script}')
