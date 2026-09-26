import csv
import os
import argparse
import shutil
import glob

def process_project(project_id, csv_directory, fastq_directory):
    csv_path = os.path.join(csv_directory, f'PRJNA{project_id}.csv')

    with open(csv_path, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            study_alias = row['sample_alias']
            run_accessions = row['run_accessions'].split()

            study_dir = os.path.join(fastq_directory, study_alias)
            os.makedirs(study_dir, exist_ok=True)

            for accession in run_accessions:
                for fastq_file in glob.glob(os.path.join(fastq_directory, f'{accession}_*.fastq.gz')):
                    shutil.move(fastq_file, study_dir)

def main():
    parser = argparse.ArgumentParser(description='Organize fastq.gz files into directories based on CSV data')
    parser.add_argument('--id', type=str, required=True, help='Project ID to process')
    parser.add_argument('--csv_directory', type=str, required=True, help='Directory containing CSV files')
    parser.add_argument('--fastq_directory', type=str, required=True, help='Directory containing fastq.gz files')

    args = parser.parse_args()

    process_project(args.id, args.csv_directory, args.fastq_directory)

if __name__ == "__main__":
    main()
