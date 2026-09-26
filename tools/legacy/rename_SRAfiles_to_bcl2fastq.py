import argparse
import glob
import os
import re

def is_present(value):
    return value and str(value).upper() != "NULL"


def files_for_suffix(directory, suffix):
    if str(suffix).upper() == "SE":
        return [
            path for path in glob.glob(f"{directory}/SRR*.fastq.gz")
            if re.match(r"^SRR\d+\.fastq\.gz$", os.path.basename(path))
        ]
    return glob.glob(f"{directory}/SRR*_{suffix}.fastq.gz")


def renamed_path(file_path, suffix, role, directory):
    base_name = os.path.basename(file_path)
    if str(suffix).upper() == "SE":
        new_file_name = base_name.replace(".fastq.gz", f"_S1_L001_{role}_001.fastq.gz")
    else:
        new_file_name = base_name.replace(f"_{suffix}.fastq.gz", f"_S1_L001_{role}_001.fastq.gz")
    return os.path.join(directory, new_file_name)


def rename_files(read1_suffix, read2_suffix=None, index1_suffix=None, index2_suffix=None, directory='.'):
    rename_history = []  # ファイル名の変更履歴を保存するリスト

    # Read1のファイル名を変更
    for file_path in files_for_suffix(directory, read1_suffix):
        new_file_path = renamed_path(file_path, read1_suffix, "R1", directory)
        os.rename(file_path, new_file_path)
        rename_history.append((file_path, new_file_path))

    # Read2のファイル名を変更
    if is_present(read2_suffix):
        for file_path in files_for_suffix(directory, read2_suffix):
            new_file_path = renamed_path(file_path, read2_suffix, "R2", directory)
            os.rename(file_path, new_file_path)
            rename_history.append((file_path, new_file_path))

    # Index1のファイル名を変更
    if is_present(index1_suffix):
        for file_path in files_for_suffix(directory, index1_suffix):
            new_file_path = renamed_path(file_path, index1_suffix, "I1", directory)
            os.rename(file_path, new_file_path)
            rename_history.append((file_path, new_file_path))

    # Index2のファイル名を変更
    if is_present(index2_suffix):
        for file_path in files_for_suffix(directory, index2_suffix):
            new_file_path = renamed_path(file_path, index2_suffix, "I2", directory)
            os.rename(file_path, new_file_path)
            rename_history.append((file_path, new_file_path))

    # 変更履歴をTSVファイルに書き出す
    output_file_name = f"rename_to_bcl_fastq_{os.path.basename(directory)}.tsv"
    output_path = os.path.join(directory, output_file_name)
    with open(output_path, 'w') as f:
        for old_name, new_name in rename_history:
            f.write(f"{old_name}\t{new_name}\n")

def main():
    parser = argparse.ArgumentParser(description="Rename FASTQ files for CellRanger.")
    parser.add_argument('--directory', required=True, help="Directory containing FASTQ files.")
    parser.add_argument('--read1', required=True, help="Suffix for Read1 files.")
    parser.add_argument('--read2', required=False, help="Suffix for Read2 files. Use NULL for single-end libraries.")
    parser.add_argument('--index1', required=False, help="Suffix for Index1 files (optional).")
    parser.add_argument('--index2', required=False, help="Suffix for Index2 files (optional).")
    args = parser.parse_args()

    rename_files(args.read1, args.read2, args.index1, args.index2, args.directory)

if __name__ == "__main__":
    main()
