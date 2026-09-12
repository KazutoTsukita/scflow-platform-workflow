#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: bash script.sh id=<id> [index1=<index1>] [index2=<index2>] Read1=<Read1> [Read2=<Read2|NULL>] [codedir=<codedir>] [filereport_dir=<filereport_dir>] [final_file_dir=<final_file_dir>]"
    echo "       -h: Show this help message"
    echo "Example: bash script.sh id=857436 index1=1 index2=2 Read1=2 Read2=3 codedir=/path/to/tools/legacy filereport_dir=/path/to/filereport final_file_dir=/path/to/ready_fastqs"
    echo "If index1 or index2 is not available, either omit the argument or provide 'NULL'."
}

# Parse arguments
for arg in "$@"; do
    case $arg in
        id=*) id=${arg#*=}
              ;;
        index1=*) index1=${arg#*=}
                  ;;
        index2=*) index2=${arg#*=}
                  ;;
        Read1=*) Read1=${arg#*=}
                 ;;
        Read2=*) Read2=${arg#*=}
                 ;;
        codedir=*) codedir=${arg#*=}
                   ;;
        filereport_dir=*) filereport_dir=${arg#*=}
                          ;;
        final_file_dir=*) final_file_dir=${arg#*=}
                          ;;
        -h) show_help
            exit 0
            ;;
        *)  echo "Invalid argument: $arg"
            show_help
            exit 1
    esac
done

# Set default values for index1 and index2 if not provided
index1=${index1:-NULL}
index2=${index2:-NULL}
Read2=${Read2:-NULL}

# Check for required arguments
if [ -z "$id" ] || [ -z "$Read1" ] || [ -z "$codedir" ] || [ -z "$filereport_dir" ] || [ -z "$final_file_dir" ]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# Define necessary variables
downloaddir="${final_file_dir}/prjna${id}"
rename_script="${codedir}/rename_SRAfiles_to_bcl2fastq.py"
csv_directory="${filereport_dir}"

# Change directory
cd "${downloaddir}"

# Execute rename_script based on the presence of index1 and index2
if [ "$index1" = "NULL" ] && [ "$index2" = "NULL" ]; then
    python3 "${rename_script}" --directory="${downloaddir}" --read1="${Read1}" --read2="${Read2}"
elif [ "$index1" != "NULL" ] && [ "$index2" = "NULL" ]; then
    python3 "${rename_script}" --directory="${downloaddir}" --index1="${index1}" --read1="${Read1}" --read2="${Read2}"
elif [ "$index1" = "NULL" ] && [ "$index2" != "NULL" ]; then
    python3 "${rename_script}" --directory="${downloaddir}" --index2="${index2}" --read1="${Read1}" --read2="${Read2}"
else
    python3 "${rename_script}" --directory="${downloaddir}" --index1="${index1}" --index2="${index2}" --read1="${Read1}" --read2="${Read2}"
fi

# Execute additional commands
python3 "${codedir}/rearrange_dir_based_on_GSM.py" --id "${id}" --csv_directory "${csv_directory}" --fastq_directory "${downloaddir}"
python3 "${codedir}/rename_SRRfiles_based_on_GSM_structures.py" --dir "${downloaddir}"
