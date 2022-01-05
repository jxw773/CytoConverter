#!/bin/bash

###############################################################################
#### Helper Functions ####
###############################################################################

## Usage description should match command line arguments defined below
usage() {
    echo "Usage: $(basename "$0")"
    echo "  --input => Input file"
    echo "  --output => Output file"
    echo "  --log => Log file"
}

###############################################################################
## SCRIPT_DIR: directory of current script
###############################################################################
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

###############################################################################
#### Parse Command-Line Arguments ####
###############################################################################
getopt --test > /dev/null
if [ $? -ne 4 ]; then
    echo "`getopt --test` failed in this environment."
    exit 1
fi

## Command line options should match usage description
OPTIONS=
LONGOPTIONS=help,input:,output:,log:

# -temporarily store output to be able to check for errors
# -e.g. use "--options" parameter by name to activate quoting/enhanced mode
# -pass arguments only via   -- "$@"   to separate them correctly
PARSED=$(\
    getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@"\
)
if [ $? -ne 0 ]; then
    # e.g. $? == 1
    #  then getopt has complained about wrong arguments to stdout
    usage
    exit 2
fi

# read getopt's output this way to handle the quoting right:
eval set -- "$PARSED"

## ****************************************************************************
## Handle each command line option. Lower-case variables, e.g., ${file}, only
## exist if they are set as environment variables before script execution.
## If the environment variable is not set, the Upper-case variable, e.g.,
## ${FILE}, is assigned from the command line parameter.
while true; do
    case "$1" in
        --help)
            usage
            exit 0
            ;;
        --input)
            if [ -z "${input}" ]; then
                INPUT=$2
            else
                INPUT=${input}
            fi
            shift 2
            ;;
        --output)
            if [ -z "${output}" ]; then
                OUTPUT=$2
            else
                OUTPUT=${output}
            fi
            shift 2
            ;;
        --log)
            if [ -z "${log}" ]; then
                LOG=$2
            else
                LOG=${log}
            fi
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option"
            usage
            exit 3
            ;;
    esac
done
## ****************************************************************************

## ****************************************************************************
## Log any variables passed as inputs
echo "Input: ${INPUT}"
echo "Output: ${OUTPUT}"
echo "Log: ${LOG}"
## ****************************************************************************

###############################################################################
#### Validate and Set Variables ####
###############################################################################

# INPUT
if [ -z "${INPUT}" ]; then
    echo "Input file required"
    echo
    usage
    exit 1
fi

INPUT_FULL=$(readlink -f "${INPUT}")
INPUT_DIR=$(dirname "${INPUT_FULL}")
INPUT_BASE=$(basename "${INPUT_FULL}")

# OUTPUT
if [ -z "${OUTPUT}" ]; then
    echo "Output file required"
    echo
    usage
    exit 1
fi

OUTPUT_FULL=$(readlink -f "${OUTPUT}")
OUTPUT_DIR=$(dirname "${OUTPUT_FULL}")
OUTPUT_BASE=$(basename "${OUTPUT_FULL}")

# LOG
if [ -z "${LOG}" ]; then
    echo "Log file required"
    echo
    usage
    exit 1
fi

LOG_FULL=$(readlink -f "${LOG}")
LOG_DIR=$(dirname "${LOG_FULL}")
LOG_BASE=$(basename "${LOG_FULL}")

# TMP location is in the same location as output file
TMP_FULL="${OUTPUT_DIR}/_tmp"
## ****************************************************************************

###############################################################################
#### Main execution logic ####
###############################################################################

# Create tmp directory
mkdir -p ${TMP_FULL}

# Split input into multiple pieces
split -l 200 ${INPUT_FULL} ${TMP_FULL}/cyto_input_

# Run CytoConverter on each piece
ls ${TMP_FULL}/cyto_input_* \
    | awk -F_input_ '{print $2};' \
    | xargs -n 1 -t -I {} -P 5 bash -c '${SCRIPT_DIR}/main_sub.R -i ${TMP_FULL}/cyto_input_{} -o ${TMP_FULL}/output_{}.txt -l ${TMP_FULL}/log_{}.txt'