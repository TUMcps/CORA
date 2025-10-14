#!/bin/bash

# example run_instance.sh script for VNNCOMP'25 for CORA
# Arguments:
# - version string "v1", 
# - benchmark identifier string, e.g., "acasxu", 
# - path to .onnx file,
# - path to .vnnlib file,
# - path to results.csv file,
# - timeout in seconds

TOOL_NAME="CORA"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
    echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
    exit 1

fi

BENCHMARK=$2
ONNX_FILE=$3
VNNLIB_FILE=$4
RESULTS_FILE=$5
TIMEOUT=$6

echo "Running $TOOL_NAME on benchmark instance $BENCHMARK with onnx file $ONNX_FILE, vnnlib file $VNNLIB_FILE, results file $RESULTS_FILE, and timeout $TIMEOUT"

# Check GPU status.
nvidia-smi

sudo matlab -nodisplay -r "run_instance('$BENCHMARK','$ONNX_FILE','$VNNLIB_FILE','$RESULTS_FILE',$TIMEOUT,true); quit;"
