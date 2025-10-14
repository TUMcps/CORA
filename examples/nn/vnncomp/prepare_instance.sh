#!/bin/bash

# Example prepare_instance.sh script for VNNCOMP'25 for CORA.
# Arguments:
# - version string "v1", 
# - benchmark identifier string, e.g., "acasxu", 
# - path to .onnx file,
# - path to .vnnlib file

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

echo "Preparing $TOOL_NAME for benchmark instance '$BENCHMARK' with onnx file '$ONNX_FILE' and vnnlib file '$VNNLIB_FILE'"

# Check GPU status.
nvidia-smi

sudo matlab -nodisplay -r "prepare_instance('$BENCHMARK','$ONNX_FILE','$VNNLIB_FILE'); quit;"

exit 0
