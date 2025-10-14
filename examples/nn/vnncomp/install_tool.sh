#!/bin/bash

# Installation script used for VNN-COMP. Ubuntu, MATLAB R2024a

TOOL_NAME="CORA"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

ip link show # get mac address (for licensing)

echo $USER # get usernme (for licensing)