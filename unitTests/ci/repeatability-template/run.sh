# run.sh - main script to run the repeatability package from command line
#
# Syntax:
#   ./run.sh <evalname> <gpu-device>
#
# Input:
#   <evalname> - name of evaluation, defaults to date and time
#   <gpu-device> - gpu device id; integer, 'all', or 'none' (default)
#
# -------------------------------------------------------------------------

# load settings
source ./settings.sh $1 $2

# RUN ---------------------------------------------------------------------

# brief summary for user
echo 'Preparing docker image. This might take several minutes to complete..'
sleep 2 # make output readable
echo

# build docker
docker build . -t "$DOCKER_NAME" --no-cache
echo
echo "Building docker image '$DOCKER_NAME' complete."
echo
echo $SEP_LINE
echo

# Matlab licensing
if [[ "$USE_LICENSE_SERVER" == "true" ]]; then
    # Option 1: License server
    LICENSE_SETTINGS="-e MLM_LICENSE_FILE=$LICENSE_SERVER"
else
    # Option 2: License file
    LICENSE_SETTINGS="--mac-address 02:42:ac:11:00:0a -e MLM_LICENSE_FILE=/home/matlab/Documents/MATLAB/license.lic"
fi

# GPU
if [[ "$GPU_DEVICE" == "device=none" ]]; then
    # no GPU
    GPU_FLAG=""
else
    # set GPU flag
    GPU_FLAG=--gpus="$GPU_DEVICE"
fi

# run docker
echo "Running docker with name '$DOCKER_NAME'.."
echo
docker run $LICENSE_SETTINGS $GPU_FLAG --name "$DOCKER_NAME" "$DOCKER_NAME" -batch "cd /; addpath(genpath('./code')); $MAIN_SCRIPT;"
echo
echo $SEP_LINE
echo

# copy results out of docker
echo "Copying results to $PWD/results/$EVALNAME.."
docker cp "$DOCKER_NAME":/results .
echo
echo $SEP_LINE
echo

# clean up
echo "Cleaning up docker.."
docker rm --force "$DOCKER_NAME"
docker image rm --force "$DOCKER_NAME"

# FINISH ------------------------------------------------------------------

echo
echo $SEP_LINE
echo
echo "Repeatability package for '$PAPERABBREV' complete."
echo "View results at: $PWD/results/$EVALNAME"
echo
echo $SEP_LINE
echo

# -------------------------------------------------------------------------
