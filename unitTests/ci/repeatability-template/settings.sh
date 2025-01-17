# settings.sh - settings for repeatability package (update as needed)
#
# Syntax:
#   ./settings.sh <evalname>
#
# Input:
#   <evalname> - name of evaluation, defaults to date and time
#
# -------------------------------------------------------------------------

# Note: Variables should be all-lowercase and words can be seperated by '-'

# general settings
NAME="myname"
PAPERABBREV="mypaper"

# setup
DATE=`date +%y%m%d-%H%M%S`
EVALNAME="${1:-$DATE}" # defaults to datetime
DOCKER_NAME=$NAME-$PAPERABBREV-$EVALNAME
SCREEN_NAME=$DOCKER_NAME

# Matlab
MAIN_SCRIPT="main('$EVALNAME')"

# Matlab licensing
USE_LICENSE_SERVER="true"
LICENSE_SERVER="28000@mlm1.rbg.tum.de"

# GPU
GPU_DEVICE="device=${2:-none}"

# system
OS_VERSION=`cat /proc/version`

# for verbose outputs
SEP_LINE="------------------------------------------------------------------"

# Output the settings -----------------------------------------------------

echo $SEP_LINE
echo
echo "Repeatability Package:"
echo "  Paper Abbreviation: $PAPERABBREV"
echo "  Name:               $NAME"
echo "  Date:               $DATE"
echo
echo "Variables:"
echo "  Evaluation Name:    $EVALNAME"
echo "  Docker Name:        $DOCKER_NAME"
echo "  Linux Screen Name:  $SCREEN_NAME (if used)"
echo
echo "Matlab:"
if [[ "$USE_LICENSE_SERVER" == "true" ]]; then
    # Option 1: License server
    echo "  License Server:     $LICENSE_SERVER"
else
    # Option 2: License file
    echo "  License File:       license.lic"
fi
echo "  Main Script:        $MAIN_SCRIPT"
echo
echo "System:"
echo "  OS:                 $OS_VERSION"
echo "  GPU:                $GPU_DEVICE"
echo
echo $SEP_LINE
echo

# -------------------------------------------------------------------------
