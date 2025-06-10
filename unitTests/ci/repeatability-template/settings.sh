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
EVAL_NAME="${1:-$DATE}" # defaults to datetime
DOCKER_NAME=$NAME-$PAPERABBREV-$EVAL_NAME
SCREEN_NAME=$DOCKER_NAME

# Matlab
MAIN_SCRIPT="main('$EVAL_NAME')"

# Matlab licensing
USE_LICENSE_SERVER="true"
LICENSE_SERVER="28000@mlm1.rbg.tum.de"

# GPU
GPU_DEVICE="device=${2:-none}"

# system
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS_NAME="Linux $(uname -r)"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS_NAME="macOS $(sw_vers -productVersion)"
elif [[ "$OSTYPE" == "cygwin" || "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    OS_NAME="Windows $(uname -r)"
else
    OS_NAME="Unknown OS"
fi

# for verbose outputs
SEP_LINE="------------------------------------------------------------------"

# Output the settings -----------------------------------------------------

echo $SEP_LINE
echo
echo "Repeatability Package:"
echo "  Paper Abbreviation:               $PAPERABBREV"
echo "  Name:                             $NAME"
echo "  Date:                             $DATE"
echo
echo "Variables:"
echo "  Evaluation Name (\$EVAL_NAME):     $EVAL_NAME"
echo "  Docker Name (\$DOCKER_NAME):       $DOCKER_NAME"
echo "  Linux Screen Name (\$SCREEN_NAME): $SCREEN_NAME (if used)"
echo "  (make available via 'source settings.sh $EVAL_NAME')"
echo
echo "Matlab:"
if [[ "$USE_LICENSE_SERVER" == "true" ]]; then
    # Option 1: License server
    echo "  License Server:                   $LICENSE_SERVER"
else
    # Option 2: License file
    echo "  License File:                     license.lic"
fi
echo "  Main Script:                      $MAIN_SCRIPT"
echo
echo "System:"
echo "  OS:                               $OS_NAME"
echo "  GPU:                              $GPU_DEVICE"
echo
echo $SEP_LINE
echo

# -------------------------------------------------------------------------
