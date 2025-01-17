# screen.sh - run repeatability package within a linux screen
#
# Syntax:   
#   ./screen.sh <evalname> <gpu-device>
#
# Input:
#   <evalname> - name of evaluation, defaults to date and time
#   <gpu-device> - gpu device id; integer, 'all', or 'none' (default)
#
# -------------------------------------------------------------------------

# load settings
source ./settings.sh $1 $2

# RUN ---------------------------------------------------------------------

# start screen to run repeatability package in
echo "Running repeatability package in linux screen '$SCREEN_NAME'.."
screen -S $SCREEN_NAME bash -c "./run.sh $1 $2"

# FINISH ------------------------------------------------------------------

echo
echo $SEP_LINE
echo
echo "Repeatability package complete."
echo "View results at: $PWD/results/$EVALNAME"
echo
echo $SEP_LINE
echo

# -------------------------------------------------------------------------