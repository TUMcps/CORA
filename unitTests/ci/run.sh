# run.sh - main script to run the test suites within docker
#
# Syntax:
#   ./run.sh <testsuite>
#
# Input:
#   <testsuite> - name of test suite: 'short' (default), 'long', 'nn', ...
#
# See also:
#   unitTests/runTestSuite
#
# -------------------------------------------------------------------------

# parse input
TESTSUITE="${1:-short}"

# build docker
echo "Building docker image.."
docker image build . -t cora-ci
echo ""

# run test suite
echo "Running test suite '$TESTSUITE'.."
docker run -v "$(pwd)/../..":/code/cora -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm cora-ci -batch "cd /code; addpath(genpath('.')); runTestSuite('$TESTSUITE')"

# -------------------------------------------------------------------------
