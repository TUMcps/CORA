
# CORA CI 

# to run (part of) the cora ci locally, run the commands below

# build docker image ------------------------------------------------------

docker image build . -t cora-ci

% run docker image --------------------------------------------------------

% Get current directory
% - Windows:    %cd%
% - Linux:      $(pwd)

# run docker image
docker run -v "%cd%/../..":/code/cora -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm cora-ci -batch "cd /code/cora/; addpath(genpath('.')); runTestSuite('short')"

% or

# run docker image interactively
docker run -v "%cd%/../..":/code/cora -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm -it -p 6080:6080 cora-ci -vnc
# open http://localhost:6080