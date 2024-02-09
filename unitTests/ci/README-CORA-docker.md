
# How to run CORA in docker

The CORA CI runs the test suites in a docker container.
If you want to run any tests (or other scripts) locally inside docker, this file provides you with a guide to do that.

## Folder structure

	./<current folder>
	- Dockerfile		(description of docker image)
	- ./code
	  - ./cora  		(containing CORA)
	  - ./scripts		(containing any other scripts)
	- ./data    		(containing additional data)
	- ./results 		(results will be stored here)

We will map the folder `./code` to `/code` within docker, 
so any matlab code (including CORA) should be placed there.
Analogous for `./data` and `./results`.
Note that the folders are shared between docker and your local device,
so any changes affect both.

## Build docker

Run the following command once in the current folder to build the docker image:

	docker image build . -t <name>

This might take a while the first time you run it.

## Run docker


Additionally, the MATLAB licence server is passed, 
which only works within the TUM network (or via [VPN](https://www.it.tum.de/en/it/faq/internet-access-eduroam-vpn-wifi/internet-access-eduroam-vpn-wifi/how-can-i-configure-vpn-access/)).

Important: Replace `%cd%` in the following commands to get the current directory depending on your system:

	- Windows:    %cd%
	- Linux:      $(pwd)

The following commands should be run within `./cora/unitTests/ci`
	
### Run docker non-interactively

Run the following command:

	docker run -v "%cd%/code/":/code -v "%cd%/data/":/data -v "%cd%/results/":/results -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm <name> -batch "cd /code; addpath(genpath('.')); rand(1)"

Here, we only run the command `rand(1)` within docker. You might want to change that to your script name.
Any results you want to save can be placed in `/results` within docker, which are then accesible at `./results` on your local device.

_or_

### Run docker interactively

	docker run -v "%cd%/code/":/code -v "%cd%/data/":/data -v "%cd%/results/":/results -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm -it -p 6080:6080 <name> -vnc

Then, open http://localhost:6080 in your browser and you should see the desktop of the running docker container.
From there, you can open MATLAB as you would on your locally.

Note: If you are using an ssh connection, you need to forward the port.

	ssh -L local_port:destination_server_ip:remote_port ssh_server_hostname

# Run CORA test suite in docker

For example, to run the CORA test suite within docker, open `./cora/unitTests/ci/` in the command line and execute the following commands:

	# build image once
	docker image build . -t cora-ci
	
	# run test suite
	docker run -v "%cd%/../..":/code/cora -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm cora-ci -batch "cd /code; addpath(genpath('.')); runTestSuite"

or interactively

	docker run -v "%cd%/../..":/code/cora -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm -it -p 6080:6080 cora-ci -vnc
	
	and open http://localhost:6080/ in the browser.