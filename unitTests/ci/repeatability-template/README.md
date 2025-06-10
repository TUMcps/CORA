
> # Creating a repeatability package
> This block contains information for creating the repeatability package.
> - Update: 
>   - In `README.md` (this file): Update lines marked with ❗. All other lines are information for reviewers (and you) to run the repeatabiltiy package.
>   - In `./settings.sh`: Update settings
>     - `NAME` is used for naming the Docker container such that we know whom to blame if it requires too much memory. Use `anonymous` for blind submissions.
>     - Use `LICENSE_SERVER=28000@mlm1.rbg.tum.de` within TUM network while testing.
>   - In `./code/main.m`: Update settings (and adapt scripts).
>   - In `Dockerfile`, update Matlab version
> - Add CORA to `./code` (code runs out of the box, but you want to do this step eventually).
> - Update [Important Notes](#important-notes) block below.
> - For blind submissions:
>   - In `settings.sh`:
>     - Set `NAME=anonymous`
>     - Set `LICENSE_SERVER=<port>@<hostname>`
>   - In `./code/scripts`: Anonymize scripts.
>   - In `./code/cora`: Anonymize (main) classes and functions relevant for your paper.
> - Delete this block before submission.

# Repeatability Package

This folder contains a repeatability package for:

- Paper: `<paper-title>`❗
- Venue: `<venue-name>`❗

To run the repeatability package, please follow these steps:
- [Step 1: Installation](#step-1-installation) *(Docker or Matlab)*
- [Step 2: Run the code](#step-2-run-the-code) *(Docker or Matlab)*
- [Step 3: View results](#step-3-view-results)

## Folder structure


- `./`                      : base path                             
  - `./code`                : path to code                          
    - `./cora`              : path to [CORA](https://cora.in.tum.de/)
    - `./scripts`           : path to auxiliary scripts                 
    - `./main.m`            : **main Matlab script**
  - `./data`                : path to data                  
  - `./results/<evalname>`  : path to results (created after execution)
    - `./evaluation`        : path to store any evaluation results  
    - `./plots`             : path to plots    
    - `./results.txt`       : logs of all outputs to command window    
  - `./Dockerfile`          : Dockerfile
  - `./license.lic`         : place [license file](#option-2-license-file) here
  - `./README.md`           : read me file (this file)
  - `./run.sh`              : **main script** to run [from command line using Docker](#run-from-command-line-recommended)
  - `./settings.sh`         : settings for scripts
  - `./screen.sh`           : script to run `run.sh` [within a linux screen](#run-from-command-line-recommended)

## Step 1: Installation
This folder contains the code as well as a Docker file to run the code in one click.

However, you need to provide a Matlab license.
You can either specify 

- a) a [license server](#a-license-server-recommended) (recommended), 
- b) a [license file](#b-license-file), or 
- c) [run the repeatability package directly in Matlab](#c-install-matlab-and-required-toolboxes-not-recommended):

### a) License server (recommended)

1. Ask your Matlab administrator if a Matlab license server is available.
2. In `settings.sh`, configure the license server : `LICENSE_SERVER=<port>@<hostname>`. 

➡️ Proceed with [Step 2a: Run from command line using Docker](#a-run-from-command-line-using-docker-recommended)

### b) License file
Download a license file `license.lic` to run the code:
1. Create a Matlab license file: 
	For the Docker container to run Matlab, one has to create a new license file for the container.
	Log in with your Matlab account at https://www.mathworks.com/licensecenter/licenses/.
	Click on your license, and then navigate to
	1. "Install and Activate"
    1. "View activated computers"
	1. "Activate a Computer"
	(...may differ depending on how your licensing is set up).
2. Choose:
	- Release: `R2024b`❗
	- Operating System: `Linux`
	- Host ID: `0242AC11000a` (= Default MAC of Docker container)
	- Computer Login Name: `matlab`
	- Activation Label: `<any name>`
3. When prompted if the software is already installed, choose "Yes".
4. Download the file and place it next to `Dockerfile`.

➡️ Proceed with [Step 2a: Run from command line using Docker](#a-run-from-command-line-using-docker-recommended).

### c) Install Matlab and required toolboxes (not recommended)

Install Matlab on your system and install all required toolboxes for CORA (see Sec. 1.3 in [CORA manual](https://cora.in.tum.de/manual)). The CORA repository is already included in `./code/cora`, so you don't have to clone it.

➡️ Proceed with [Step 2b: Run from Matlab](#b-run-from-matlab).


## Step 2: Run the code

Run the code depending on the chosen option in [step 1](#step-1-installation):
- [a) Run from command line using Docker](#a-run-from-command-line-using-docker-recommended) (recommended).
- [b) Run from Matlab](#b-run-from-matlab)

---

### a) Run from command line using Docker (recommended)

Make all bash scripts executable using 

    chmod +x *.sh

(see also [bug fix: windows/linux line breaks](#known-error-messages) below).

#### i) Run the package within a linux screen (recommended)

If you are using [linux screens](https://www.howtogeek.com/662422/how-to-use-linuxs-screen-command/),
you can run the package in one click in a Docker container using

    ./screen.sh <evalname> <gpu-device>

where the argument `<evalname>` is used to name the evaluation run (defaults to datetime),
and the optional argument `<gpu-device>` is used to select the GPU (see [GPU settings](#gpu-settings) below).

Linux screens are helpful when running the script on a server to ensure it finishes correctly even if your connection is interrupted.
You can always detach from the screen using `CTRL+A+D` and reattach using 

    screen -rd $SCREEN_NAME

where `SCREEN_NAME` is as in `settings.sh` (see [variables](#variables) below) or using `screen -ls`.

➡️ Procceed with [Step 3: View results](#step-3-view-results)

#### ii) Run the package directly

You can also run the evaluation in one click in a Docker container using the `run.sh` script:

	./run.sh <evalname> <gpu-device>

where the arguments `<evalname>` and `<gpu-device>` are as above.

➡️ Procceed with [Step 3: View results](#step-3-view-results)

#### Variables

To set the variables `DOCKER_NAME` and `SCREEN_NAME` automatically, you can call

    source settings.sh <evalname>

which makes the variables available in the current terminal instance,
where `<evalname>` is again the name of the evaluation.

For example, you can then stop a run using

    docker stop $DOCKER_NAME

#### GPU Settings

For Docker to use the GPU, you have to specify the `<gpu-device>` Docker should use.
You can find your available GPUs using the command `nvidia-smi`.
Possible options are the GPU id (e.g., `0`), `all`, and `none` (default).
Read more about it here: https://docs.docker.com/desktop/features/gpu/.

Please note that this setting might not be necessary for this repeatability package.

---

### b) Run from Matlab

Alternatively, open this directory in Matlab and run:

	addpath(genpath('./code')); 
    main('<evalname>');

where the optional argument `<evalname>` is used to name the evaluation run (defaults to datetime).
The results will be stored to `./results/<evalname>`.
	
**Note:** Please ensure that all required toolboxes for CORA are installed (see [Step 1c:](#c-install-matlab-and-required-toolboxes-not-recommended) above).

➡️ Procceed with [Step 3: View results](#step-3-view-results)

---

## Step 3: View results
	
The results will be stored to `./results/<evalname>`.

If you run the repeatability package from the command line within Docker,
you can view intermediate results by copying the current `results` folder out of the Docker container using

    docker cp "$DOCKER_NAME":/results .

where `DOCKER_NAME` is as in `settings.sh` (see [variables](#variables) above) or using `Docker ps`.


## Important notes

- *add any information here, e.g., settings to reduce evaluation time.* ❗
- When running the evaluation in Docker, Docker might randomly stop if not enough memory is available.


## Known error messages

If running `run.sh`/`screen.sh` results in obscure error messages (`$'\r': command not found`), 
it might be due to different line breaks in `run.sh`/`screen.sh` using windows/linux. 
You can fix it using:

    sed -i 's/\r$//' *.sh

