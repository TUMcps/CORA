
# CORA Docker Image

This files describes the process to update the docker image for the CI.

The `Dockerfile` in this folder pulls the desired docker matlab version
and installs all dependencies for CORA automatically (see Sec. 1.3 in the [CORA manual](https://cora.in.tum.de/manual)).

1. Update desired Matlab version in `Dockerfile`:
    - Update `ARG MATLAB_RELEASE=<VERSION>` (e.g., 2023b)
    
2. Copy `Dockerfile` from `./cora/unitTests/ci/build` to `./cora`
3. Run within `./cora`:

       docker build -t tobiasladnertum/cora:r<VERSION> --build-arg LICENSE_SERVER=28000@mlm1.rbg.tum.de .

4. Test docker image

       docker run -v "%cd%":/code/cora -e MLM_LICENSE_FILE=28000@mlm1.rbg.tum.de --rm tobiasladnertum/cora:r<VERSION> -batch "cd /code; addpath(genpath('.')); runTestSuite"

5. Push to docker explore

       docker image push tobiasladnertum/cora:r<VERSION>

6. Update `./cora/unitTests/ci/Dockerfile` <u>and</u> `./cora/unitTests/ci/.gitlab-ci.yml` to use new image

    

