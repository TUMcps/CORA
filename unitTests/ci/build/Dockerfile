# Copyright 2023 The MathWorks, Inc.
# https://github.com/mathworks-ref-arch/matlab-dockerfile

# To specify which MATLAB release to install in the container, edit the value of the MATLAB_RELEASE argument.
# Use lower case to specify the release, for example: ARG MATLAB_RELEASE=r2023b
ARG MATLAB_RELEASE=2023b

# Specify the extra products to install into the image. These products can either be toolboxes or support packages.
ARG ADDITIONAL_PRODUCTS="Symbolic_Math_Toolbox Optimization_Toolbox Statistics_and_Machine_Learning_Toolbox Deep_Learning_Toolbox Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format Parallel_Computing_Toolbox"

# This Dockerfile builds on the Ubuntu-based mathworks/matlab image.
# To check the available matlab images, see: https://hub.docker.com/r/mathworks/matlab
FROM mathworks/matlab:r$MATLAB_RELEASE

# Declare the global argument to use at the current build stage
ARG MATLAB_RELEASE
ARG ADDITIONAL_PRODUCTS

# By default, the MATLAB container runs as user "matlab". To install mpm dependencies, switch to root.
USER root

# Install mpm dependencies
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install --no-install-recommends --yes \
    wget \
    unzip \
    ca-certificates \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Run mpm to install MathWorks products into the existing MATLAB installation directory,
# and delete the mpm installation afterwards.
# Modify it by setting the ADDITIONAL_PRODUCTS defined above,
# e.g. ADDITIONAL_PRODUCTS="Statistics_and_Machine_Learning_Toolbox Parallel_Computing_Toolbox MATLAB_Coder".
# If mpm fails to install successfully then output the logfile to the terminal, otherwise cleanup.

# Switch to user matlab, and pass in $HOME variable to mpm,
# so that mpm can set the correct root folder for the support packages.
WORKDIR /tmp
USER matlab
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && EXISTING_MATLAB_LOCATION=$(dirname $(dirname $(readlink -f $(which matlab)))) \
    && sudo HOME=${HOME} ./mpm install \
        --destination=${EXISTING_MATLAB_LOCATION} \
        --release=r${MATLAB_RELEASE} \
        --products ${ADDITIONAL_PRODUCTS} \
    || (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) \
    && sudo rm -f mpm /tmp/mathworks_root.log

# When running the container a license file can be mounted,
# or a license server can be provided as an environment variable.
# For more information, see https://hub.docker.com/r/mathworks/matlab

# Alternatively, you can provide a license server to use
# with the docker image while building the image.
# Specify the host and port of the machine that serves the network licenses 
# if you want to bind in the license info as an environment variable.
# You can also build with something like --build-arg LICENSE_SERVER=27000@MyServerName,
# in which case you should uncomment the following two lines.
# If these lines are uncommented, $LICENSE_SERVER must be a valid license
# server or browser mode will not start successfully.
ARG LICENSE_SERVER
ENV MLM_LICENSE_FILE=$LICENSE_SERVER

# The following environment variables allow MathWorks to understand how this MathWorks 
# product is being used. This information helps us make MATLAB even better. 
# Your content, and information about the content within your files, is not shared with MathWorks. 
# To opt out of this service, delete the environment variables defined in the following line.
# See the Help Make MATLAB Even Better section in the accompanying README to learn more: 
# https://github.com/mathworks-ref-arch/matlab-dockerfile#help-make-matlab-even-better
ENV MW_DDUX_FORCE_ENABLE=true MW_CONTEXT_TAGS=$MW_CONTEXT_TAGS,MATLAB:TOOLBOXES:DOCKERFILE:V1

WORKDIR /home/matlab
# Inherit ENTRYPOINT and CMD from base image.

# -------------------------------------------------------------------------
# CORA CUSTOM INSTALLATIONS

WORKDIR /tmp
USER matlab

# run installCORA non-interactively; requires write access to startup
RUN sudo chmod 777 /home/matlab/Documents/MATLAB/startup.m
RUN --mount=type=bind,target=/tmp/cora matlab -batch "cd /tmp/cora; installCORA(false,true,'/home/matlab')"

# reading ONNX networks within docker can cause exceptions
# due to some gui issue (see neuralNetwork/readONNXNetwork)
# fixing it on-the-fly requires writing permission
RUN sudo chmod 777 /home/matlab/Documents/MATLAB/SupportPackages/R$MATLAB_RELEASE/toolbox/nnet/supportpackages/onnx/+nnet/+internal/+cnn/+onnx/+fcn/ModelTranslation.m \
   && sudo chmod 777 /home/matlab/Documents/MATLAB/SupportPackages/R$MATLAB_RELEASE/toolbox/nnet/supportpackages/onnx/+nnet/+internal/+cnn/+onnx/CustomLayerManager.m

# -------------------------------------------------------------------------
