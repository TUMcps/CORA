# CORA VNN-COMP Scripts

This folder contains all scripts for submitting CORA to the VNN-COMP (see [Repository Link](https://github.com/kollerlukas/cora-vnncomp2025)).

## Run Locally

Run `./main_vnncomp.m` to start CORA on the bechmarks specified in the file. The benchmarks are stored in `data/`; place additional benchmarks there.

## VNN-COMP Submission

1. Create a new repository that contains the required scripts, e.g., [https://github.com/kollerlukas/cora-vnncomp2025](https://github.com/kollerlukas/cora-vnncomp2025). You might need to adapte the file structure.
2. Check if the configuration `config.yaml` is correct; in particular `ami`.
3. Create a license file for the correct mac address. This can be fixed in the submission system.
4. Ensure that a server restart is triggert after the installation (`post_install.sh`) to fix GPU drivers. Check the logs to make sure the benchmarks are run with the GPU.

<hr style="height: 1px;">

<img src="../app/images/coraLogo_readme.svg"/>