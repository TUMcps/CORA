<img src="./app/images/coraLogo_readme.svg" alt="CORA"/>

A Tool for COntinuous Reachability Analysis.

<a href='https://cora.in.tum.de' target='_blank'>cora.in.tum.de</a> - <a href='https://cora.in.tum.de/manual' target='_blank'>Manual</a>

![TUMcps - CORA](https://img.shields.io/static/v1?label=TUMcps&message=CORA&color=4596FF&logo=github&link=https://github.com/TUMcps/CORA)
![CORA version](https://img.shields.io/github/tag/TUMcps/CORA?include_prereleases=&sort=semver&color=4596FF&link=https://github.com/TUMcps/CORA/tags/)
![CORA CI:passing](https://img.shields.io/static/v1?label=CI&message=passing&color=78C57F)

<hr style="height: 1px;">

The Continuous Reachability Analyzer (CORA) is a MATLAB-based toolbox designed for the formal verification of cyber-physical systems through reachability analysis. 
It offers a comprehensive suite of tools for modeling and analyzing various system dynamics, including linear, nonlinear, and hybrid systems. 
CORA supports both continuous and discrete-time systems, accommodating uncertainties in system inputs and parameters. 
These uncertainties are captured by a diverse range of set representations such as intervals, zonotopes, Taylor models, and polytopes. 
Additionally, CORA provides functionalities for the formal verification of neural networks as well as data-driven system identification with reachset conformance. 
Various converters are implemented to easily model a system in CORA such as the well-established SpaceEx format for dynamic systems and ONNX format for neural networks. 
CORA ensures the seamless integration of different reachability algorithms without code modifications and aims for a user-friendly experience through automatic parameter tuning, 
making it a versatile tool for researchers and engineers in the field of cyber-physical systems.

### Reachability Analysis for Continuous Systems

CORA computes reachable sets for linear systems, nonlinear systems as well as for systems with constraints. 
Continuous as well as discrete time models are supported. 
Uncertainty in the system inputs as well as uncertainty in the model parameters can be explicitly considered. 
In addition, CORA also provides capabilities for the simulation of dynamical models.

### Reachability Analysis for Hybrid Systems

The toolbox is also capable of computing the reachable sets for hybrid systems. 
All implemented dynamic system classes can be used to describe the different continuous flows for the discrete system states. 
Further, multiple different methods for the computation of the intersections with guard sets are implemented in CORA.

### Geometric Sets

CORA has a modular design, making it possible to use the capabilities of the various set representations for other purposes besides reachability analysis. 
The toolbox implements vector set representation, e.g., intervals, zonotopes, Taylor models, and polytopes, 
as well as matrix set representations such as matrix zonotopes and interval matrices.

### Neural Network Verification and Robust Training

CORA enables the formal verification of neural networks, both in open-loop and in closed-loop scenarios. 
Open-loop verification refers to the task where properties of the output set of a neural network are verified, e.g., correctly classified images given noisy input. 
In closed-loop scenarios, the neural network is used as a controller of a dynamic system and is neatly integrated in the reachability algorithms above, e.g., controlling a car while keeping a safe distance. 
Additionally, one can train verifiably robust neural networks in CORA.


### Installation

Please check Section 1.3 in the <a target='_blank' href="https://cora.in.tum.de/manual">CORA manual</a>.

The installation of all required toolboxes can be checked individually by running `test_requiredToolboxes` in MATLAB. 
To check whether the core functionality of CORA has been set up correctly, 
run the standard test suite `runTestSuite` which should take about 10 minutes.

Furthermore, if you clone CORA using git, please also 
i) install <a href="https://git-lfs.com/" target="_blank">git lfs</a> (large file storage) 
sand ii) run the `git lfs pull` in the command line to ensure all data files are downloaded correctly.
Then, iii) you can check the correct download of the data files in MATLAB using `test_gitlfs`.

**Repeatability package:**
We also provide a template for a repeatability package to run CORA within Docker in one click.
This can also be used to use CORA without installing Matlab directly on your device.
Please visit `./unitTests/ci/repeatability-template` for details.

### Folder Structure

| Folder | Description |
|---|---|
|`./app` | CORA app, graphical user interface |
| `./contDynamics`| continuous dynamics classes |
| `./contSet`| continuous set classes and operations |
| `./converter`| converter from various formats to CORA |
| `./discrDynamics`| discrete dynamics |
| `./examples`| examples demonstrating the capabilities of CORA |
| `./global`| global classes, functions, and macros |
| `./hybridDynamics`| hybrid dynamics classes |
| `./matrixSet`| matrix set classes |
| `./models`| model files |
| `./nn`| neural network verification and robust training |
| `./specification`| specification classes for verification |
| `./unitTests`| unit tests |

<hr style="height: 1px;">

<img src="./app/images/coraLogo_readme.svg"/>