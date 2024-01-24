<img src="./app/images/coraLogo_readme.svg" alt="CORA"/>

A Tool for COntinuous Reachability Analysis.

<a href='https://cora.in.tum.de' target='_blank'>cora.in.tum.de</a> - <a href='https://cora.in.tum.de/manual' target='_blank'>Manual</a>

![TUMcps - CORA](https://img.shields.io/static/v1?label=TUMcps&message=CORA&color=4596FF&logo=github&link=https://github.com/TUMcps/CORA)
![CORA version](https://img.shields.io/github/tag/TUMcps/CORA?include_prereleases=&sort=semver&color=4596FF&link=https://github.com/TUMcps/CORA/tags/)
![CORA CI:passing](https://img.shields.io/static/v1?label=CI&message=passing&color=78C57F)

<hr style="height: 1px;">

The COntinuous Reachability Analyzer (CORA) is a collection of MATLAB classes for the formal verification of cyber-physical systems using reachability analysis. CORA integrates various vector and matrix set representations and operations on them as well as reachability algorithms of various dynamic system classes. The software is designed such that set representations can be exchanged without having to modify the code for reachability analysis. CORA is designed using the object-oriented paradigm, such that users can safely use methods without concerning themselves with detailed information hidden inside the objects. Since the toolbox is written in MATLAB, the installation and use is platform independent. From Release 2018 on, the direct import of SpaceEx models into CORA is also supported. The following points summarize the main capabilities of the CORA toolbox:


### Reachability Analysis for Continuous Systems

CORA computes reachable sets for linear systems, nonlinear systems as well as for systems with constraints. Continuous as well as discrete-time models are supported. Uncertainty in the system inputs as well as uncertainty in the model parameters can be explicitly considered. In addition, CORA also provides capabilities for the simulation of dynamical models.


### Reachability Analysis for Hybrid Systems

The toolbox is also capable of computing the reachable sets for hybrid systems. All implemented dynamic system classes can be used to describe the different continuous flows for the discrete system states. Furthermore, various methods for the calculation of the intersections with guard sets are implemented in CORA.


### Geometric Sets

CORA has a modular design, making it possible to use the capabilities of the various set representations for other purposes besides reachability analysis. The toolbox implements vector set representation, e.g., zonotopes, Taylor models and intervals, as well as matrix set representations such as matrix zonotope and interval matrices.


### Installation

Please check Section 1.3 in the <a target='_blank' href="https://cora.in.tum.de/manual">CORA manual</a>.


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
| `./specification`| specification classes for verification |
| `./unitTests`| unit tests |

<hr style="height: 1px;">

<img src="./app/images/coraLogo_readme.svg"/>