# CORA

The COntinuous Reachability Analyzer (CORA) is a collection of MATLAB classes for the formal verification of cyber-physical systems using reachability analysis. CORA integrates various vector and matrix set representations and operations on them as well as reachability algorithms of various dynamic system classes. The software is designed such that set representations can be exchanged without having to modify the code for reachability analysis. CORA is designed using the object oriented paradigm, such that users can safely use methods without concerning themselves with detailed information hidden inside the object. Since the toolbox is written in MATLAB, the installation and use is platform independent. From Release 2018 on, the direct import of SpaceEx models into CORA is also supported. The following points summerize the main capabilities of the CORA toolbox:


### Reachability Analysis for Continious Systems

CORA computes reachable sets for linear systems, nonlinear systems as well as for systems with constraints. Continious as well as discrete time models are supported. Uncertainty in the system inputs as well as uncertainty in the model parameters can be explicitely considered. In addition, CORA also provides capabilities for the simulation of dynamical models.


### Reachability Analysis for Hybrid Systems

The toolbox is also capable to calculate the reachable sets for hybrid systems. All implemented dynamic system classes can be used to describe the different continious flows for the discrete system states. Further, multiple different methods for the calculation of the intersections with guard sets are implemented in CORA.


### Geometric Sets

CORA has a modular design, making it possible to use the capabilities of the various set representations for other purposes besides reachability analysis. The toolbox implements vector set representation, e.g., zonotopes, Taylor models and intervals, as well as matrix set representations such as matrix zonotope and interval matrices.
