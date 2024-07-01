function res = example_linear_verify_randomGeneration
% example_linear_verify_randomGeneration - example for the random
%    generation of verification benchmarks in 2D from [1, Sec. 4.1]
%
% Syntax:
%    res = example_linear_verify_randomGeneration
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] M. Wetzlinger and M. Althoff. "Randomized Generation of
%        Arbitrarily Difficult Verification Tasks for Linear Time-Invariant
%        Systems", ARCH 2024.

% Authors:       Mark Wetzlinger
% Written:       06-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% specify this directory for the loading of the saved parameters
dataDirectory = [CORAROOT filesep 'examples' filesep 'contDynamics' ...
    filesep 'linearSys' filesep 'data'];

% random generation - satisfiable
% ... provides us with sys, params, spec
load([dataDirectory filesep 'ARCH2024_verify.mat']);

% call verification algorithm
options.verifyAlg = 'reachavoid:zonotope';
res_verify = verify(sys,params,options,spec);
disp("Verifiable verification benchmark verified? " + res_verify);


% random generation - unsatisfiable
% [sys,params,spec,sat] = 0;
load([dataDirectory filesep 'ARCH2024_falsify.mat']);
% ... provides us with sys, params, spec

% call verification algorithm
options.verifyAlg = 'reachavoid:supportFunc';
res_falsify = verify(sys,params,options,spec);
disp("Falsifiable verification benchmark falsified? " + res_falsify);


% ------------------------------ END OF CODE ------------------------------
