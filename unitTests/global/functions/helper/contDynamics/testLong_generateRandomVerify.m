function res = testLong_generateRandomVerify
% testLong_generateRandomVerify - unit test for the random generation of
%    verification benchmarks
%
% Syntax:
%    res = testLong_generateRandomVerify
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       06-May-2024
% Last update:   10-May-2024 (TL, replaced 'Difficulty' with 'MaxError')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% dimensions
n = 2;
m = 1;
r = 2;

% difficulty
maxerror = 0.2;

% random generation - satisfiable
setRepSpec = 'interval';
nrSpecs = 1;
sat = true;
[sys,params,spec,sat] = generateRandomVerify(...
    "StateDimension",n,...
    "InputDimension",m,...
    "OutputDimension",r,...
    "NrSpecs",nrSpecs,...
    "SetRepSpec",setRepSpec,...
    "Satisfiable",sat,...
    "MaxError",maxerror);

% call verification algorithm
options.verifyAlg = 'reachavoid:zonotope';
assert(verify(sys,params,options,spec) == sat);

% random generation - unsatisfiable
setRepSpec = 'halfspace';
nrSpecs = 1;
sat = false;
[sys,params,spec,sat] = generateRandomVerify(...
    "StateDimension",n,...
    "InputDimension",m,...
    "OutputDimension",r,...
    "NrSpecs",nrSpecs,...
    "SetRepSpec",setRepSpec,...
    "Satisfiable",sat,...
    "MaxError",maxerror);

% call verification algorithm
options.verifyAlg = 'reachavoid:supportFunc';
assert(verify(sys,params,options,spec) == sat);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
