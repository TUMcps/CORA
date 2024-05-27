function res = testLong_linearSys_simulateRandom_04
% testLong_linearSys_simulateRandom_04 - unit test for simulation of linearSys,
%    with different simulation types
%
% Syntax:
%    res = testLong_linearSys_simulateRandom_04
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       22-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% Model Parameters --------------------------------------------------------

tFinal = 5;
params = {};
% case: only R0
params.tFinal = tFinal;
params.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);
params.U = zonotope(interval([0.9; -0.25; -0.1], [1.1; 0.25; 0.1]));

% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
dim_x = length(A);          % n
dim_u = dim(params.U);   % m
B = randn(dim_x,dim_u);
sys = linearSys('fiveDimSys',A,B);

% Random Simulations ------------------------------------------------------

N = 10;
simOpt.points = N;

% Check for completion ----------------------------------------------------

types = {'standard','gaussian','rrt','constrained'};

for i=1:numel(types)
    simRes = simulateRandom(sys, params, simOpt);
    resvec(end+1) = all(size(simRes) == [N,1]);
end

% combine results
res = all(resvec);


end

% ------------------------------ END OF CODE ------------------------------
