function res = test_linearSys_reach_08_adaptive
% test_linearSys_reach_08_adaptive - unit test for adaptive reachability 
%   analysis; reduced version of testLong_linearSys_reach_08
%
% Syntax:
%    res = test_linearSys_reach_08_adaptive
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

% Authors:       Mark Wetzlinger
% Written:       08-October-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% System Dynamics----------------------------------------------------------

A = [-0.1 -2; 2 -0.1];
dim_x = length(A);
B = 1;
sys = linearSys('sys',A,B);


% Parameters and Reachability Settings ------------------------------------

params.tFinal = log(0.5)/real(eigs(A,1,'smallestabs'));
options.linAlg = 'adaptive';
options.error = 0.1;


% Reachability Analysis ---------------------------------------------------

% 1. homogeneous system
params.R0 = zonotope(10*ones(dim_x,1),0.5*eye(dim_x));
params.U = zonotope(zeros(dim_x,1));
R = reach(sys,params,options);


% 2. inhomogeneous system
params.U = zonotope(ones(dim_x,1),0.25*eye(dim_x));
R = reach(sys,params,options);


% 3. system with output matrix
C = [2 -1];
sys = linearSys('sys',A,B,[],C);
R = reach(sys,params,options);


% 4. with specification
sys = linearSys('sys',A,B);
spec = specification(polytope([-1,0],-15));
R = reach(sys,params,options,spec);

% 5. with disturbances
sys = linearSys('sys',1,1,[],1,[],[],[],1);
params.V = zonotope(interval(-1,1));
params.R0 = zonotope(interval(-1,1));
params.U = zonotope(interval(-1,1));
% high error tolerance to account for disturbances
options.error = 2;
R = reach(sys,params,options);


% ------------------------------ END OF CODE ------------------------------
