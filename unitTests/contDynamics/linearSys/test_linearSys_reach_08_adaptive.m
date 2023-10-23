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
tic
R = reach(sys,params,options);
tComp = toc;
disp("Computation time (homogeneous case): " + tComp);


% 2. inhomogeneous system
params.U = zonotope(ones(dim_x,1),0.25*eye(dim_x));
tic
R = reach(sys,params,options);
tComp = toc;
disp("Computation time (inhomogeneous case): " + tComp);


% 3. system with output matrix
C = [2 -1];
sys = linearSys('sys',A,B,[],C); %initialize system
tic
R = reach(sys,params,options);
tComp = toc;
disp("Computation time (output matrix case): " + tComp);


% 4. with specification
sys = linearSys('sys',A,B); %initialize system
spec = specification(halfspace([-1,0],-15));
tic
[R,res] = reach(sys,params,options,spec);
tComp = toc;
disp("Computation time (specification): " + tComp);


end

% ------------------------------ END OF CODE ------------------------------
