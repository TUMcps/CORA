function res = test_levelSet_tightenDomain
% test_levelSet_tightenDomain - unit test function of tightenDomain
%
% Syntax:
%    res = test_levelSet_tightenDomain
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
% See also: tightenDomain.m

% Authors:       Maximilian Perschl
% Written:       08-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Define problem
zono = zonotope([0;0],[0 0.5 0 0.5;2 0 0.5 0.5]);

% Nonlinear equality levelSet
syms x y
eq = y*sin(x)^2 + x*y^2 - 4;
ls = levelSet(eq,[x;y],'==');

results_nonlinear{1} = tightenDomain(ls,zono,'forwardBackward');
results_nonlinear{2} = tightenDomain(ls,zono,'split');
results_nonlinear{3} = tightenDomain(ls,zono,'taylor',5);
results_nonlinear{4} = tightenDomain(ls,zono,'linear',5);
expectedSolution_nonlinear{1} = interval([-1;-3],[1;3]);
expectedSolution_nonlinear{2} = interval([-0.1086;-3],[1;3]);
expectedSolution_nonlinear{3} = interval([-1;-3],[1;3]);
expectedSolution_nonlinear{4} = interval([-1;-3],[1;3]);

for i = 1:length(expectedSolution_nonlinear) 
    if ~isequal(results_nonlinear{i},expectedSolution_nonlinear{i},1e-4)
        res = false;
    end
end


% ------------------------------ END OF CODE ------------------------------
