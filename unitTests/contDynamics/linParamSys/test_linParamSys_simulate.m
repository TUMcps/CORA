function res = test_linParamSys_simulate
% test_linParamSys_simulate - unit test for simulate
%
% Syntax:
%    res = test_linParamSys_simulate
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

% Authors:       Laura Luetzow
% Written:       28-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define linearSys --------------------------------------------------------

states = 2;
inputs = 1;
A = [-2 0; 1.5 -3];
Aw = [0 0; 0.5 0];
A_int = intervalMatrix(A,Aw);
A_zon = matZonotope(A_int);

sys_zono = linParamSys(A_zon, 1,'varParam'); %instantiate system
sys_int = linParamSys(A_int, 1,'varParam'); %instantiate system
list = {sys_zono, sys_int};


% model parameters --------------------------------------------------------

% time horizon
params.tFinal = 1;
dt_steps = 10;
% initial set
params.R0 = zonotope(10*ones(states,1),0.5*diag(ones(states,1)));
% initial state
params.x0 = center(params.R0);
% input vector
u_m = randn(inputs,dt_steps+1);


% simulate ----------------------------------------------------------------

for j = 1:length(list)
    sys = list{j};

    % no input set, disturbance, or sensor noise
    simulate(sys,params);

    % u
    params.u = u_m;
    [t,x,ind,y] = simulate(sys,params);
    assert(size(t,2) == size(x,2))
    params = rmfield(params,'u');
end

% all checks ok
res = true;
end


% ------------------------------ END OF CODE ------------------------------
