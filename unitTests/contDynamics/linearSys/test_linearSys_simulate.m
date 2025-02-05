function res = test_linearSys_simulate
% test_linearSys_simulate - unit test for simulate
%
% Syntax:
%    res = test_linearSys_simulate
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
% Written:       19-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define linearSys --------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
states = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
inputs = size(B,2);

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];

% constant input: q x 1
k = [0; 0.02];

% disturbance matrix: n x r
E = [1 0.5; 0 -0.5; 1 -1; 0 1];
dists = size(E,2);

% noise matrix: q x s
F = [1; 0.5];
noises = size(F,2);

% initialize different linearSys-objects
sys_A = linearSys(A,1);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);
sys_ABcCDkE = linearSys(A,B,c,C,D,k,E);
sys_ABcCDkEF = linearSys(A,B,c,C,D,k,E,F);
list = {sys_A, sys_AB, sys_ABC, sys_ABCD, sys_ABcCDk, ...
    sys_ABcCDkE, sys_ABcCDkEF};


% model parameters --------------------------------------------------------

% time horizon
params.tFinal = 1;
dt_steps = 10;
% initial set
params.R0 = zonotope(10*ones(states,1),0.5*diag(ones(states,1)));
% initial state
params.x0 = center(params.R0);
% input vector
u_n = randn(states,dt_steps+1);
u_n_noD = randn(states,dt_steps);
u_m = randn(inputs,dt_steps+1);
u_m_noD = randn(inputs,dt_steps);
% disturbance set
W = zonotope(0.02+zeros(dists,1),0.02*diag(ones(dists,1)));
% sensor noise set
V = zonotope(-0.01+zeros(noises,1),0.01*diag(ones(noises,1)));


% simulate ----------------------------------------------------------------

for j = 1:length(list)
    sys = list{j};

    % vectors for u, w, and v
    u_sys = aux_get_u(sys,states,inputs,u_n,u_m,u_n_noD,u_m_noD);
    w_sys = aux_get_w(sys,W,dt_steps);
    v_sys = aux_get_v(sys,V,dt_steps);

    % no input set, disturbance, or sensor noise
    simulate(sys,params);

    % u, no w, no v
    params.u = u_sys;
    simulate(sys,params);
    params = rmfield(params,'u');

    % no u, w, no v
    params.w = w_sys;
    simulate(sys,params);
    params = rmfield(params,'w');

    % no u, no w, v
    params.v = v_sys;
    simulate(sys,params);
    params = rmfield(params,'v');

    % u, w, no v
    params.u = u_sys;
    params.w = w_sys;
    simulate(sys,params);
    params = rmfield(params,'u');
    params = rmfield(params,'w');

    % u, no w, v
    params.u = u_sys;
    params.v = v_sys;
    simulate(sys,params);
    params = rmfield(params,'u');
    params = rmfield(params,'v');

    % no u, w, v
    params.w = w_sys;
    params.v = v_sys;
    simulate(sys,params);
    params = rmfield(params,'w');
    params = rmfield(params,'v');

    % u, w, v
    params.u = u_sys;
    params.w = w_sys;
    params.v = v_sys;
    simulate(sys,params);
    params = rmfield(params,'u');
    params = rmfield(params,'w');
    params = rmfield(params,'v');
end


% all checks ok
res = true;

end


% Auxiliary functions -----------------------------------------------------

function u = aux_get_u(sys,n,m,u_n,u_m,u_n_noD,u_m_noD)
% return corresponding u based on system

    if sys.nrOfInputs == n
        if any(any(sys.D))
            u = u_n;
        else
            u = u_n_noD;
        end
    elseif sys.nrOfInputs == m
        if any(any(sys.D))
            u = u_m;
        else
            u = u_m_noD;
        end
    end

end

function w_sys = aux_get_w(sys,W,dt_steps)
% return corresponding w based on system

    if sys.nrOfDisturbances == 1
        w_sys = randn(1,dt_steps);
    else
        w_sys = randPoint(W,dt_steps);
    end

end

function v_sys = aux_get_v(sys,V,dt_steps)
% return corresponding v based on system

    if sys.nrOfDisturbances == 1
        v_sys = randn(1,dt_steps+1);
    else
        v_sys = randPoint(V,dt_steps+1);
    end

end

% ------------------------------ END OF CODE ------------------------------
