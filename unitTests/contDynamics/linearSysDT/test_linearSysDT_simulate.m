function res = test_linearSysDT_simulate
% test_linearSysDT_simulate - unit test for simulate
%
% Syntax:
%    res = test_linearSysDT_simulate
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

res = false;

% define linearSysDT ------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
n = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
m = size(B,2);

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];
y = 2;

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];

% constant input: q x 1
k = [0; 0.02];

% initialize different linearSys-objects
sys_A = linearSys(A,1);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);

% initialize different linearSysDT-objects
dt = 0.05;
sysDT_A = linearSysDT(sys_A,dt);
sysDT_AB = linearSysDT(sys_AB,dt);
sysDT_ABC = linearSysDT(sys_ABC,dt);
sysDT_ABCD = linearSysDT(sys_ABCD,dt);
sysDT_ABcCDk = linearSysDT(sys_ABcCDk,dt);
list = {sysDT_A, sysDT_AB, sysDT_ABC, sysDT_ABCD, sysDT_ABcCDk};

% model parameters --------------------------------------------------------

% time horizon
params.tFinal = 1;
dt_steps = params.tFinal / dt;
% initial set
params.R0 = zonotope(10*ones(n,1),0.5*diag(ones(n,1)));
% initial state
params.x0 = center(params.R0);
% input vector
u_n = randn(n,dt_steps+1);
u_n_noD = randn(n,dt_steps);
u_m = randn(m,dt_steps+1);
u_m_noD = randn(m,dt_steps);
% disturbance set
W = zonotope(0.02+zeros(n,1),0.02*diag(ones(n,1)));
% sensor noise set
V_n = zonotope(zeros(n,1),diag(ones(n,1)));
V_y = zonotope(-0.01+zeros(y,1),0.01*diag(ones(y,1)));


% simulate ----------------------------------------------------------------

for j = 1:length(list)
    sys = list{j};

    % vectors for u, w, and v
    u_sys = aux_get_u(sys,n,m,u_n,u_m,u_n_noD,u_m_noD);
    w_sys = aux_get_w(sys,W,dt_steps);
    v_sys = aux_get_v(sys,n,y,V_n,V_y,dt_steps);

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
    params = rmfield(params,'w');

    % u, no w, v
    params.v = v_sys;
    simulate(sys,params);
    params = rmfield(params,'u');

    % no u, w, v
    params.w = w_sys;
    simulate(sys,params);

    % u, w, v
    params.u = u_sys;
    simulate(sys,params);

    % remove all
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

function w = aux_get_w(sys,W,dt_steps)

    w = randPoint(W,dt_steps);

end

function v = aux_get_v(sys,n,y,V_n,V_y,dt_steps)

    if sys.nrOfOutputs == n
        v = randPoint(V_n,dt_steps+1);
    elseif sys.nrOfOutputs == y
        v = randPoint(V_y,dt_steps+1);
    end

end

% ------------------------------ END OF CODE ------------------------------
