function res = testLong_linearSys_simulate
% testLong_linearSys_simulate - unit test for simulate
%
% Syntax:
%    res = testLong_linearSys_simulate
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

% Authors:       Mark Wetzlinger, Laura Luetzow
% Written:       16-February-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system matrix
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
dim_x = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
dim_u = size(B,2);

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];
dim_y = size(D,1);

% initialize linearSys object
sys_lin = linearSys(A,B,[],C,D);


% Model parameters --------------------------------------------------------

% time horizon
params.tFinal = 1;

% initial set
params.R0 = zonotope(10*ones(dim_x,1),0.5*diag(ones(dim_x,1)));

% initial state
params.x0 = center(params.R0);

% input set
params.U = zonotope(zeros(dim_u,1),0.25*diag(ones(dim_u,1)));

% disturbance set
params.W = zonotope(0.02+zeros(dim_x,1),0.02*diag(ones(dim_x,1)));

% sensor noise set
params.V = zonotope(-0.01+zeros(dim_y,1),0.01*diag(ones(dim_y,1)));


% Reachability analysis ---------------------------------------------------

options.linAlg = 'adaptive';
options.error = 0.01;
R = reach(sys_lin,params,options);


% Simulation --------------------------------------------------------------

% sampling time for time-varying vectors
dt = 0.05;
dt_steps = params.tFinal / dt;

% vectors for u, w, and v
u = randPoint(params.U,dt_steps+1);
w = randPoint(params.W,dt_steps);
v = randPoint(params.V,dt_steps+1);

% u, w, v
params.u = u;
params.w = w;
params.v = v;
[t_lin,~,~,y_lin] = simulate(sys_lin,params);


% visualization of output set
% figure; hold on; box on;
% plot(R);
% plot(C*params.R0,[1,2],'k');
% plot(y_lin(:,1),y_lin(:,2));
% close;


% ensure that simulation is inside reachable set
res = true;
for k=1:length(y_lin)
    % read out reachable set at given time in simulation
    S = find(R,'time',t_lin(k));

    if isempty(S.timeInterval.set)
        % could not find a reachable set for the given time...
        continue
    end
    
    for j=1:length(S.timeInterval.set)
        if ~contains(S.timeInterval.set{j},y_lin(k,:)')
            res = false; break
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
