function res = example_linearARX_simulate
% example_linearARX_simulate - example to compare the simulation of
%   linearSysDT, ARX and ARX with tvp
%
% Syntax:
%    res = example_linearARX_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in 
%       Proc. of the 62nd IEEE Conference on Decision and Control, pp. 
%       7021–7028, 2023.

% Authors:       Laura Luetzow
% Written:       16-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% user specifications
nrStates = 4;
nrInputs = 2;
nrOutputs = 2;
nrDists = nrStates;
nrNoises = nrOutputs;

% tolerance for equality comparison
tol_ARX = 1e-5; 
tol_time = 1e-10;

% time horizon
tFinal = 2; % time horizon
dt = 0.1; % or 0.05 or 0.1 % sampling time

%% create linearSys -------------------------------------------------------

% system matrix: n x n
A = 0.001* randn(nrStates);

% input matrix: n x m
B = 0.002* randn(nrStates,nrInputs);

% output matrix: q x n
C = 0.003 * randn(nrOutputs,nrStates);

% throughput matrix: q x m
D = 0.001 * randn(nrOutputs,nrInputs);

% disturbance matrix: n x r
E = eye(nrDists);

% noise matrix: q x s
F = eye(nrNoises);

% initialize linearSys object
sys_linDT = linearSysDT(A,B,[],C,D,[],E,F,dt);
sys_ARX = linearARX(sys_linDT);

% model parameters --------------------------------------------------------

params.tFinal = tFinal;
dt_steps = tFinal / dt;
% initial state
params.x0 = 10*ones(nrStates,1);
% disturbance sets
W = zonotope(0.02+zeros(nrDists,1),0.02*diag(ones(nrDists,1)));
V = zonotope(-0.01+zeros(nrNoises,1),0.01*diag(ones(nrNoises,1)));


%% simulate ---------------------------------------------------------------

% vectors for u, w, and v
params.u(1,:) = 80*sin(8*[0:dt:tFinal]);
if nrInputs > 1
    params.u(2,:) = 10*sin(5*[0:dt:tFinal]+1);
end
if nrInputs > 2
    params.u(3:nrInputs,:) = 20*randn(nrInputs-2,dt_steps+1);
end
params.w  = randPoint(W,dt_steps);
params.v = randPoint(V,dt_steps+1);

% simulate the linear discrete-time system
[t_linDT,~,~,y_linDT] = simulate(sys_linDT,params);

% simulate the linear ARX system with normal parametrization
params.u = [params.u; params.w zeros(nrDists,1); params.v];
start_idx = sys_ARX.n_p + 1;
params.y0 = y_linDT(1:start_idx-1,:)';
[t_ARX,~,~,y_ARX] = simulate(sys_ARX,params);

% simulate the linear ARX system with time-varying parametrization
setTVP(sys_ARX);
[t_tvpARX,~,~,y_tvpARX] = simulate(sys_ARX,params);

%% results

% check equality of the simulated outputs 
res_DT_ARX = withinTol(y_linDT,y_ARX,tol_ARX*dt_steps);
res_DT_ARX = res_DT_ARX & withinTol(t_linDT,t_ARX,tol_time);

res_ARX_tvpARX = withinTol(y_ARX,y_tvpARX,tol_ARX*dt_steps);
res_ARX_tvpARX = res_ARX_tvpARX & withinTol(t_ARX,t_tvpARX,tol_time);

% all checks ok
res = all(res_DT_ARX & res_ARX_tvpARX,'all');

% plot the outputs
aux_plot(nrOutputs, t_linDT, y_linDT, t_ARX, y_ARX, t_tvpARX, y_tvpARX);

end


% Auxiliary functions -----------------------------------------------------

function aux_plot(dim_y, t_linDT, y_linDT, t_ARX, y_ARX, t_tvpARX, y_tvpARX)
% plot the simulated outputs

col = colormap("lines");
figure; 
for i_y = 1:dim_y
    subplot(dim_y,1,i_y);hold on,
    box on
    grid on
    ax = gca;
    ax.LineWidth = 1;
    plot(t_linDT,y_linDT(:,i_y), '-o','Color', col(1,:), ...
        'DisplayName', 'SS', 'lineWidth', 1.5);
    plot(t_ARX,y_ARX(:,i_y), ':*','Color', col(3,:), ...
        'DisplayName', 'ARX', 'lineWidth', 1.5);
    plot(t_tvpARX,y_tvpARX(:,i_y), '--x','Color', col(2,:),...
        'DisplayName', 'ARXtvp', 'lineWidth', 1.5);
    ylabel(sprintf('Output y_%d', i_y))
    if i_y == 1
        legend('Location', 'best')
    end
    if i_y == dim_y
        xlabel('Time [s]')
    end
    fontname(gca,"times");
end

end

% ------------------------------ END OF CODE ------------------------------
