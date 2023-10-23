function res = example_linearARMAX_simulate
% example_linearARMAX_simulate - example to compare the simulation of
%   linearSysDT, ARMAX and ARMAX with tvp
%
% Syntax:
%    res = example_linearARMAX_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

% Authors:       Laura Luetzow
% Written:       16-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% user specifications
dim_x = 4;
dim_u = 2;
dim_y = 2;

% set random number stream
rng('default')

% tolerance for equality comparison
tol_ARMAX = 1e-5; 
tol_time = 1e-10;

% time horizon
tFinal = 2; % time horizon
dt = 0.1; % or 0.05 or 0.1 % sampling time

res = false;

%% create linearSys -------------------------------------------------------

% stable system matrix: n x n
A = 0.001* randn(dim_x);
dim_x = length(A);

% input matrix: n x m
B = 0.002* randn(dim_x,dim_u);
dim_u = size(B,2);

% output matrix: q x n
C = 0.003 * randn(dim_y,dim_x);

% throughput matrix: q x m
D = 0.001 * randn(dim_y,dim_u);
dim_y = size(D,1);

% true initial state
x0 = 10*ones(dim_x,1);

% initialize linearSys object
sys_lin = linearSys(A,B,[],C,D);
sys_linDT = linearSysDT(sys_lin, dt);
sys_ARMAX = linearARMAX(sys_linDT);

% model parameters --------------------------------------------------------

params.tFinal = tFinal;
dt_steps = tFinal / dt;
% initial state
params.x0 = x0; 
% disturbance sets
W = zonotope(0.02+zeros(dim_x,1),0.02*diag(ones(dim_x,1)));
V = zonotope(-0.01+zeros(dim_y,1),0.01*diag(ones(dim_y,1)));


%% simulate ---------------------------------------------------------------

% vectors for u, w, and v
params.u(1,:) = 80*sin(8*[0:dt:tFinal]);
if dim_u > 1
    params.u(2,:) = 10*sin(5*[0:dt:tFinal]+1);
end
if dim_u > 2
    params.u(3:dim_u,:) = 20*randn(dim_u-2,dt_steps+1);
end
params.w  = randPoint(W,dt_steps);
params.v = randPoint(V,dt_steps+1);

% simulate the linear discrete-time system
[t_linDT,~,~,y_linDT] = simulate(sys_linDT,params);

% simulate the linear ARMAX system with normal parametrization
params.u = [params.u; params.w zeros(dim_x,1); params.v];
start_idx = sys_ARMAX.dim +1;
params.y0 = y_linDT(1:start_idx-1, : )';
[t_ARMAX,~,~,y_ARMAX] = simulate(sys_ARMAX,params);

% simulate the linear ARMAX system with time-varying parametrization
setTVP(sys_ARMAX);
[t_tvpARMAX,~,~,y_tvpARMAX] = simulate(sys_ARMAX,params);

%% results

% check equality of the simulated outputs 
res_DT_ARMAX = withinTol(y_linDT,y_ARMAX,tol_ARMAX*dt_steps);
res_DT_ARMAX = res_DT_ARMAX & withinTol(t_linDT,t_ARMAX,tol_time);

res_ARMAX_tvpARMAX = withinTol(y_ARMAX,y_tvpARMAX,tol_ARMAX*dt_steps);
res_ARMAX_tvpARMAX = res_ARMAX_tvpARMAX & withinTol(t_ARMAX,t_tvpARMAX,tol_time);

% all checks ok
res = all(res_DT_ARMAX & res_ARMAX_tvpARMAX, 'all');

% plot the outputs
aux_plot(dim_y, t_linDT, y_linDT, t_ARMAX, y_ARMAX, t_tvpARMAX, y_tvpARMAX);

end


% Auxiliary functions -----------------------------------------------------

function aux_plot(dim_y, t_linDT, y_linDT, t_ARMAX, y_ARMAX, t_tvpARMAX, y_tvpARMAX)
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
    plot(t_ARMAX,y_ARMAX(:,i_y), ':*','Color', col(3,:), ...
        'DisplayName', 'ARMAX', 'lineWidth', 1.5);
    plot(t_tvpARMAX,y_tvpARMAX(:,i_y), '--x','Color', col(2,:),...
        'DisplayName', 'ARMAXtvp', 'lineWidth', 1.5);
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
