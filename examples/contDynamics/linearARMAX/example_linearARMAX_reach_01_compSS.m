function res = example_linearARMAX_reach_01_compSS
% example_linearARMAX_reach_01_compSS - example script for reach where the 
%   reachable set of ARMAX models is compared with the reachable set of a 
%   state space model
%
% Syntax:
%    res = example_linearARMAX_reach_01_compSS
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
% Written:       10-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% user specifications
num_samples = 500;
path_fig = ""; % save figure here
k_save = []; % time steps which are saved

% choose reachability algorithm:
% - 'exactAddition' [1, Proposition 3]
% - 'tvpGeneral' [1, Theorem 2]
% - 'tvpEfficient' [1, Theorem 3]
options_ARMAX.armaxAlg = 'exactAddition'; 

% specifications for plotting
dims_y = [1 2];
col = colormap("lines");
lines = ["-", "--", ":", "-."];
name_fig = "reach_ped";

% set random number stream
rng('default')

%% define linearSys -------------------------------------------------------

load([CORAROOT filesep 'unitTests' filesep 'contDynamics' filesep 'linearSysDT' filesep 'models' filesep 'pedestrianModel.mat'], "pedestrian");
A = pedestrian.A;
B = pedestrian.B;
C = pedestrian.C;
D = zeros(size(C,1), size(B,2));
dt = pedestrian.dt;
sys_lin = linearSysDT(A,B,[],C,D, dt);
dim_x = length(sys_lin.A);
dim_u = size(sys_lin.B,2);
dim_y = size(sys_lin.C,1);

% create linear ARMAX object (transformation with deadbeat observer gain 
% M must be used! Higher model order p leads to bigger reachable sets and 
% incorrect results for state estimation)
sys_ARMAX = linearARMAX(sys_lin, 'M'); 
p = sys_ARMAX.dim;

%% model parameters and reachability settings -----------------------------
% general parameters ------------------------------------------------------

% disturbance sets
U = zonotope(zeros(dim_u,1));
W = 10*zonotope(0.01+zeros(dim_x,1),0.02*diag(ones(dim_x,1)));
V = 10*zonotope(-0.005+zeros(dim_y,1),0.01*[diag(ones(dim_y,1)) ones(dim_y,1)]);

% true initial state (unknown)
x0 = 10*ones(dim_x,1);

% time horizon
dt_steps = 11;
tFinal = dt * dt_steps;
tStart = 0;

% input and distuabance vectors
u = randn(dim_u,dt_steps+1);
w = randPoint(W,dt_steps);
v = randPoint(V,dt_steps+1);

% options for order reduction
options.zonotopeOrder = 50;
options.reductionTechnique = 'pca';

%% initialization ---------------------------------------------------------

% compute the initial measurements y(1),..., y(p) -------------------------
% paramters for computing the initial measurements
params_init.x0 = x0;
params_init.u = u(:,1:p);
params_init.w = w(:,1:p-1);
params_init.v = v(:,1:p);
params_init.tFinal = (p-1)*dt;

% compute the initial measurements
[~,~,~,y_init] = simulate(sys_lin,params_init);

% estimate the initial state from the measurements y(1),..., y(p) ---------
% parameters for estimating X0
params_lin_estX.u = params_init.u;
params_lin_estX.y = y_init';
params_lin_estX.V = V;
params_lin_estX.W = W;

% compute X0
[X] = aux_estimateXfromY(sys_lin, params_lin_estX);
X0_est = X{end};


%% simulation of sample points --------------------------------------------

% parameters for linearARMAX
params_ARMAX_sim.tStart = tStart;
params_ARMAX_sim.tFinal = tFinal;
params_ARMAX_sim.y0 = y_init';

% simulate initial point
y_sim = cell(num_samples,1);
params_ARMAX_sim.u = [u; w zeros(dim_x,1); v];
[t_sim,~,~,y_sim{1}] = simulate(sys_ARMAX,params_ARMAX_sim);
setTVP(sys_ARMAX);

% simulate random points (sample new disturbances for each time step)
for i_sample=1:num_samples
    if i_sample < num_samples/2 
        % 300 points with maximum disturbances
        v = randPoint(V,dt_steps+1, 'extreme');
        w = randPoint(W,dt_steps, 'extreme');
    else 
        % random points
        v = randPoint(V,dt_steps+1);
        w = randPoint(W,dt_steps);
    end
    params_ARMAX_sim.u(dim_u+1:end,:)= [w zeros(dim_x,1); v]; 
    [~,~,~,y_sim{i_sample}] = simulate(sys_ARMAX,params_ARMAX_sim);    
end


%% reachability analysis --------------------------------------------------

% reachability with linearSysDT -------------------------------------------

% reachability parameters for linearSysDT
params_lin_reach.tStart = params_init.tFinal;
params_lin_reach.tFinal = tFinal;
params_lin_reach.R0 = X0_est;
params_lin_reach.W = W;
params_lin_reach.V = V;
params_lin_reach.u = u(:, p:end);

% compute the reachable set with SS model
R_linDT = reach(sys_lin,params_lin_reach,options);
R_linDT = [cell(p-1,1); R_linDT.timePoint.set];

% reachability with linearARMAX -------------------------------------------

% reachability parameters for linearARMAX
params_ARMAX_reach.tStart = tStart;
params_ARMAX_reach.tFinal = tFinal;
params_ARMAX_reach.y0 = y_init';
params_ARMAX_reach.U = zonotope(cartProd(U, cartProd(W,V)));
params_ARMAX_reach.u = [u; zeros(dim_x+dim_y, size(u,2))];

% reachable set with ARMAX & algorithms from [1]
R_ARMAX = reach(sys_ARMAX,params_ARMAX_reach,options_ARMAX);

% reachable set with ARMAX & Minkowski sum (dependency problem)
R_ARMAX_DP = aux_reachNaive(sys_ARMAX,params_ARMAX_reach);

%% Visualization
% plot the sample points
ax = [aux_plotMeas(t_sim, y_sim, "Reachable Output Sets", {'black', '*'}, dims_y)];

% plot the reachable sets
aux_plotSets(R_linDT, y_sim, ...
    "SS", {col(1,:), lines(1)}, dims_y, ax);
aux_plotSets(R_ARMAX.timePoint.set, y_sim, ...
    "ARMAX", {col(2,:), lines(2)}, dims_y, ax);
aux_plotSets(R_ARMAX_DP.timePoint.set, y_sim, ...
    "ARMAX-DP", {col(3,:), lines(3)}, dims_y, ax);

% save the results
if path_fig ~=""
    aux_savePlots(k_save,path_fig, name_fig, ax)
end

% example completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function X = aux_estimateXfromY(sys, params)
% Estimate the state sets X(k) from past measurements of a discrete-time 
% linear system
%   Y(k) = CA^k X(0) + sum_{i=0}^{k-1} CA^i (B u(k-1-i) + W) + Du(k) + V
%    -> X(0) = [C; CA^1; ...; CA^k]^{-1} ([Y(0); Y(1); ...; Y(k)]
%           - [Du(0); Du(1); ...; Du(k)] - [V; V; ...; V]
%           - [0; C(B u(0) + W); sum_{i=0}^{k-1} CA^i (B u(k-1-i) + W)]);

yu = params.y; % dim_y x p
dim_y = size(params.y,1);
p = size(params.y,2); % p must be minimal ARMAX model order

u = params.u;
E = {sys.C};
for i=1:p
    E{i+1} = sys.C*sys.A^i;
end

% Initialization for k=1
yu_stacked = yu(:,1)-sys.D*u(:,1);
V_stacked = params.V;
O = sys.C;

for k=2:p
    % compute yu-term
    for i_sum=0:k-2
        yu(:,k) = yu(:,k) - E{i_sum+1}*sys.B*u(:,k-1-i_sum);
    end
    yu_stacked = [yu_stacked; yu(:,k)-sys.D*u(:,k)];

    %compute v-term
    V_stacked = cartProd(V_stacked, params.V);

    % compute observability matrix
    O = [O; E{k}];
end

% compute w term (consider the dependencies)
W_stacked = 0;
zero_matrix = zeros(dim_y, size(params.W.G,2)+1);
for i_sum=1:p
    W_i_stacked = zonotope(zero_matrix);
    for k=2:p
        if i_sum >= k
            N_new = dim(W_i_stacked)+dim_y;
            W_i_stacked = projectHighDim(W_i_stacked,N_new);
        else
            W_i_new = E{k-i_sum}*params.W;
            W_i_stacked = zonotope([W_i_stacked.c; W_i_new.c], [W_i_stacked.G; W_i_new.G]);
        end
    end
    W_stacked = W_stacked + W_i_stacked;
end
yu_stacked = yu_stacked + (-1)*V_stacked + (-1)*W_stacked;

% compute X0
X{1} = pinv(O) * yu_stacked;

% predict state sets X1, X2, ...
if (isempty(params.W.G) || sum(abs(params.W.G), 'all') == 0) && ...
    (isempty(params.V.G) || sum(abs(params.V.G), 'all') == 0)
    % no disturbances
    for i = 1:p-1
        X{i+1,1} = sys.A*X{i,1} + sys.B*u(:,1);
    end
else
    % predict X1
    X1 = sys.A*X{1} + sys.B*u(:,1) + params.W;

    params.R0 = X1;
    params.uTransVec = u(:,2:end);
    params.y = params.y(:,2:end);

    for i = 1:size(params.y,2)
        params.y(:,i) = params.y(:,i) - sys.D * params.uTransVec(:,i);
    end

    D_new = zeros(size(sys.C,1), size(sys.B, 2));
    sys = linearSysDT(sys.A, sys.B, [], sys.C, D_new, sys.dt);

    for k = 1:p-1
        X{k+1,1} = sys.A^(k)*X{1};
        for i_sum=0:k-1
            X{k+1} = X{k+1} + sys.A^(i_sum)*(sys.B*u(:,k-i_sum) + params.W);
        end
    end
end

end

function R = aux_reachNaive(sys,params)
% compute the reachable set of an ARMAX model with Minkowski sum (without
% considering dependencies)

% time period and number of steps
p = sys.dim;
tVec = params.tStart:sys.dt:params.tFinal;
N = length(tVec)-p;

% initialize
Y = cell(N,1);
for i=1:p
    Y{i} = params.y0(:, i);
end

% loop over all reachability steps
for k=p:p+N-1
    Y{k+1,1} = sys.B_bar{1}*(params.u(:, k+1) + params.U);
    for i=1:p
        Y{k+1,1} = Y{k+1,1} + sys.A_bar{i} * Y{k+1-i,1} ...
            + sys.B_bar{i+1} * (params.u(:, k+1-i) + params.U);
    end
end

% construct reachable set object
timePoint.set = Y;
timePoint.time = num2cell(tVec');
R = reachSet(timePoint);
end

function ax = aux_plotMeas(t, y, name, lineSpec, dims)
% plot sample points

% define plotting specifications
col = lineSpec{1}; % color
marker = lineSpec{2}; % marker
num_row = 4; % number of subplots in vertical direction
num_col = 3; % number of subplots in horizontal direction
num_samples = length(y); % number of sample points
msize = min(max(1, 100/num_samples), 4); % marker size

% create plot for each timestep
for k = 1:size(t,1)
    % create figure
    if mod(k,num_row * num_col) == 1
        figure('units','normalized','outerposition',[0 0 1 1])
        sgtitle(name);
    end
    subplot(num_row,num_col,mod(k-1,num_row * num_col)+1)
    hold on; grid on;
    box on
    axis = gca;
    axis.LineWidth = 1;

    % plot each sample point
    for i_sample = 1:num_samples
        y_i = y{i_sample};
        if i_sample > 1
            plot(y_i(k, dims(1)),y_i(k, dims(2)),marker,'Color', ...
                col, 'MarkerSize', msize, 'HandleVisibility','off','LineWidth',1);
        else
            plot(y_i(k, dims(1)),y_i(k, dims(2)),marker,'Color', ...
                col, 'MarkerSize', msize, 'DisplayName',"Samples",'LineWidth',1);
        end
    end

    % create legend and labels
    xlabel(sprintf("y_%d", dims(1)))
    ylabel(sprintf("y_%d", dims(2)))
    lg = legend('Location', 'best', 'NumColumns', 2);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
    fontsize(gca,"scale", 1.1);
    fontname(gca,"times");
    ax{k} = gca;
    drawnow;
end
end

function aux_plotSets(R, y, name, lineSpec, dims, ax)
% plot reachable sets

% define tolerance for containment check
tol_contains = 1e-3;
tol_equal = 1e-5;

% define plotting specifications
col = lineSpec{1};
line = lineSpec{2};
num_samples = length(y);

for k=1:length(R)
    if (isnumeric(R{k}) && isempty(R{k})) || ...
        (isa(R{k},'contSet') && representsa(R{k},'emptySet'))
        continue
    end
    axes(ax{k});

    % convert the set to a zonotope
    if isa(R{k}, 'double') || isa(R{k}, 'polyZonotope')
        R{k} = zonotope(R{k});
    end

    % plot sets and check containment of the sample points
    if representsa(R{k},'point')
        % reachable set is a point

        plot(R{k}, 1:2,...
            'Color', col, 'DisplayName',name,'LineWidth',1);
        for i_sample = 1:min(num_samples,500)
            y_i = y{1};
            if ~all(withinTol(R{k}.c,y_i(k,:)',tol_equal))
                warning("Sample not equal to the computed result.");
            end
        end
        
    else
        if name == "ARMAX-DP"
            lw = 2;
            xl = xlim;
            yl = ylim;
        else
            lw = 1.5;
        end
        plot(R{k}, dims, 'LineStyle', line, 'Color', col, ...
            'DisplayName',name,'LineWidth',lw);
        if name == "ARMAX-DP" && k>3
            xlim(xl);
            ylim(yl);
        end
        for i_sample = 1:min(num_samples,500)
            y_i = y{i_sample};
            if ~contains(R{k}, y_i(k,:)', 'exact', tol_contains)
                warning("Sample not contained in reachable set!");
            end
        end
    end

    drawnow;
end
end

function aux_savePlots(k_save,path_fig, name_fig, ax)
% save the subplots

for k=k_save
    axes(ax{k});
    if k~= k_save(1)
        legend("off")
    end
    
    % optimize the ticks
    xlim_old = xlim;
    xadd = (xlim_old(2) - xlim_old(1))/20;
    xlim([xlim_old(1)-xadd, xlim_old(2)+xadd]);
    ylim_old = ylim;
    yadd = (ylim_old(2) - ylim_old(1))/10;
    ylim([ylim_old(1)-yadd, ylim_old(2)+yadd]);

    ticks = num2cell(xticks());
    if length(ticks) > 5
        ticks(1:2:length(ticks)) = {""};
    end
    xticklabels(ticks);
    xtickangle(0);
    ticks = num2cell(yticks());
    if length(ticks) > 5
        ticks(1:2:length(ticks)) = {""};
    end
    yticklabels(ticks);

    % save as pdf
    exportgraphics(ax{k},path_fig+name_fig+sprintf('_k%d.pdf',k-1),'ContentType','vector')
end 
end

% ------------------------------ END OF CODE ------------------------------
