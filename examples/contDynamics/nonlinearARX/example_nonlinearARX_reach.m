function completed = example_nonlinearARX_reach
% example_nonlinearARX_reach - example of nonlinear reachability 
%    analysis with conservative linearization with NARX models.
%
% Syntax:
%    completed = example_nonlinearARX_reach
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 

% Authors:       Laura Luetzow
% Written:       08-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set random number stream
rng('default')

% System Dynamics ---------------------------------------------------------

f = @(y,u) [0.5*y(1,1) + u(1,1) - cos(u(2,1)); ...
    0.4*y(3,1) + u(2,1)*cos(y(1,1)); 0.6*y(5,1) + u(4,1)*sin(y(1,1))];
dt = 0.1;
dim_y = 3;
dim_u = 2;
p = 2;
sys = nonlinearARX(f,dt,dim_y, dim_u, p);

% Parameters --------------------------------------------------------------

% time horizon
N_k = 21;
params.tStart = 0;
params.tFinal = dt * (N_k-1);

% initilization
params.Y0{1} = zonotope([[2; 4; 4] 0.1*eye(3)]);
params.Y0{2} = zonotope([[2; 10; 4] 0.1*eye(3)]);

% input
params.U = zonotope([0;0.05], 0.02*eye(2));
params.u = 0.1*randn(2,N_k+1);

% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 2;
options.errorOrder = 1;
options.lagrangeRem.simplify = 'simplify';


% Reachability Analysis ---------------------------------------------------

tic

options.armaxAlg = 'exactAddition';
R_exact = reach(sys, params, options);
options.armaxAlg = 'tvpGeneral';
R_gen = reach(sys, params, options);

tComp = toc;
disp(['computation time of reachablity analysis: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

sim_points = 100;
y_init1 = randPoint(params.Y0{1},sim_points, 'extreme');
y_init2 = randPoint(params.Y0{2},sim_points, 'extreme');
y_sim = cell(sim_points,1);
u = params.u;
for i = 1:sim_points
    params.y0 = [y_init1(:,i) y_init2(:,i)];
    params.u = u + randPoint(params.U, size(params.u, 2), 'extreme');
    [tVec,~,~,y_sim{i}] = simulate(sys, params);
end


% Visualization -----------------------------------------------------------

dims = {[1 2],[2 3],[1 3]};
col = colormap("lines");

for k = 1:length(dims)
    
    figure; hold on; box on;
    projDim = dims{k};
    useCORAcolors("CORA:contDynamics")    

    % plot simulation results     
    for i = 1:sim_points
        if i == 1
            plot(y_sim{i}(:, projDim(1)),y_sim{i}(:,projDim(2)),'Color',[.7 .7 .7],'Marker','*', 'LineWidth',0.05, 'MarkerEdgeColor','black', 'DisplayName', 'Simulations');
        else
            plot(y_sim{i}(:, projDim(1)),y_sim{i}(:,projDim(2)),'Color',[.7 .7 .7],'Marker','*', 'LineWidth',0.05, 'MarkerEdgeColor','black', 'HandleVisibility','off');
        end
    end

    % plot reachable sets
    for i = 1:length(R_exact.timePoint.set)
        if i <= 2
            col_1 = col(2,:);
            col_2 = col(2,:);
            name_1 = 'Y0, exactAdd';
            name_2 = 'Y0, general';
        else
            col_1 = col(1,:);
            col_2 = col(3,:);
            name_1 = 'Ypred, exactAdd';
            name_2 = 'Ypred, general';
        end
        if i == 1 || i == 3
            plot(R_exact.timePoint.set{i},projDim, 'Color', col_1, 'LineWidth',2, 'DisplayName', name_1);
            plot(R_gen.timePoint.set{i},projDim, 'Color', col_2, 'LineWidth',2, 'DisplayName', name_2);
        else
            plot(R_exact.timePoint.set{i},projDim, 'Color', col_1, 'LineWidth',2, 'HandleVisibility','off');
            plot(R_gen.timePoint.set{i},projDim, 'Color', col_2, 'LineWidth',2, 'HandleVisibility','off');
        end
    end

    % label plot
    title('Reachability of NARX Model')
    xlabel(['y_{',num2str(projDim(1)),'}']);
    ylabel(['y_{',num2str(projDim(2)),'}']);
    legend('Location', 'best');
end


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
