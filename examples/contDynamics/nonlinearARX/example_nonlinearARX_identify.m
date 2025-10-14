function res = example_nonlinearARX_identify
% example_nonlinearARX_identify - example for identification of a 
%   nonlinear ARX model from trajectory data
%
% Syntax:
%    res = example_nonlinearARX_identify
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-Conformant System
%        Identification," arXiv, 2025. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearARX/identify

% Authors:       Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Original System ---------------------------------------------------------

% system from Equation (9) in [1]
f = @(y,u) [y(1,1)^2 + y(2,1)^2 + u(3,1); ...
    y(2,1)^2 + u(6,1)];
dt = 0.1;
dim_y = 2;
dim_u = 2;
p_dim = 2;

sysOrig = nonlinearARX('NARX',f,dt,dim_y, dim_u, p_dim);


% Trajectory Data ---------------------------------------------------------

params.R0 = zonotope(zeros(4,1), 0.01*eye(4));
params.U = zonotope(zeros(2,1), 0.1*eye(2));


% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 1;
options_reach.tensorOrderOutput = 2;
options_reach.verbose = false;

% Black-box approximation options
options = options_reach;
options.id.gp_parallel = true;
options.id.gp_pop_size = 50;
options.id.gp_num_gen = 30;
options.id.gp_func_names = {'times','plus', 'square'};
options.id.gp_max_genes = 2;
options.id.gp_max_depth = 2;
options.id.gp_parallel = false;
options.id.cgp_num_gen = 5;
options.id.cgp_pop_size_base = 5;
options.id.save_res = false;
options.id.p = sysOrig.n_p;
options.params.tFinal = dt*3;

% identification and validation data
n_m_id = 100;
n_s_id = 10;
n_k_id = 4;
n_m_val = 5;
n_s_val = 10;
n_k_val = 4;
traj = createTestSuite(sysOrig, params, ...
    n_k_id, n_m_id, n_s_id);
options.id.testSuite_val = createTestSuite(sysOrig, params, ...
    n_k_val, n_m_val, n_s_val);


% System Identification ---------------------------------------------------

% identify using gp algorithm
options.idAlg = 'gp';
sys_gp = nonlinearARX.identify(traj, options);

% identify using cgp algorithm
options.idAlg = 'cgp';
options.params.testSuite = options.id.testSuite_val; % data for conformance synthesis
options.params.R0 = zonotope(zeros(4,1), eye(4));
options.params.U = zonotope(zeros(2,1), eye(2));
sys_cgp = nonlinearARX.identify(traj, options);

% Simulation --------------------------------------------------------------

traj_gp(length(traj),1) = trajectory();
traj_cgp(length(traj),1) = trajectory();
for i = 1:length(traj)

    simOpts.x0 = traj(i).x(:,1,1);
    simOpts.u = traj(i).u;
    simOpts.tFinal = dt*3;

    % simulate model identified with gp
    [t,~,~,y] = simulate(sys_gp,simOpts);
    traj_gp(i) = trajectory(simOpts.u, simOpts.x0, y, t);

    % simulate model identified with cgp
    [t,~,~,y] = simulate(sys_cgp,simOpts);
    traj_cgp(i) = trajectory(simOpts.u, simOpts.x0, y, t);
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;
plotOverTime(traj, 1, 'r','Traj','y', 'DisplayName', 'true');
plotOverTime(traj_gp, 1, 'b', 'Traj','y', 'DisplayName', 'gp');
plotOverTime(traj_cgp, 1, 'g', 'Traj','y', 'DisplayName', 'cgp');
legend;

res = true;

end

% ------------------------------ END OF CODE ------------------------------
