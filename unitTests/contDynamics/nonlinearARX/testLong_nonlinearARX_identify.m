function res = testLong_nonlinearARX_identify
% testLong_nonlinearARX_identify - unit test for system identification
%
% Syntax:
%    res = testLong_nonlinearARX_identify
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
% See also: nonlinearARX/identify

% Authors:       Laura Luetzow
% Written:       04-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define algorithms to test
alg = {'gp','cgp'};

% load dynamics of the test system
[sys, params.R0, params.U] = loadDynamics('Square');
dt = sys.dt;

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
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
options.id.p = sys.n_p;
options.params.tFinal = dt*3;

% identification and validation data
n_m_id = 100;
n_s_id = 10;
n_m_val = 5;
n_s_val = 10;
n_k = 4;
params.tFinal = sys.dt * n_k - sys.dt;
traj = createTestSuite(sys, params, ...
    n_k, n_m_id, n_s_id);
options.id.testSuite_val = createTestSuite(sys, params, ...
    n_k, n_m_val, n_s_val);

% identify using gp algorithm
options.idAlg = 'gp';
sys_gp = nonlinearARX.identify(traj, options);
assert(isa(sys_gp,'nonlinearARX'));

% identify using cgp algorithm
options.idAlg = 'cgp';
options.params.testSuite = options.id.testSuite_val; % data for conformance synthesis
options.params.R0 = zonotope(zeros(sys.n_p*sys.nrOfOutputs,1), eye(sys.n_p*sys.nrOfOutputs));
options.params.U = zonotope(zeros(sys.nrOfInputs,1), eye(sys.nrOfInputs));
sys_cgp = nonlinearARX.identify(traj, options);
assert(isa(sys_cgp,'nonlinearARX'));

% test successfull if it runs through without errors
res = true;

% ------------------------------ END OF CODE ------------------------------
