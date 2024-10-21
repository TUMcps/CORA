function res = testLong_linearSysDT_conform_03_gray
% testLong_linearSysDT_conform_03_gray - example script for reachset 
%   conformance synthesis using a gray-box model. This method identifies 
%   parameters of a gray-box model using fmincon in an outer loop as in 
%   Sec. III.B of [1,2].
%
% Syntax:
%    res = testLong_linearSysDT_conform_03_gray
%
% Inputs:
%    None
%
% Outputs:
%    res - boolean
%
% References:
%    [1] S. B. Liu and M. Althoff, "Online Verification of 
%        Impact-Force-Limiting Control for Physical Human-Robot 
%        Interaction," 2021 IEEE/RSJ International Conference on 
%        Intelligent Robots and Systems (IROS), 2021, pp. 777-783.
%    [2] L. Luetzow and M. Althoff, "Reachset-conformant System 
%        Identification," Transactions on Automatic Control, 2024. (?) 
%

% Authors:       Laura Luetzow, Matthias Althoff
% Written:       14-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')
methods_gray = ["graySim","graySeq", "grayLS"];

% define state space system 
dt = 0.1;

% system matrices
A = [0.5, -0.1;
      0.03, 0.4];
B = [ -0.01, -0.03;
      -0.02, 0.03];
C = [ 24, -110;
     150,  50];
D = [ 0, -2;
      2,  0];
  
% instantiate system  
sys = linearSysDT(linearSys(A,B,[],C,D), dt);

% values for conformance
n_m = 20; 
n_s = 10;
n_k = 6; % number of prediction time steps for one test case

% reachability and conformance options
options = struct();
options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 3;
options.tensorOrderOutput = 3;
options.errorOrder = 3;
options.cs.cost = 'interval';
options.cs.constraints = 'half';
options.cs.timeout = 300; %timeout for fmincon when estimating c
options.cs.set_p = @(p,params) aux_set_p(p, params, [2 2 2], dt);
options.cs.p0 = 0.01*randn(16);

% options for creating the test data
options_testS.p_extr = 0.4;
options_testS.p_biased = 0;


% Parameters ----------------------------------------------------------
params_true.tStart = 0;
params_true.tFinal = dt * (n_k-1);
params_true.U = zonotope([zeros(2,1),diag([1 0.9])]);
params_true.R0 = zonotope([randn(2,1),diag([0. 0.12])]);

% uncertainty template for gray-box identification
params_gray = params_true;
params_gray.R0 = zonotope([zeros(2,1),eye(2)]);
params_gray.U = zonotope([zeros(2,1),eye(2)]);

% test case options for creating the data
options_testS.stateSet = interval(-10 * ones(2,1), 10 * ones(2,1));
params_gray.testSuite = createTestSuite(sys, params_true, n_k, n_m, ...
        n_s, options_testS);

% approximate the dynamics using different methods
options_isconform = rmfield(options, ["tensorOrder", "tensorOrderOutput", "errorOrder", "cs"]);
for i = 1:length(methods_gray)
    [params_conform, ~] = conform(sys,params_gray,options,methods_gray(i));

    % Check conformance of obtained model
    confAlg = 'BF';
    params_conform.U = enlarge(params_conform.U,1.001); % slightly enlarge V for numerical robustness
    assert(isconform(sys,params_conform,options_isconform,confAlg));
end

% example completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function [sys,params] = aux_set_p(p, params, dims, dt)
[A, B, C, D] = getMatricesFromP_twoDimExample(p,dims);
sys = linearSysDT(A,B,[],C,D,[],dt);
end


% ------------------------------ END OF CODE ------------------------------
