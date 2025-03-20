function res = test_linearSysDT_reach_04_withDmatrix()
% test_linearSysDT_reach_04_withDmatrix - unit test function of linear reachability 
%    discrete-time analysis with uncertain inputs and D matrix
%
% Compare the reachability of a linearSysDT object with random A,B,C, and D 
% matrices computed with reach and with naive reachability formulas
%
% Syntax:
%    res = test_linearSysDT_reach_04_withDmatrix()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Laura Luetzow
% Written:       24-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;


% System Dynamics ---------------------------------------------------------

n_x = 4;
n_u = 2;
n_y = 3;

% system matrices
A = randn(n_x,n_x);
B = randn(n_x,n_u);
C = randn(n_y,n_x);
D = randn(n_y,n_u);

% number of time steps
n_k = 10;

% sampling time
dt = 0.1;

% linearSysDT object
sys = linearSysDT('randSys',A,B,[],C,D,dt); %initialize system


% Parameters --------------------------------------------------------------

params.tFinal = (n_k-1)*dt; %final time
params.R0 = zonotope(ones(n_x,1),0.1*diag(randn(n_x))); %initial set
params.U = zonotope(zeros(n_u,1),0.1*diag(randn(n_u))); %input
params.u = randn(n_u,n_k);

% Reachability Settings ---------------------------------------------------

options.zonotopeOrder=200; %zonotope order

% Reachability Analysis ---------------------------------------------------

% compute reachability with reach function
R = reach(sys, params, options);

% compute reachbility with naive reachability algorithm
X{1} = params.R0;
Y{1} = sys.C*X{1} + sys.D*(params.u(:,1)+params.U);
for k=2:n_k
    X{k} = sys.A*X{k-1} + sys.B*(params.u(:,k-1)+params.U);
    Y{k} = sys.C*X{k} + sys.D*(params.u(:,k)+params.U);
end

% compare results
for k=1:n_k
    if ~isequal(Y{k},R.timePoint.set{k},1e-8)
        res = false;
    end
    assert(isequal(Y{k},R.timePoint.set{k},1e-8))
end

end

% ------------------------------ END OF CODE ------------------------------
