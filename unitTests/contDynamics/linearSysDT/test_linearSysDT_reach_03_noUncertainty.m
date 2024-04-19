function res = test_linearSysDT_reach_03_noUncertainty()
% test_linearSysDT_reach_03_noUncertainty - unit test function of linear 
%    reachability discrete-time analysis without uncertainties (compare the
%    results with the output trajectory from simulate)
%
% Syntax:
%    res = test_linearSysDT_reach_03_noUncertainty()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Laura Luetzow
% Written:       17-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% system matrix
dim_x = 4;
dim_u = 3;
dim_y = 3;
A = randn(dim_x);
B = randn(dim_x, dim_u);
C = randn(dim_y,dim_x);
D = randn(dim_y, dim_u);

% sampling time
dt = 0.04;
sys = linearSysDT('randSys',A,B,[],C,D,dt); %initialize system


% Parameters --------------------------------------------------------------

n_k = 20;
params.tFinal = n_k*dt; %final time
params.R0 = zonotope(zeros(dim_x,1)); %initial set without uncertainty
params.U = zonotope(zeros(dim_u,1)); %input set without uncertainty
params.u = randn(dim_u,n_k+1);


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder=200; %zonotope order


% Reachability Analysis (zonotope) ----------------------------------------

R = reach(sys, params, options);
params.x0 = center(params.R0); 
[~,~,~,y] = simulate(sys, params);

% compare simulation and reach results
res = true;
for k=1:n_k
    if sum(abs(center(R.timePoint.set{k})-y(k,:)'))>1e-6
        res = false;
    end
end

% ------------------------------ END OF CODE ------------------------------
