function res = testMP_Krylov_homogeneousSolution_iss(~)
% testMP_Krylov_homogeneousSolution_iss - unit_test_function for checking
%    the Krylov method for the homegeneous solution of an initial set;
%    larger ISS system used.
%    This test requires the multiple precision toolbox.
%
% Syntax:
%    res = testMP_Krylov_homogeneousSolution_iss(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       23-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "initReach_Krylov"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSys','private','initReach_Krylov.m');
target = fullfile(path,'contDynamics','@linearSys','initReach_Krylov.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));

% load system matrices
load('iss.mat');
dim = length(A);

% initial set
options.R0 = zonotope(interval(-0.0001*ones(dim,1),0.0001*ones(dim,1)));


%set options --------------------------------------------------------------
options.timeStep = 0.01; %time step size for reachable set computation
options.tFinal = options.timeStep;

options.x0 = center(options.R0);
options.U = zonotope(zeros(3,2));
options.uTrans = zeros(3,1);
options.taylorTerms = 6;
options.krylovError = eps;
options.reductionTechnique = 'girard';
options.zonotopeOrder = 1;
options.krylovStep = 1;
%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i)= options.timeStep^(i)/factorial(i);    
end
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
linDyn = linearSys('iss',A,B);
%--------------------------------------------------------------------------

% compute overapproximation
expFactor = 2+1e-17;
[Rnext, options] = initReach_Krylov(linDyn, options.R0, options);
Rnext_box = interval(Rnext.tp);
Rnext_box = enlarge(Rnext_box,expFactor);

% compute exact solution
R_exact = expm(A*options.timeStep)*options.R0;
R_exact_box = interval(R_exact);

% revoke access to private function "initReach_Krylov"
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% Is exact solution in zonotope?
res = (R_exact_box <= Rnext_box);

% ------------------------------ END OF CODE ------------------------------
