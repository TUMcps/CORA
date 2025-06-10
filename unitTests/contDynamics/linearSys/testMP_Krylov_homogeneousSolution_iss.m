function res = testMP_Krylov_homogeneousSolution_iss(~)
% testMP_Krylov_homogeneousSolution_iss - unit_test_function for checking
%    the Krylov method for the homogeneous solution of an initial set;
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
% Last update:   02-June-2025 (MP, fix interaction with new init_reach_Krylov)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% enable access to private Krylov functions
currDir = pwd;
cd(strcat(CORAROOT,filesep,"contDynamics",filesep,"@linearSys",filesep,"private"));

% load system matrices
load('iss.mat');
dim = length(A);

% initial set
params.R0 = zonotope(interval(-0.0001*ones(dim,1),0.0001*ones(dim,1)));


%set options --------------------------------------------------------------
options.timeStep = 0.01; %time step size for reachable set computation
params.tFinal = options.timeStep;

params.U = zonotope(0);
options.taylorTerms = 6;
options.krylovError = eps;
options.reductionTechnique = 'girard';
options.zonotopeOrder = 1;
options.krylovStep = 1;
options.linAlg = 'krylov';
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
linDyn = linearSys('iss',A);
%--------------------------------------------------------------------------

% compute overapproximation
expFactor = 2+1e-17;
Rnext = reach(linDyn, params, options);
Rnext_box = interval(Rnext.timePoint.set{end});
Rnext_box = enlarge(Rnext_box,expFactor);

% compute exact solution
R_exact = expm(A*options.timeStep)*params.R0;
R_exact_box = interval(R_exact);

% revoke access to private functions
cd(currDir);

% Is exact solution in zonotope?
assert(contains(Rnext_box,R_exact_box));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
