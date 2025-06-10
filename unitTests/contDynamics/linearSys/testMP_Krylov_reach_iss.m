function res = testMP_Krylov_reach_iss
% testMP_Krylov_reach_iss - unit test for checking the Krylov method
%    for the solution of the first time interval using the ISS model.
%    This test requires the multiple precision toolbox.
%
% Syntax:
%    res = testMP_Krylov_reach_iss
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
% See also: none

% Authors:       Matthias Althoff
% Written:       13-November-2018
% Last update:   02-June-2025 (MP, rename + fix interaction with new init_reach_Krylov)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private Krylov functions
currDir = pwd;
cd(strcat(CORAROOT,filesep,"contDynamics",filesep,"@linearSys",filesep,"private"));

% load system matrices
load('iss.mat');
n = length(A);

%set options --------------------------------------------------------------
R0 = interval(-0.0001*ones(n,1),0.0001*ones(n,1));

options.taylorTerms = 1;
options.zonotopeOrder = 1;
options.reductionTechnique = 'girard';
options.timeStep = 0.02;

params.U = zonotope(interval([0;0.8;0.9],[0.1;1;1]));
params.R0 = zonotope(R0);
params.tFinal = options.timeStep;

options.krylovError = eps;
options.krylovStep = 10;

%specify continuous dynamics-----------------------------------------------
linDyn_std = linearSys('iss',full(A),full(B),[],C);
linDyn_krylov = linearSys('iss',A,B,[],C);

%--------------------------------------------------------------------------

% copy for Krylov techniques
options_Krylov = options;
options_Krylov.linAlg = 'krylov';

% compute overapproximation using standard methods
R_std = reach(linDyn_std, params, options);

% compute overapproximation using Krylov methods
R_krylov = reach(linDyn_krylov, params, options_Krylov);

% check whether solutions of Krylov method enclose those of the standard 
% method; to save computational time, the results of the Krylov method are
% boxed
% time point solution
R_tp_Krylov_boxed = interval(R_krylov.timePoint.set{end});
R_tp_boxed = interval(R_std.timePoint.set{end});
assert(contains(R_tp_Krylov_boxed,R_tp_boxed,'approx',1e-4));
% check the other way so the sets are close
assert(contains(R_tp_boxed,R_tp_Krylov_boxed,'approx',1e-4));


% time interval solution
R_ti_Krylov_boxed = interval(R_krylov.timeInterval.set{end});
R_ti_boxed = interval(R_std.timeInterval.set{end});
assert(contains(R_ti_Krylov_boxed,R_ti_boxed,'approx',1e-4));

% revoke access to private functions
cd(currDir);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
