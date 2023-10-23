function res = testMP_Krylov_errorBounds(~)
% testMP_Krylov_errorBounds - unit_test_function for checking
%    the obtained error bounds of the Krylov method.
%    The exact error bounds are taken from Tab. 1 of [1].
%    This test requires the multiple precision toolbox.
%
% Syntax:
%    res = testMP_Krylov_errorBounds(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Saad, Y. Analysis of Some Krylov Subspace Approximations to the Matrix 
%        Exponential Operator, SIAM Journal on Numerical Analysis, 1992, 29, 209-228.

% Authors:       Matthias Althoff
% Written:       13-November-2018
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

N = 100;

% example matrix 1
i = 1:N;
diagVec = (i+1)./(N+1);
A = diag(diagVec);

% % example matrix 2
% c = 0.5;
% %init example matrix 2
% matrix_2 = [];
% for j = 1:N/2
%     % update a
%     a = (2*j - 1)/(N+1);
%     % block matrix
%     B = [a c; -c a];
%     % integrate in example matrix 2
%     matrix_2(end+1:end+2,end+1:end+2) = B;
% end

% exact errors from Table 1
example{1}.order = 3;
example{1}.error = 0.301e-1;
example{2}.order = 5;
example{2}.error = 0.937e-4;
example{3}.order = 6;
example{3}.error = 0.388e-5;
example{4}.order = 7;
example{4}.error = 0.137e-6;
example{5}.order = 8;
example{5}.error = 0.424e-8;
example{6}.order = 9;
example{6}.error = 0.119e-9;
example{7}.order = 10;
example{7}.error = 0.220e-10;

% set options -------------------------------------------------------------
options.timeStep = 1; %time step size for reachable set computation
options.tFinal = options.timeStep;
options.x0 = exp(A)\ones(N,1);
options.R0 = zonotope(options.x0);
options.U = zonotope(0*options.x0);
options.uTrans = 0*options.x0;
options.taylorTerms = 4;
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
linDyn = linearSys('KrylovTest',A,1); %initialize quadratic dynamics
%--------------------------------------------------------------------------

%init res
res = zeros(length(example),1);
for i = 1:length(example)
    % set Krylov error
    options.krylovError = example{i}.error;
    
    % compute overapproximation
    initReach_Krylov(linDyn, options.R0, options);
    
    % obtain Krylov order
    V = linDyn.krylov.state.V_c;
    KrylovOrder = length(V(1,:));
    
    % Is Krylov order larger or equal?
    res(i) = (KrylovOrder >= example{i}.order);
end

% revoke access to private function "initReach_Krylov"
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% All orders larger or equal?
res = all(res);

% ------------------------------ END OF CODE ------------------------------
